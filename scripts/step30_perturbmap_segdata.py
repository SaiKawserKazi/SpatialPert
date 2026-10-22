"""Priority 5: download Perturb-map imaging segmentation tables (S-BIAD267) and reduce them.

For each mouse-lobe and panel, streams <mouse>_<lobe>_<panel>_segData.csv into its own raw
directory, keeps centroids, cell area, per-cell marker means, procode and Gene, writes
data/interim/perturbmap/<mouse>_<lobe>_<panel>.parquet, keeps the raw CSV unchanged and records
its size and SHA-256 in data/raw/S-BIAD267_seg/MANIFEST.tsv. Resumable: lobes with both the
parquet and the raw CSV are skipped; with --raw-only, raw CSVs deleted by earlier versions are
fetched again (parquet untouched). Usage: step30_perturbmap_segdata.py [Lymphoid|Myeloid] [--raw-only].
"""

import hashlib
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/267/S-BIAD267/Files/bia_submission"
RAW = ROOT / "data/raw/S-BIAD267_seg"
OUT = ROOT / "data/interim/perturbmap"


PARTS = 12  # EBI throttles each connection to ~30 KB/s; parallel byte ranges add up


def fetch(url: str, dst: Path) -> None:
    """Download with PARTS parallel HTTP range requests, then concatenate."""
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=120) as r:  # noqa: S310
        size = int(r.headers["Content-Length"])
    step = -(-size // PARTS)

    def part(i: int) -> Path:
        out = dst.with_suffix(f".part{i}")
        lo, hi = i * step, min(size, (i + 1) * step) - 1
        have = out.stat().st_size if out.exists() else 0
        if lo + have > hi:
            return out
        req = urllib.request.Request(url, headers={"Range": f"bytes={lo + have}-{hi}"})
        with urllib.request.urlopen(req, timeout=300) as r, out.open("ab") as fh:  # noqa: S310
            if r.status != 206:  # server ignored the range: never append a whole file
                raise RuntimeError(f"range not honoured (HTTP {r.status})")
            while chunk := r.read(1 << 16):
                fh.write(chunk)
        return out

    for attempt in range(5):
        try:
            with ThreadPoolExecutor(PARTS) as ex:
                parts = list(ex.map(part, range(PARTS)))
            break
        except Exception as exc:  # noqa: BLE001 - resume the missing ranges
            print("retry", attempt, exc, flush=True)
    else:
        raise RuntimeError(f"failed {url}")
    got = sum(q.stat().st_size for q in parts)
    if got != size:
        for q in parts:
            q.unlink()
        raise RuntimeError(f"size mismatch {url}: expected {size}, got {got}")
    with dst.open("wb") as fh:
        for q in parts:
            fh.write(q.read_bytes())
            q.unlink()


def combos() -> list[str]:
    lines = (ROOT / "data/raw/S-BIAD267/p.txt").read_text().replace("\r", "").splitlines()
    return sorted({ln.split("/")[1] for ln in lines if ln.count("/") >= 2})


def main() -> None:
    panel = sys.argv[1] if len(sys.argv) > 1 else "Lymphoid"
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for ml in combos():
        dst = OUT / f"{ml}_{panel}.parquet"
        url = f"{BASE}/{ml}/{panel}/{ml}_{panel}_segData.csv"
        raw = RAW / f"{ml}_{panel}_segData.csv"
        raw_only = "--raw-only" in sys.argv
        if dst.exists() and (raw.exists() or not raw_only):
            continue
        t0 = time.time()
        if not raw.exists():
            try:
                fetch(url, raw)
            except Exception as exc:  # noqa: BLE001 - report and continue with the next lobe
                print(ml, panel, "download failed:", exc, flush=True)
                continue
            h = hashlib.sha256()
            with raw.open("rb") as fh:
                while chunk := fh.read(1 << 20):
                    h.update(chunk)
            with (RAW / "MANIFEST.tsv").open("a") as fh:
                fh.write(f"{raw.name}\t{raw.stat().st_size}\t{h.hexdigest()}\t{url}\n")
        if dst.exists():
            print(f"{ml} {panel}: raw restored, {raw.stat().st_size / 1e6:.0f} MB", flush=True)
            continue
        head = pd.read_csv(raw, nrows=0).columns
        keep = ["Centroid X px", "Centroid Y px", "Cell: Area", "procode", "delta", "Gene"]
        keep += [c for c in head if c.startswith("Cell: ") and c.endswith(" mean")]
        df = pd.read_csv(raw, usecols=keep, low_memory=False)
        df.columns = [c.replace("Cell: ", "").replace(" mean", "").replace("Centroid ", "").replace(" px", "").replace(" ", "_") for c in df.columns]
        df.insert(0, "mouse", ml.split("_")[0])
        df.insert(1, "lobe", ml.split("_")[1])
        df.to_parquet(dst, index=False)
        mb = raw.stat().st_size / 1e6
        print(f"{ml} {panel}: {len(df):,} cells, {mb:.0f} MB raw, {time.time() - t0:.0f}s", flush=True)
    print("done", panel)


if __name__ == "__main__":
    main()
