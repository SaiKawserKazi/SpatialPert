"""Integrity check of the Perturb-map extraction (S-BIAD267).

For every lobe with a raw CSV: (1) size and SHA-256 equal MANIFEST.tsv; (2) the analysis
parquet equals the raw CSV in every kept column (same row count and order; numbers equal to
float64 precision; text columns identical including missing values).
Writes reports/step37_verify_perturbmap.json; exits 1 on any mismatch.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/S-BIAD267_seg"
OUT = ROOT / "data/interim/perturbmap"
RENAME = {"Centroid X px": "X", "Centroid Y px": "Y", "Cell: Area": "Area", "procode": "procode", "delta": "delta", "Gene": "Gene"}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    man = {}
    if (RAW / "MANIFEST.tsv").exists():
        for ln in (RAW / "MANIFEST.tsv").read_text().splitlines():
            name, size, digest, _ = ln.split("\t")
            man[name] = (int(size), digest)
    res, bad = {}, 0
    for raw in sorted(RAW.glob("*_Lymphoid_segData.csv")):
        ml = raw.name.replace("_Lymphoid_segData.csv", "")
        rec: dict[str, object] = {}
        size, digest = raw.stat().st_size, sha(raw)
        rec["manifest_ok"] = man.get(raw.name) == (size, digest)
        head = pd.read_csv(raw, nrows=0).columns
        means = [c for c in head if c.startswith("Cell: ") and c.endswith(" mean")]
        cols = {**RENAME, **{c: c.replace("Cell: ", "").replace(" mean", "") for c in means}}
        src = pd.read_csv(raw, usecols=list(cols), low_memory=False).rename(columns=cols)
        pq = pd.read_parquet(OUT / f"{ml}_Lymphoid.parquet")
        rec["rows_raw"], rec["rows_parquet"] = len(src), len(pq)
        ok = len(src) == len(pq)
        for c in cols.values():
            a, b = src[c], pq[c]
            if pd.api.types.is_numeric_dtype(a):
                same = bool(np.allclose(a.to_numpy(float), b.to_numpy(float), rtol=0, atol=0, equal_nan=True))
            else:
                same = bool((a.isna() == b.isna()).all() and (a.dropna().astype(str).to_numpy() == b.dropna().astype(str).to_numpy()).all())
            ok &= same
            if not same:
                rec.setdefault("mismatched_columns", []).append(c)  # type: ignore[union-attr]
        rec["columns_identical"] = ok
        bad += not (ok and rec["manifest_ok"])
        res[ml] = rec
        print(ml, rec, flush=True)
    out = {"lobes_checked": len(res), "failed": bad, "per_lobe": res}
    (ROOT / "reports/step37_verify_perturbmap.json").write_text(json.dumps(out, indent=2))
    print(f"verified {len(res) - bad}/{len(res)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
