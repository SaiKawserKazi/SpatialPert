"""A4: lesion benchmark in a second tissue - Xenium mouse brain coronal subset (Zenodo 11209134,
CC-BY 4.0). Transcripts are binned into 10 µm squares (pseudo-cells; no segmentation available).

Writes reports/step22_brain.csv / .json.
"""

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402

from spatialpert.benchmark import cell_permutation_p, unit_rms
from spatialpert.domains import neighbour_counts, normalise_log1p
from spatialpert.effects import clone_patches, fit_control_model, patch_domains, patch_permutation_test, residuals
from spatialpert.pipeline import NicheModel, XeniumAnalysis
from spatialpert.smallcluster import wild_cluster_bootstrap_p, within_lesion_signflip_p

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
BIN = 10.0
N_PERM = 250
DELTA = 0.5


def load_brain() -> XeniumAnalysis:
    with gzip.open(ROOT / "data/raw/xenium_mouse_brain/Xenium_Mouse_Brain_Coronal_Sub.txt.gz", "rt") as fh:
        tx = pd.read_csv(fh, sep="\t")
    tx = tx[~tx["gene"].str.startswith(("NegControl", "BLANK", "Unassigned", "Deprecated"))]
    bx = np.floor(tx["x"].to_numpy() / BIN).astype(np.int64)
    by = np.floor(tx["y"].to_numpy() / BIN).astype(np.int64)
    key = bx * 1_000_003 + by
    ukey, binv = np.unique(key, return_inverse=True)
    genes, ginv = np.unique(tx["gene"].to_numpy(), return_inverse=True)
    counts = sparse.csr_matrix((np.ones(binv.size), (binv, ginv)), shape=(ukey.size, genes.size))
    tot = np.asarray(counts.sum(1)).ravel()
    keep = tot >= 10
    counts = counts[keep]
    xy = np.column_stack([(ukey[keep] // 1_000_003 + 0.5) * BIN, (ukey[keep] % 1_000_003 + 0.5) * BIN])
    y = normalise_log1p(counts)
    dens, _ = neighbour_counts(xy, np.ones(xy.shape[0], dtype=bool), radius=50.0)
    cells = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "label": "control", "domain": 0})
    cov = np.column_stack([np.log(np.asarray(counts.sum(1)).ravel()), dens])
    return XeniumAnalysis(cells=cells, y=y, genes=genes.tolist(), cov=cov)


def perm_p(r: np.ndarray, units: np.ndarray, niche: np.ndarray, seed: int) -> float:
    ud = patch_domains(units, niche)
    return float("nan") if np.unique(ud).size < 2 else patch_permutation_test(r, units, ud, n_perm=N_PERM, seed=seed).p_value


def main() -> None:
    a = load_brain()
    nm = NicheModel.fit(a, n_pcs=20)
    ok = nm.labels >= 0
    model = fit_control_model(a.y[ok], nm.labels[ok], a.cov[ok])
    ci = np.where(ok)[0]
    cxy = nm.coords[ci]
    tree = cKDTree(cxy)
    gsd = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    rng = np.random.default_rng(22)
    info = {"bins": int(a.y.shape[0]), "genes": len(a.genes), "niche_sizes": np.bincount(nm.labels[ok]).tolist(),
            "extent_um": [float(np.ptp(cxy[:, 0])), float(np.ptp(cxy[:, 1]))]}
    print(info, flush=True)
    rows = []
    for n_les in (1, 3, 5, 10, 20):
        for lam, reps in ((0.0, 300), (0.5, 150)):
            for rep in range(reps):
                il, les = grow_lesions(tree, cxy, 400, n_les, rng)
                idx = ci[il]
                d = nm.leave_out_labels(idx, idx)
                keep = d >= 0
                idx, d, les = idx[keep], d[keep], les[keep]
                r = residuals(model, a.y[idx], d, a.cov[idx])
                seed = int(rng.integers(1 << 31))
                if lam > 0:
                    r = r + DELTA * gsd * unit_rms(rng.normal(size=r.shape[1]))
                    tgt = int(np.unique(d)[seed % np.unique(d).size])
                    r[d == tgt] += lam * DELTA * gsd * unit_rms(rng.normal(size=r.shape[1]))
                rows.append({"lesions": n_les, "lambda": lam, "rep": rep,
                             "cell": cell_permutation_p(r, d, n_perm=N_PERM, seed=seed),
                             "patch30": perm_p(r, clone_patches(nm.coords[idx], link=30.0), d, seed),
                             "lesion": perm_p(r, les, d, seed),
                             "signflip": within_lesion_signflip_p(r, d, les, seed=seed),
                             "wcb": wild_cluster_bootstrap_p(r, d, les, n_boot=299, seed=seed)})
            print("L", n_les, "lam", lam, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(R / "step22_brain.csv", index=False)
    summ = []
    for L in (1, 3, 5, 10, 20):
        nul, alt = df[(df.lesions == L) & (df["lambda"] == 0)], df[(df.lesions == L) & (df["lambda"] == 0.5)]
        rec = {"lesions": L}
        for t in ("cell", "patch30", "lesion", "signflip", "wcb"):
            pn, pa = nul[t].dropna(), alt[t].dropna()
            crit = float(np.quantile(pn, 0.05)) if len(pn) else float("nan")
            rec[t] = {"type1": float(np.mean(pn < 0.05)) if len(pn) else None, "applicable": round(len(pn) / len(nul), 3),
                      "power_raw": float(np.mean(pa < 0.05)) if len(pa) else None,
                      "power_size_adjusted": float(np.mean(pa <= crit)) if len(pa) else None}
        summ.append(rec)
    (R / "step22_brain.json").write_text(json.dumps({"data": info, "summary": summ}, indent=2))
    print("done")


if __name__ == "__main__":
    main()
