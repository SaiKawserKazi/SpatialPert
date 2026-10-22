"""A6: core lesion benchmark repeated with niches from external callers.

Labels are read from data/interim/a6/<caller>_labels.csv (control tumour cells). External
callers saw all control cells, so their labels are circular for pseudo-knockouts; the k-means
row is run the same way (circular) for a like-for-like comparison.
Writes reports/step24_callers.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402

from spatialpert.benchmark import unit_rms
from spatialpert.effects import fit_control_model, patch_domains, patch_permutation_test, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import wild_cluster_bootstrap_p, within_lesion_signflip_p

ROOT = Path(__file__).resolve().parents[1]
A6 = ROOT / "data/interim/a6"
DELTA = 0.5


def main() -> None:
    a = load_xenium_analysis(ROOT)
    ids = a.cells["cell_id"].to_numpy()
    wanted = sys.argv[1:] or ["kmeans", "banksy", "spagcn", "cellcharter"]
    callers: dict[str, np.ndarray] = {"kmeans": NicheModel.fit(a).labels} if "kmeans" in wanted else {}
    for name in [w for w in wanted if w != "kmeans"]:
        f = A6 / f"{name}_labels.csv"
        if f.exists():
            m = pd.read_csv(f).set_index("cell_id")[name]
            callers[name] = np.array([int(m[c]) if c in m.index else -1 for c in ids])
    out: dict[str, object] = {}
    for name, lab in callers.items():
        ok = a.is_control & (lab >= 0)
        model = fit_control_model(a.y[ok], lab[ok], a.cov[ok])
        ci = np.where(ok)[0]
        xy = a.cells[["x", "y"]].to_numpy()[ci]
        tree = cKDTree(xy)
        gsd = residuals(model, a.y[ci[:20000]], lab[ci[:20000]], a.cov[ci[:20000]]).std(0)
        rng = np.random.default_rng(24)
        rows = []
        for n_les in (3, 5, 10, 20):
            for lam, reps in ((0.0, 200), (0.5, 100)):
                for rep in range(reps):
                    il, les = grow_lesions(tree, xy, 400, n_les, rng)
                    idx = ci[il]
                    d = lab[idx]
                    r = residuals(model, a.y[idx], d, a.cov[idx])
                    seed = int(rng.integers(1 << 31))
                    if lam > 0:
                        r = r + DELTA * gsd * unit_rms(rng.normal(size=r.shape[1]))
                        tgt = int(np.unique(d)[seed % np.unique(d).size])
                        r[d == tgt] += lam * DELTA * gsd * unit_rms(rng.normal(size=r.shape[1]))
                    ud = patch_domains(les, d)
                    rows.append({"lesions": n_les, "lambda": lam,
                                 "lesion": patch_permutation_test(r, les, ud, n_perm=250, seed=seed).p_value if np.unique(ud).size > 1 else np.nan,
                                 "signflip": within_lesion_signflip_p(r, d, les, seed=seed),
                                 "wcb": wild_cluster_bootstrap_p(r, d, les, n_boot=199, seed=seed)})
            print(name, "L", n_les, flush=True)
        df = pd.DataFrame(rows)
        summ = []
        for L in (3, 5, 10, 20):
            nul, alt = df[(df.lesions == L) & (df["lambda"] == 0)], df[(df.lesions == L) & (df["lambda"] == 0.5)]
            rec = {"lesions": L}
            for t in ("lesion", "signflip", "wcb"):
                pn, pa = nul[t].dropna(), alt[t].dropna()
                crit = float(np.quantile(pn, 0.05)) if len(pn) else float("nan")
                rec[t] = {"type1": float(np.mean(pn < 0.05)) if len(pn) else None, "applicable": round(len(pn) / max(len(nul), 1), 3),
                          "power_size_adjusted": float(np.mean(pa <= crit)) if len(pa) and len(pn) else None}
            summ.append(rec)
        out[name] = {"n_niches": int(np.unique(lab[lab >= 0]).size), "largest_share": float(np.bincount(lab[ok]).max() / ok.sum()), "summary": summ}
    suffix = "" if not sys.argv[1:] else "_" + "_".join(sys.argv[1:])
    (ROOT / f"reports/step24_callers{suffix}.json").write_text(json.dumps(out, indent=2))
    print("done", list(out))


if __name__ == "__main__":
    main()
