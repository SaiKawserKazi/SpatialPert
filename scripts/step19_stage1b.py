"""AMENDMENT-008 A1-A3: circular vs leave-out niche labels, small-cluster methods, size-adjusted power.

Lesion design, 400 cells, L in {1,3,5,10,20}; 300 null + 150 (lambda = 0.5) replicates.
Writes reports/step19_stage1b.csv (one row per draw with every p-value) and .json.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402

from spatialpert.baselines import patch_means, studentized_permutation_p
from spatialpert.benchmark import cell_permutation_p, unit_rms
from spatialpert.effects import clone_patches, fit_control_model, patch_domains, patch_permutation_test, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import cr2_two_niche_min_q, mixedlm_min_q, wild_cluster_bootstrap_p, within_lesion_signflip_p

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
N_PERM = 250
DELTA = 0.5
LS = [int(x) for x in sys.argv[1:] if x.isdigit()] or [1, 3, 5, 10, 20]
REPS = {0.0: 300, 0.5: 150}


def perm_p(r: np.ndarray, units: np.ndarray, niche: np.ndarray, seed: int) -> float:
    ud = patch_domains(units, niche)
    if np.unique(ud).size < 2:
        return float("nan")
    return patch_permutation_test(r, units, ud, n_perm=N_PERM, seed=seed).p_value


def main() -> None:
    t0 = time.time()
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    ok = nm.labels >= 0
    ctrl = a.is_control & ok
    model = fit_control_model(a.y[ctrl], nm.labels[ctrl], a.cov[ctrl])
    ci = np.where(ctrl)[0]
    cxy = nm.coords[ci]
    tree = cKDTree(cxy)
    gsd = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    seeds_map: list[list[float]] = []
    for n_les in LS:
        rng = np.random.default_rng(1900 + n_les)  # one stream per lesion level (parallel runs)
        rows = []
        for lam, reps in REPS.items():
            for rep in range(reps):
                il, les = grow_lesions(tree, cxy, 400, n_les, rng)
                idx = ci[il]
                lab = {"circ": nm.labels[idx], "loo": nm.leave_out_labels(idx, idx)}
                if lam == 0 and n_les == 20 and rep < 60:
                    seeds_map.extend(nm.coords[idx][les == k].mean(0).tolist() for k in np.unique(les))
                u = unit_rms(rng.normal(size=a.y.shape[1]))
                v = unit_rms(rng.normal(size=a.y.shape[1]))
                seed = int(rng.integers(1 << 31))
                row = {"lesions": n_les, "lambda": lam, "rep": rep, "cells_per_lesion": 400 / n_les}
                for tag, d in lab.items():
                    keep = d >= 0
                    ix, dd = idx[keep], d[keep]
                    ll = np.unique(les[keep], return_inverse=True)[1]  # contiguous after filtering
                    r = residuals(model, a.y[ix], dd, a.cov[ix])
                    if lam > 0:
                        r = r + DELTA * gsd * u
                        tgt = int(np.unique(dd)[seed % np.unique(dd).size])
                        r[dd == tgt] += lam * DELTA * gsd * v
                    xy = nm.coords[ix]
                    row[f"{tag}_cell"] = cell_permutation_p(r, dd, n_perm=N_PERM, seed=seed)
                    row[f"{tag}_patch30"] = perm_p(r, clone_patches(xy, link=30.0), dd, seed)
                    row[f"{tag}_lesion"] = perm_p(r, ll, dd, seed)
                    if tag == "loo":
                        ud = patch_domains(ll, dd)
                        mns, _, g = patch_means(r, ll, ud)
                        row["loo_stud_lesion"] = studentized_permutation_p(mns, g, n_perm=N_PERM, seed=seed)
                        row["loo_wcb"] = wild_cluster_bootstrap_p(r, dd, ll, n_boot=299, seed=seed)
                        row["loo_cr2"] = cr2_two_niche_min_q(r, dd, ll)
                        row["loo_signflip"] = within_lesion_signflip_p(r, dd, ll, seed=seed)
                        sub = rep < (40 if lam == 0 else 20)  # AMENDMENT-008 deviation: runtime
                        row["loo_mixedlm"] = mixedlm_min_q(r, dd, ll) if sub else float("nan")
                        row["loo_frac_relabelled"] = float(np.mean(lab["circ"][keep] != dd))
                rows.append(row)
            print(f"L={n_les} lam={lam} done {time.time() - t0:.0f}s", flush=True)
            pd.DataFrame(rows).to_csv(R / f"step19_part_L{n_les}.csv", index=False)  # incremental
        if seeds_map:
            pd.DataFrame(seeds_map, columns=["x", "y"]).to_csv(R / "step19_seed_centroids_L20.csv", index=False)


def summarise() -> None:
    """Merge per-level parts and apply the AMENDMENT-008 rules."""
    t0 = time.time()
    df = pd.concat([pd.read_csv(R / f"step19_part_L{L}.csv") for L in (1, 3, 5, 10, 20)], ignore_index=True)
    df.to_csv(R / "step19_stage1b.csv", index=False)
    LS_ALL = [1, 3, 5, 10, 20]

    tests = [c for c in df.columns if c.startswith(("circ_", "loo_")) and c != "loo_frac_relabelled"]
    summ = []
    for L in LS_ALL:
        nul = df[(df.lesions == L) & (df["lambda"] == 0)]
        alt = df[(df.lesions == L) & (df["lambda"] == 0.5)]
        rec = {"lesions": L, "frac_relabelled_loo": float(nul["loo_frac_relabelled"].mean())}
        for t in tests:
            pn, pa = nul[t].dropna(), alt[t].dropna()
            if len(pn) == 0:
                rec[t] = {"type1": None, "applicable_null": 0}
                continue
            crit = float(np.quantile(pn, 0.05))
            rec[t] = {
                "type1": float(np.mean(pn < 0.05)), "m_null": int(len(pn)), "applicable_null": round(len(pn) / len(nul), 3),
                "power_raw": float(np.mean(pa < 0.05)) if len(pa) else None,
                "power_size_adjusted": float(np.mean(pa <= crit)) if len(pa) else None, "m_alt": int(len(pa)), "crit_p": crit,
            }
        summ.append(rec)
    worst = {t: max((s[t]["type1"] for s in summ if s[t].get("type1") is not None), default=None) for t in tests}
    qual = [t for t in tests if t.startswith("loo_") and worst[t] is not None and worst[t] <= 0.06]
    pw = {t: float(np.mean([s[t]["power_size_adjusted"] for s in summ if s[t].get("power_size_adjusted") is not None])) for t in qual}
    out = {"summary": summ, "max_type1": worst, "qualifying": qual, "mean_size_adjusted_power": pw,
           "adopted": max(pw, key=lambda k: pw[k]) if pw else None, "runtime_s": time.time() - t0}
    (R / "step19_stage1b.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ("max_type1", "qualifying", "adopted")}, indent=2))


if __name__ == "__main__":
    if sys.argv[1:] == ["summarise"]:
        summarise()
    else:
        main()
