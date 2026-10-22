"""Priority 1: why does lesion permutation inflate in the brain? (null draws only)

Variants on the same pseudo-knockout draws (lesion permutation, leave-out labels unless noted):
  v0      control model fitted on all control bins (as step22, pseudo-knockout included)
  v1      control model refitted without the pseudo-knockout bins (matches a real screen,
          where knockout cells never enter the control fit)
  v2      v1 + local-control residuals (50 nearest controls, pseudo-knockout excluded)
  v3      v0 with circular labels
Per draw also records the fraction of bins relabelled by leave-out and the mean squared norm
of residuals of relabelled vs kept bins. Usage: step29_brain_diagnosis.py [tumour|brain].
Writes reports/step29_<tissue>_diagnosis.csv / .json.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402
from step22_brain import load_brain  # noqa: E402

from spatialpert.effects import (
    clone_patches,
    covariate_adjusted,
    fit_control_model,
    local_residuals,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.baselines import patch_means, studentized_permutation_p
from spatialpert.pipeline import NicheModel, load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
N_PERM = 250
TAG = sys.argv[4] if len(sys.argv) > 4 else ""
LS = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (5, 10, 20)
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 200


def perm_p(r: np.ndarray, units: np.ndarray, niche: np.ndarray, seed: int) -> float:
    ud = patch_domains(units, niche)
    return float("nan") if np.unique(ud).size < 2 else patch_permutation_test(r, units, ud, n_perm=N_PERM, seed=seed).p_value


def stud(r: np.ndarray, units: np.ndarray, niche: np.ndarray, seed: int) -> float:
    ud = patch_domains(units, niche)
    if np.unique(ud).size < 2:
        return float("nan")
    mns, _, g = patch_means(r, units, ud)
    labs, cnt = np.unique(g, return_counts=True)
    use = np.isin(g, labs[cnt >= 2])  # Welch needs >= 2 lesions per niche
    if np.unique(g[use]).size < 2:
        return float("nan")
    return float(studentized_permutation_p(mns[use], g[use], n_perm=N_PERM, seed=seed))


def sd_ratio(r: np.ndarray, units: np.ndarray, niche: np.ndarray) -> float:
    """Max/min across niches of the SD of lesion-mean scores (heteroscedasticity of lesion baselines)."""
    ud = patch_domains(units, niche)
    mns, _, g = patch_means(r, units, ud)
    sc = mns.mean(1)
    sds = [sc[g == k].std(ddof=1) for k in np.unique(g) if (g == k).sum() >= 2]
    return float(max(sds) / max(min(sds), 1e-12)) if len(sds) >= 2 else float("nan")


def pure(r: np.ndarray, units: np.ndarray, niche: np.ndarray, seed: int) -> float:
    """Lesion permutation using only cells whose niche equals their lesion's majority niche."""
    keep = niche == patch_domains(units, niche)[units]
    u = np.unique(units[keep], return_inverse=True)[1]
    return perm_p(r[keep], u, niche[keep], seed)


def maj_resid(model, a, idx: np.ndarray, units: np.ndarray, niche: np.ndarray) -> np.ndarray:  # noqa: ANN001
    """Residuals of every cell against its lesion's majority-niche mean (not its own niche)."""
    return residuals(model, a.y[idx], patch_domains(units, niche)[units], a.cov[idx])


def main() -> None:
    tissue = sys.argv[1] if len(sys.argv) > 1 else "brain"
    a = load_brain() if tissue == "brain" else load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a, n_pcs=20 if tissue == "brain" else 30)
    ok = nm.labels >= 0
    ctrl = (a.is_control & ok) if tissue != "brain" else ok
    ci = np.where(ctrl)[0]
    model0 = fit_control_model(a.y[ci], nm.labels[ci], a.cov[ci])
    tree = cKDTree(nm.coords[ci])
    rng = np.random.default_rng(29)
    rows = []
    for n_les in LS:
        for rep in range(REPS):
            il, les = grow_lesions(tree, nm.coords[ci], 400, n_les, rng)
            idx = ci[il]
            circ = nm.labels[idx]
            loo = nm.leave_out_labels(idx, idx)
            keep = loo >= 0
            idx, circ, loo, les = idx[keep], circ[keep], loo[keep], np.unique(les[keep], return_inverse=True)[1]
            seed = int(rng.integers(1 << 31))
            excl = np.zeros(a.y.shape[0], dtype=bool)
            excl[idx] = True
            rest = ci[~excl[ci]]
            model1 = fit_control_model(a.y[rest], nm.labels[rest], a.cov[rest])
            r0 = residuals(model0, a.y[idx], loo, a.cov[idx])
            r1 = residuals(model1, a.y[idx], loo, a.cov[idx])
            z_ko = covariate_adjusted(model1, a.y[idx], a.cov[idx])
            z_ctrl = covariate_adjusted(model1, a.y[rest], a.cov[rest])
            r2 = local_residuals(z_ko, nm.coords[idx], z_ctrl, nm.coords[rest], k=50)
            rel = circ != loo
            nrm = (r1**2).sum(1)
            rows.append({
                "lesions": n_les, "rep": rep, "frac_relabelled": float(rel.mean()),
                "norm_relabelled": float(nrm[rel].mean()) if rel.any() else np.nan,
                "norm_kept": float(nrm[~rel].mean()) if (~rel).any() else np.nan,
                "n_niches_loo": int(np.unique(patch_domains(les, loo)).size),
                "v0": perm_p(r0, les, loo, seed), "v1": perm_p(r1, les, loo, seed),
                "v2": perm_p(r2, les, loo, seed), "v3": perm_p(residuals(model0, a.y[idx], circ, a.cov[idx]), les, circ, seed),
                "patch30_v1": perm_p(r1, clone_patches(nm.coords[idx], link=30.0), loo, seed),
                "v4": stud(r1, les, loo, seed), "v5": stud(r2, les, loo, seed),
                "lesion_sd_ratio": sd_ratio(r1, les, loo),
                "v6": pure(r1, les, loo, seed), "v7": perm_p(maj_resid(model1, a, idx, les, loo), les, loo, seed),
                "v8": pure(r2, les, loo, seed), "frac_minority": float(np.mean(loo != patch_domains(les, loo)[les])),
            })
        print(tissue, "L", n_les, flush=True)
        pd.DataFrame(rows).to_csv(R / f"step29_{tissue}_diagnosis{TAG}.csv", index=False)
    df = pd.DataFrame(rows)
    out = {}
    for L in LS:
        s = df[df.lesions == L]
        rec = {"frac_relabelled_mean": float(s.frac_relabelled.mean()),
               "norm_ratio_relabelled_vs_kept": float((s.norm_relabelled / s.norm_kept).median())}
        for v in ("v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "patch30_v1"):
            p = s[v].dropna()
            rec[v] = {"type1": float((p < 0.05).mean()) if len(p) else None, "m": int(len(p))}
        # does inflation track relabelling? type-I (v1) in draws below / above the median relabelled fraction
        med = s.frac_relabelled.median()
        for tag, sub in (("low_relabel", s[s.frac_relabelled <= med]), ("high_relabel", s[s.frac_relabelled > med])):
            p = sub["v1"].dropna()
            rec[f"v1_{tag}"] = float((p < 0.05).mean()) if len(p) else None
        out[L] = rec
    (R / f"step29_{tissue}_diagnosis{TAG}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
