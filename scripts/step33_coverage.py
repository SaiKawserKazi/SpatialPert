"""Review 3, Priority 1: bias, RMSE and 95% CI coverage of the knockout-by-niche contrast Psi.

Pseudo-knockouts (400 cells, L lesions, leave-out labels, control model without the
pseudo-knockout) get a shared effect delta*sd*u in all cells plus lambda*delta*sd*v in niche b,
so the true contrast is Psi_g = lambda*delta*sd_g*v_g between niche b and niche a (the two
niches with most lesions). Intervals: CR2 clustered on lesions (Imbens-Kolesar df) and naive
cell-level OLS. Usage: step33_coverage.py tumour|brain. Writes reports/step33_<tissue>_coverage.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402
from step22_brain import load_brain  # noqa: E402

from spatialpert.benchmark import unit_rms
from spatialpert.effects import fit_control_model, patch_domains, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import niche_contrast

ROOT = Path(__file__).resolve().parents[1]
DELTA = 0.5
REPS = 200


def naive_ci(r: np.ndarray, nic: np.ndarray, a: int, b: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ra, rb = r[nic == a], r[nic == b]
    est = rb.mean(0) - ra.mean(0)
    se = np.sqrt(ra.var(0, ddof=1) / len(ra) + rb.var(0, ddof=1) / len(rb))
    z = stats.norm.ppf(0.975)
    return est, est - z * se, est + z * se


def main() -> None:
    tissue = sys.argv[1]
    a_ = load_brain() if tissue == "brain" else load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a_, n_pcs=20 if tissue == "brain" else 30)
    ok = nm.labels >= 0
    ci = np.where((a_.is_control & ok) if tissue == "tumour" else ok)[0]
    model0 = fit_control_model(a_.y[ci], nm.labels[ci], a_.cov[ci])
    sd = residuals(model0, a_.y[ci[:20000]], nm.labels[ci[:20000]], a_.cov[ci[:20000]]).std(0)
    tree = cKDTree(nm.coords[ci])
    rng = np.random.default_rng(33)
    out: dict[str, object] = {}
    for n_les in (5, 10, 20):
        for lam in (0.0, 0.5):
            rec: dict[str, list[float]] = {k: [] for k in ("cov_cr2", "cov_naive", "bias", "rmse", "width_cr2", "width_naive")}
            tried = 0
            while len(rec["cov_cr2"]) < REPS and tried < 3 * REPS:
                tried += 1
                il, les = grow_lesions(tree, nm.coords[ci], 400, n_les, rng)
                idx = ci[il]
                d = nm.leave_out_labels(idx, idx)
                k = d >= 0
                idx, d, les = idx[k], d[k], np.unique(les[k], return_inverse=True)[1]
                ud = patch_domains(les, d)
                labs, cnt = np.unique(ud, return_counts=True)
                if (cnt >= 2).sum() < 2:
                    continue
                a, b = labs[np.argsort(-cnt)[:2]]
                rest = np.setdiff1d(ci, idx)
                model = fit_control_model(a_.y[rest], nm.labels[rest], a_.cov[rest])
                r = residuals(model, a_.y[idx], d, a_.cov[idx]) / sd
                v = unit_rms(rng.normal(size=r.shape[1]))
                truth = lam * DELTA * v
                r = r + DELTA * unit_rms(rng.normal(size=r.shape[1]))
                r[d == b] += truth
                fit = niche_contrast(r, d, les, int(a), int(b))
                if fit is None:
                    continue
                est, nlo, nhi = naive_ci(r, d, int(a), int(b))
                rec["cov_cr2"].append(float(np.mean((fit["lo"] <= truth) & (truth <= fit["hi"]))))
                rec["cov_naive"].append(float(np.mean((nlo <= truth) & (truth <= nhi))))
                rec["bias"].append(float(np.mean(fit["estimate"] - truth)))
                rec["rmse"].append(float(np.sqrt(np.mean((fit["estimate"] - truth) ** 2))))
                rec["width_cr2"].append(float(np.median(fit["hi"] - fit["lo"])))
                rec["width_naive"].append(float(np.median(nhi - nlo)))
            out[f"L{n_les}_lam{lam}"] = {kk: float(np.mean(vv)) for kk, vv in rec.items()} | {
                "draws": len(rec["cov_cr2"]), "estimable_fraction": len(rec["cov_cr2"]) / tried}
            print(tissue, n_les, lam, out[f"L{n_les}_lam{lam}"], flush=True)
    (ROOT / f"reports/step33_{tissue}_coverage.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
