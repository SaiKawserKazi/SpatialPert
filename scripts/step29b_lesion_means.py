"""Priority 1 follow-up: are lesion-mean residuals centred and exchangeable across niches?

Pools lesions from many null draws (control model without the pseudo-knockout, leave-out labels)
and reports per niche: number of lesions, mean and SD of the lesion-mean score along random
directions, and the between-lesion share of variance; plus the spatial range of residual
correlation (lesion-mean correlation versus seed distance). Usage: step29b_lesion_means.py tissue.
Writes reports/step29b_<tissue>_lesion_means.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402
from step22_brain import load_brain  # noqa: E402

from spatialpert.effects import fit_control_model, patch_domains, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    tissue = sys.argv[1]
    a = load_brain() if tissue == "brain" else load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a, n_pcs=20 if tissue == "brain" else 30)
    ok = nm.labels >= 0
    ci = np.where((a.is_control & ok) if tissue != "brain" else ok)[0]
    model = fit_control_model(a.y[ci], nm.labels[ci], a.cov[ci])
    tree = cKDTree(nm.coords[ci])
    rng = np.random.default_rng(291)
    ref = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]])
    sd = ref.std(0)
    dirs = rng.normal(size=(20, a.y.shape[1]))
    dirs /= np.sqrt((dirs**2).mean(1, keepdims=True))
    lm, ln, lxy, wvar = [], [], [], []
    for _ in range(60):
        il, les = grow_lesions(tree, nm.coords[ci], 400, 20, rng)
        idx = ci[il]
        loo = nm.leave_out_labels(idx, idx)
        keep = loo >= 0
        idx, loo, les = idx[keep], loo[keep], np.unique(les[keep], return_inverse=True)[1]
        r = residuals(model, a.y[idx], loo, a.cov[idx]) / sd
        sc = r @ dirs.T / r.shape[1]  # cells x directions
        ud = patch_domains(les, loo)
        for k in range(ud.size):
            m = les == k
            lm.append(sc[m].mean(0))
            wvar.append(sc[m].var(0, ddof=1).mean() / m.sum())
            ln.append(int(ud[k]))
            lxy.append(nm.coords[idx][m].mean(0))
    lm_a, ln_a, lxy_a = np.array(lm), np.array(ln), np.array(lxy)
    out: dict[str, object] = {"lesions": int(ln_a.size), "per_niche": {}}
    for d in np.unique(ln_a):
        x = lm_a[ln_a == d]
        out["per_niche"][int(d)] = {  # type: ignore[index]
            "n_lesions": int(x.shape[0]),
            "mean_of_lesion_means": float(x.mean()),
            "t_mean_vs_0": float(x.mean(0).mean() / (x.std(ddof=1) / np.sqrt(x.shape[0]))),
            "sd_lesion_means": float(x.std(ddof=1)),
        }
    out["sampling_var_of_lesion_mean"] = float(np.mean(wvar))
    out["between_lesion_var"] = float(lm_a.var(0, ddof=1).mean())
    # spatial correlation of lesion means: pairs binned by centroid distance (lesions from any draw)
    tr = cKDTree(lxy_a)
    z = (lm_a - lm_a.mean(0)) / lm_a.std(0)
    bins = [(0, 200), (200, 500), (500, 1000), (1000, 2000)]
    corr = {}
    for lo, hi in bins:
        pairs = np.array([(i, j) for i, j in tr.query_pairs(hi) if np.hypot(*(lxy_a[i] - lxy_a[j])) >= lo])
        if len(pairs) > 20:
            corr[f"{lo}-{hi}um"] = float(np.mean(z[pairs[:, 0]] * z[pairs[:, 1]]))
    out["lesion_mean_correlation_by_distance"] = corr
    (ROOT / f"reports/step29b_{tissue}_lesion_means.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
