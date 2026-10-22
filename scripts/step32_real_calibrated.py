"""Review 3, Priority 1 on real data: calibrated p-values for the Xenium knockouts, and the
lesion-exchangeability diagnostic in tumour and brain.

Writes reports/step32_real_calibrated.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402
from step22_brain import load_brain  # noqa: E402

from spatialpert.api import location_dependence
from spatialpert.calibration import calibrate_p, lesion_exchangeability, matched_null_p
from spatialpert.effects import fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]
KOS = ["CIITA", "DIO3", "ESAM", "PRSS23"]
N_NULL = 200


def diagnostic(a, nm, tissue: str, rng: np.random.Generator) -> dict[str, float]:  # noqa: ANN001
    ok = nm.labels >= 0
    ci = np.where((a.is_control & ok) if tissue == "tumour" else ok)[0]
    model = fit_control_model(a.y[ci], nm.labels[ci], a.cov[ci])
    sd = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    tree = cKDTree(nm.coords[ci])
    vr, sr = [], []
    for _ in range(40):
        il, les = grow_lesions(tree, nm.coords[ci], 400, 20, rng)
        idx = ci[il]
        d = nm.leave_out_labels(idx, idx)
        k = d >= 0
        r = residuals(model, a.y[idx[k]], d[k], a.cov[idx[k]]) / sd
        res = lesion_exchangeability(r, les[k], d[k])
        vr.append(res["variance_ratio"])
        sr.append(res["niche_sd_ratio"])
    return {"variance_ratio_median": float(np.nanmedian(vr)), "niche_sd_ratio_median": float(np.nanmedian(sr))}


def main() -> None:
    rng = np.random.default_rng(32)
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    out: dict[str, object] = {"diagnostic": {"tumour": diagnostic(a, nm, "tumour", rng)}}
    print(out, flush=True)
    ok = nm.labels >= 0
    lab = a.cells["label"].to_numpy().astype(object)
    keep = ok & np.isin(lab, ["control", *KOS])
    kidx = np.where(keep)[0]
    y, xy, lab_k, dom_k, cov = a.y[kidx], nm.coords[kidx], lab[keep], nm.labels[keep], a.cov[keep]

    def relabel(pick: np.ndarray) -> np.ndarray:
        g = kidx[pick]
        return nm.leave_out_labels(g, g)

    real: dict[str, object] = {}
    for unit, kw in (("lesion", {"unit": "lesion"}), ("patch30", {"unit": "patch", "link_um": 30.0})):
        res = location_dependence(y, xy, lab_k, dom_k, covariates=cov, knockouts=KOS, n_perm=1000, **kw)  # type: ignore[arg-type]
        pobs = dict(zip(res.table["knockout"], res.table["p"], strict=True))
        for k in KOS:
            nul = matched_null_p(y, xy, lab_k, dom_k, k, relabel=relabel, n_null=N_NULL, seed=7, covariates=cov, n_perm=250, **kw)
            p = pobs.get(k, float("nan"))
            real[f"{k}_{unit}"] = {"p_permutation": p, "p_calibrated": calibrate_p(p, nul, seed=1),
                                   "null_testable": float(np.mean(~np.isnan(nul))), "null_type1_raw": float(np.nanmean(nul < 0.05))}
            print(k, unit, real[f"{k}_{unit}"], flush=True)
    out["real"] = real
    b = load_brain()
    out["diagnostic"]["brain"] = diagnostic(b, NicheModel.fit(b, n_pcs=20), "brain", rng)  # type: ignore[index]
    (ROOT / "reports/step32_real_calibrated.json").write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
