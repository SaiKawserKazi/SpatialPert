"""A7: within-lesion ICC of the tested score in control pseudo-lesions (tumour tissue).

Score per cell = mean over genes of (residual / gene SD) * v, v a random unit-RMS direction,
i.e. the projection along which a niche effect would be injected. Writes reports/step23_icc.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402

from spatialpert.benchmark import unit_rms
from spatialpert.design import icc_oneway
from spatialpert.effects import fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    ok = a.is_control & (nm.labels >= 0)
    model = fit_control_model(a.y[ok], nm.labels[ok], a.cov[ok])
    ci = np.where(ok)[0]
    cxy = nm.coords[ci]
    tree = cKDTree(cxy)
    gsd = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    rng = np.random.default_rng(23)
    out = {}
    for m in (20, 40, 80, 133, 400):
        iccs, sig, pg = [], [], []
        for _ in range(40):
            n_les = max(2, 2000 // m)
            il, les = grow_lesions(tree, cxy, n_les * m, n_les, rng)
            idx = ci[il]
            z = residuals(model, a.y[idx], nm.labels[idx], a.cov[idx]) / gsd
            s = z @ unit_rms(rng.normal(size=z.shape[1])) / z.shape[1]
            iccs.append(icc_oneway(s, les))
            sig.append(float(s.std()))
            pg.append(float(np.median([icc_oneway(z[:, g], les) for g in range(z.shape[1])])))
        out[str(m)] = {"icc_median": float(np.median(iccs)), "icc_iqr": np.quantile(iccs, [0.25, 0.75]).tolist(), "sigma_median": float(np.median(sig)), "per_gene_icc_median": float(np.median(pg))}
        print(m, out[str(m)], flush=True)
    (ROOT / "reports/step23_icc.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
