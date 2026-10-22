"""Review 3, Priority 1: is the CR2 contrast test calibrated where lesion permutation is not?

Null pseudo-knockouts in the brain (and tumour for reference), leave-out labels, control model
without the pseudo-knockout; CR2 two-niche test (minimum BH q over genes) at 5, 10, 20 lesions.
Usage: step34_brain_cr2.py tissue. Writes reports/step34_<tissue>_cr2.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402
from step22_brain import load_brain  # noqa: E402

from spatialpert.effects import fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import cr2_two_niche_min_q

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    tissue = sys.argv[1]
    a = load_brain() if tissue == "brain" else load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a, n_pcs=20 if tissue == "brain" else 30)
    ok = nm.labels >= 0
    ci = np.where((a.is_control & ok) if tissue == "tumour" else ok)[0]
    tree = cKDTree(nm.coords[ci])
    rng = np.random.default_rng(34)
    out = {}
    for n_les in (5, 10, 20):
        q = []
        for _ in range(200):
            il, les = grow_lesions(tree, nm.coords[ci], 400, n_les, rng)
            idx = ci[il]
            d = nm.leave_out_labels(idx, idx)
            k = d >= 0
            idx, d, les = idx[k], d[k], les[k]
            rest = np.setdiff1d(ci, idx)
            model = fit_control_model(a.y[rest], nm.labels[rest], a.cov[rest])
            q.append(cr2_two_niche_min_q(residuals(model, a.y[idx], d, a.cov[idx]), d, les))
        qa = np.array(q)
        qa = qa[~np.isnan(qa)]
        out[n_les] = {"type1": float((qa < 0.05).mean()) if qa.size else None, "applicable": qa.size / 200}
        print(tissue, n_les, out[n_les], flush=True)
    (ROOT / f"reports/step34_{tissue}_cr2.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
