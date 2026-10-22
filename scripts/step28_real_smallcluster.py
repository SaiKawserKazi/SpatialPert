"""Review 2 (A1 on real data): small-cluster tests on the four Xenium knockouts.

Lesions are 100 um single-linkage clusters; niches are the benchmark caller's labels (NicheModel:
k-means K = 4 on control neighbourhood features), the same caller as the stage-1b benchmark. The
permutation tests are re-run on these labels so every real-data p-value uses one niche definition.
Writes reports/step28_real_smallcluster.json.
"""

import json
from pathlib import Path

import numpy as np

from spatialpert.api import location_dependence
from spatialpert.benchmark import cell_permutation_p
from spatialpert.effects import clone_patches, fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import cr2_two_niche_min_q, mixedlm_min_q, wild_cluster_bootstrap_p, within_lesion_signflip_p

ROOT = Path(__file__).resolve().parents[1]
KOS = ["CIITA", "DIO3", "ESAM", "PRSS23"]


def main() -> None:
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    lab = a.cells["label"].to_numpy().astype(object)
    ok = nm.labels >= 0
    ctrl = a.is_control & ok
    model = fit_control_model(a.y[ctrl], nm.labels[ctrl], a.cov[ctrl])
    out: dict[str, object] = {}
    for k in KOS:
        ix = np.where((lab == k) & ok)[0]
        dd = nm.labels[ix]
        les = np.unique(clone_patches(nm.coords[ix], link=100.0), return_inverse=True)[1]
        r = residuals(model, a.y[ix], dd, a.cov[ix])
        per = {int(d): int(np.unique(les[dd == d]).size) for d in np.unique(dd)}
        cells = {int(d): int((dd == d).sum()) for d in np.unique(dd)}
        out[k] = {
            "cells": int(ix.size), "lesions": int(np.unique(les).size),
            "lesions_per_niche": per, "cells_per_niche": cells,
            "cell_p": cell_permutation_p(r, dd, n_perm=999, seed=0),
            "wcb_p": wild_cluster_bootstrap_p(r, dd, les, n_boot=999, seed=0),
            "cr2_min_q": cr2_two_niche_min_q(r, dd, les),
            "signflip_p": within_lesion_signflip_p(r, dd, les, seed=0),
            "mixedlm_min_q": mixedlm_min_q(r, dd, les),
        }
        print(k, out[k], flush=True)
    keep = ok & np.isin(lab, ["control", *KOS])
    for unit, kw in (("patch30", {"unit": "patch", "link_um": 30.0}), ("lesion", {"unit": "lesion"}),
                     ("stud_lesion", {"unit": "lesion", "statistic": "studentized"})):
        res = location_dependence(a.y[keep], nm.coords[keep], lab[keep], nm.labels[keep], covariates=a.cov[keep],
                                  genes=a.genes, knockouts=KOS, n_perm=1000, **kw)  # type: ignore[arg-type]
        for row in res.table.to_dict(orient="records"):
            out[row["knockout"]][f"{unit}_p"] = row["p"]  # type: ignore[index]
        print(unit, res.table[["knockout", "units", "p"]].to_dict(orient="records"), flush=True)
    (ROOT / "reports/step28_real_smallcluster.json").write_text(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
