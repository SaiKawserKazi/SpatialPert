"""SpatialPert tutorial: a 2-minute end-to-end run on simulated tissue.

Run:  uv run python docs/tutorial.py
For real data, build the same five inputs from your AnnData (see `spatialpert.api.from_anndata`).
"""

import numpy as np

from spatialpert.api import location_dependence

rng = np.random.default_rng(0)

# 1. A tissue of 3,000 control cells in two niches (left / right half), 8 genes.
n_ctrl, genes = 3000, 8
coords = rng.uniform(0, 3000, size=(n_ctrl, 2))
niche = (coords[:, 0] > 1500).astype(int)
y = rng.normal(size=(n_ctrl, genes)) + niche[:, None] * 0.5
labels = np.array(["control"] * n_ctrl, dtype=object)

# 2. Two knockouts growing as clones (patches of 10 cells), 15 clones per niche each.
def clones(effect_in_niche1: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xy, nic, expr = [], [], []
    for d in (0, 1):
        centres = rng.uniform([200 + 1500 * d, 200], [1300 + 1500 * d, 2800], size=(15, 2))
        for c in centres:
            xy.append(c + rng.normal(0, 6, size=(10, 2)))
            nic.append(np.full(10, d))
            e = rng.normal(size=(10, genes)) + d * 0.5 + 0.8      # shared core effect 0.8
            e[:, :3] += effect_in_niche1 * d                     # niche-specific part
            expr.append(e)
    return np.vstack(xy), np.concatenate(nic), np.vstack(expr)

xy_a, n_a, e_a = clones(1.5)   # KO_A: effect depends on niche
xy_b, n_b, e_b = clones(0.0)   # KO_B: same effect everywhere

Y = np.vstack([y, e_a, e_b])
XY = np.vstack([coords, xy_a, xy_b])
LAB = np.concatenate([labels, ["KO_A"] * len(n_a), ["KO_B"] * len(n_b)])
NICHE = np.concatenate([niche, n_a, n_b])

# 3. One call. Here each simulated clone is its own lesion, so patch units are used; for real
#    data the default lesion unit (100 um single linkage) is recommended.
res = location_dependence(Y, XY, LAB, NICHE, genes=[f"g{i}" for i in range(genes)], unit="patch", n_perm=500)
print(res.table.to_string(index=False))
print("pooled p:", round(res.pooled_p, 4))
print(res.effects["KO_A"].round(2))
# Expected: KO_A p < 0.01 (niche-dependent), KO_B p well above 0.05.
