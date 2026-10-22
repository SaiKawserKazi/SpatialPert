"""A6: export control tumour cells (full panel, log-normalised) as h5ad for external niche callers."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from spatialpert.domains import normalise_log1p
from spatialpert.xenium import read_xenium_h5

ROOT = Path(__file__).resolve().parents[1]
cells = pd.read_parquet(ROOT / "data/interim/xenium_cells.parquet")
t = cells[cells["tumour"] & (cells["domain"] >= 0) & (cells["label"] == "control")].reset_index(drop=True)
x = read_xenium_h5(ROOT / "data/raw/GSE319123/GSM9511110_Xenium_lung_cell_feature_matrix.h5")
m = x.matrix
eidx = [i for i, (fid, ft) in enumerate(zip(m.feature_ids, x.feature_types, strict=True)) if ft == "Gene Expression" and fid.startswith("ENS")]
pos = {b: i for i, b in enumerate(m.barcodes)}
rows = np.array([pos[c] for c in t["cell_id"]])
counts = m.counts[rows][:, eidx].tocsr().astype(np.float32)
a = ad.AnnData(X=normalise_log1p(counts).astype(np.float32), obs=t[["cell_id", "domain"]].set_index("cell_id"),
               var=pd.DataFrame(index=[m.feature_names[i] for i in eidx]))
a.layers["counts"] = counts
a.obsm["spatial"] = t[["x", "y"]].to_numpy()
a.write_h5ad(ROOT / "data/interim/a6/control_tumour.h5ad")
print(a)
