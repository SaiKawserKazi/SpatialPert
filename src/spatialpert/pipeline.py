"""Shared loaders that rebuild the exact Step 4 inputs for later steps."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA

from spatialpert.domains import (
    CONTROL,
    MULTI,
    neighbourhood_mean,
    neighbourhood_mean_query,
    normalise_log1p,
)
from spatialpert.xenium import read_xenium_h5

COVARIATES = ["log_total", "host_frac50", "density50"]


@dataclass(frozen=True)
class XeniumAnalysis:
    """Tumour cells with a domain, their analysis-gene expression and covariates.

    Attributes:
        cells: Per-cell table (Step 4 columns), tumour cells with domain >= 0 only.
        y: Cells x analysis genes, log-normalised (CSR).
        genes: Analysis gene names.
        cov: Cells x covariates (NaN host fraction filled with its median).
    """

    cells: pd.DataFrame
    y: sparse.csr_matrix
    genes: list[str]
    cov: np.ndarray

    @property
    def is_control(self) -> np.ndarray:
        """Boolean mask of guide-negative cells."""
        mask: np.ndarray = (self.cells["label"] == CONTROL).to_numpy()
        return mask

    def knockout_mask(self, gene: str) -> np.ndarray:
        """Boolean mask of single-guide cells for one knockout."""
        mask: np.ndarray = (self.cells["label"] == gene).to_numpy()
        return mask

    @property
    def is_perturbed(self) -> np.ndarray:
        """Cells with exactly one guide."""
        mask: np.ndarray = (~self.cells["label"].isin([CONTROL, MULTI])).to_numpy()
        return mask


def neighbourhood_features(a: XeniumAnalysis, *, n_pcs: int = 30) -> np.ndarray:
    """Recompute control-defined neighbourhood features for the analysis cells.

    Same recipe as Step 4 (30 PCs of control cells, SpaGCN kernel, l = 50 µm, 100 µm
    radius, self excluded), restricted to the analysis cells.

    Args:
        a: Analysis bundle.
        n_pcs: Principal components.

    Returns:
        Cells x n_pcs features (NaN where no control neighbour).
    """
    ctrl = a.is_control
    pcs = PCA(n_components=n_pcs, svd_solver="arpack", random_state=0).fit_transform(a.y[ctrl])
    coords = a.cells[["x", "y"]].to_numpy()
    feats, _ = neighbourhood_mean(coords, pcs, ctrl, radius=100.0, length=50.0)
    return feats


@dataclass
class NicheModel:
    """Control-defined niches that can be re-derived with some cells left out (review A2).

    Attributes:
        pcs: Analysis cells x PCs (control PCA applied to every cell).
        is_ref: Cells that serve as niche references (controls with valid PCs).
        coords: Cells x 2.
        centroids: K-means centroids in neighbourhood-feature space.
        labels: Niche per cell using all references (-1 where undefined).
    """

    pcs: np.ndarray
    is_ref: np.ndarray
    coords: np.ndarray
    centroids: np.ndarray
    labels: np.ndarray

    @classmethod
    def fit(cls, a: XeniumAnalysis, *, k: int = 4, n_pcs: int = 30, seed: int = 0) -> NicheModel:
        """Fit PCA and k-means on control cells (same recipe as Step 4)."""
        from sklearn.cluster import KMeans  # noqa: PLC0415

        ctrl = a.is_control
        pca = PCA(n_components=n_pcs, svd_solver="arpack", random_state=seed).fit(a.y[ctrl])
        pcs = pca.transform(a.y)
        coords = a.cells[["x", "y"]].to_numpy()
        feats, _ = neighbourhood_mean(coords, pcs[ctrl], ctrl, radius=100.0, length=50.0)
        ok = ~np.isnan(feats).any(axis=1)
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(feats[ctrl & ok])
        labels = np.full(len(coords), -1)
        labels[ok] = km.predict(feats[ok])
        return cls(
            pcs=pcs, is_ref=ctrl, coords=coords, centroids=km.cluster_centers_, labels=labels
        )

    def leave_out_labels(self, query: np.ndarray, exclude: np.ndarray) -> np.ndarray:
        """Niche of ``query`` cells with ``exclude`` cells removed from the reference.

        Args:
            query: Indices of cells to label.
            exclude: Indices of cells not allowed as references (e.g. the pseudo-knockout).

        Returns:
            Niche per query cell (-1 where no reference neighbour remains).
        """
        ref = self.is_ref.copy()
        ref[exclude] = False
        feats = neighbourhood_mean_query(self.coords[query], self.coords[ref], self.pcs[ref])
        out = np.full(query.size, -1)
        ok = ~np.isnan(feats).any(axis=1)
        d2 = ((feats[ok][:, None, :] - self.centroids[None]) ** 2).sum(-1)
        out[ok] = d2.argmin(1)
        return out


def tile_ids(coords: np.ndarray, size: float = 250.0) -> np.ndarray:
    """Integer id of the square tile each cell falls in.

    Args:
        coords: Cells x 2 (µm).
        size: Tile edge (µm).

    Returns:
        Consecutive tile ids (0..T-1) for the given cells.
    """
    ij = np.floor(coords / size).astype(np.int64)
    key = ij[:, 0] * 1_000_003 + ij[:, 1]
    _, ids = np.unique(key, return_inverse=True)
    out: np.ndarray = ids
    return out


def load_xenium_analysis(root: Path) -> XeniumAnalysis:
    """Rebuild log-normalised expression for Step 4 tumour cells.

    Normalisation uses all cells kept in Step 4 (same median target), then subsets.

    Args:
        root: Project root.

    Returns:
        The analysis bundle.
    """
    interim = root / "data" / "interim"
    cells_all = pd.read_parquet(interim / "xenium_cells.parquet")
    genes = pd.read_csv(interim / "analysis_genes.csv", header=None)[0].astype(str).tolist()

    x = read_xenium_h5(
        root / "data" / "raw" / "GSE319123" / "GSM9511110_Xenium_lung_cell_feature_matrix.h5"
    )
    m = x.matrix
    eidx = [
        i
        for i, (fid, t) in enumerate(zip(m.feature_ids, x.feature_types, strict=True))
        if t == "Gene Expression" and fid.startswith("ENS")
    ]
    pos = {b: i for i, b in enumerate(m.barcodes)}
    rows = np.array([pos[c] for c in cells_all["cell_id"]])
    expr = m.counts[rows][:, eidx].tocsr()
    x_log = normalise_log1p(expr)
    names = [m.feature_names[i] for i in eidx]
    col = {g: i for i, g in enumerate(names)}
    x_log = x_log[:, [col[g] for g in genes]].tocsr()

    keep = (cells_all["tumour"] & (cells_all["domain"] >= 0)).to_numpy()
    cells = cells_all.loc[keep].reset_index(drop=True)
    cov = cells[COVARIATES].to_numpy(dtype=np.float64).copy()
    for j in range(cov.shape[1]):
        bad = np.isnan(cov[:, j])
        if bad.any():
            cov[bad, j] = np.nanmedian(cov[:, j])
    return XeniumAnalysis(cells=cells, y=x_log[keep], genes=genes, cov=cov)
