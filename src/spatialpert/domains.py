"""Tumour-cell calling, knockout calls, neighbourhood features and domains (protocol §2-3)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

CONTROL = "control"
MULTI = "multi"


def normalise_log1p(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    """Scale each cell to the median total count, then apply log1p (protocol §2.1).

    Args:
        counts: Cells x genes raw counts.

    Returns:
        Normalised, log-transformed matrix (CSR, float64).
    """
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    target = float(np.median(totals[totals > 0]))
    scale = np.divide(target, totals, out=np.zeros_like(totals), where=totals > 0)
    out = sparse.diags(scale) @ counts.astype(np.float64)
    out = sparse.csr_matrix(out)
    out.data = np.log1p(out.data)
    return out


def marker_score(x_log: sparse.csr_matrix, columns: Sequence[int]) -> np.ndarray:
    """Mean z-score of marker genes per cell.

    Args:
        x_log: Normalised log matrix.
        columns: Marker gene columns.

    Returns:
        One score per cell.
    """
    sub = np.asarray(x_log[:, list(columns)].todense(), dtype=np.float64)
    sd = sub.std(axis=0)
    sd[sd == 0] = 1.0
    z = (sub - sub.mean(axis=0)) / sd
    return np.asarray(z.mean(axis=1)).ravel()


def call_tumour(score: np.ndarray, *, posterior: float = 0.9, seed: int = 0) -> np.ndarray:
    """Two-component Gaussian mixture on the tumour score (protocol §2.2).

    Args:
        score: Tumour score per cell.
        posterior: Minimum posterior of the high component.
        seed: Random seed.

    Returns:
        Boolean tumour flag per cell.
    """
    gmm = GaussianMixture(n_components=2, random_state=seed).fit(score.reshape(-1, 1))
    high = int(np.argmax(gmm.means_.ravel()))
    post = gmm.predict_proba(score.reshape(-1, 1))[:, high]
    flag: np.ndarray = post > posterior
    return flag


def knockout_calls(guides: sparse.csr_matrix, targets: Sequence[str]) -> np.ndarray:
    """Label cells by guide content (protocol §2.3).

    Args:
        guides: Cells x guide-probe counts.
        targets: Target gene per guide probe column.

    Returns:
        Per cell: the target gene (exactly one guide detected), ``"control"`` (none) or
        ``"multi"`` (two or more different targets).
    """
    pos = sparse.csr_matrix(guides > 0)
    n_pos = np.asarray(pos.sum(axis=1)).ravel()
    labels = np.full(guides.shape[0], CONTROL, dtype=object)
    labels[n_pos >= 2] = MULTI
    single = np.where(n_pos == 1)[0]
    if single.size:
        cols = pos[single].indices
        labels[single] = np.asarray(targets, dtype=object)[cols]
    return labels


def neighbourhood_mean(
    coords: np.ndarray,
    values: np.ndarray,
    reference: np.ndarray,
    *,
    radius: float = 100.0,
    length: float = 50.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Kernel-weighted mean of reference cells' values around every cell, excluding itself.

    Weight w = exp(-d^2 / 2 l^2), as in SpaGCN (protocol §3.2).

    Args:
        coords: Cells x 2 positions (µm).
        values: Reference cells x features (rows aligned with ``np.where(reference)``).
        reference: Boolean mask of cells that may contribute.
        radius: Neighbourhood radius (µm).
        length: Kernel length scale l (µm).

    Returns:
        (features, total weight) per cell. Cells with no reference neighbour get NaN features.
    """
    ref_idx = np.where(reference)[0]
    tree_all = cKDTree(coords)
    tree_ref = cKDTree(coords[ref_idx])
    dist = tree_all.sparse_distance_matrix(tree_ref, radius, output_type="coo_matrix")
    rows, cols, d = dist.row, dist.col, dist.data
    keep = rows != ref_idx[cols]  # drop self-pairs (zero distances are absent anyway)
    rows, cols, d = rows[keep], cols[keep], d[keep]
    w = np.exp(-(d**2) / (2 * length**2))
    wmat = sparse.csr_matrix((w, (rows, cols)), shape=(coords.shape[0], ref_idx.size))
    total = np.asarray(wmat.sum(axis=1)).ravel()
    feats = np.asarray(wmat @ values, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        feats = feats / total[:, None]
    feats[total == 0] = np.nan
    return feats, total


def neighbourhood_mean_query(
    query_xy: np.ndarray,
    ref_xy: np.ndarray,
    ref_values: np.ndarray,
    *,
    radius: float = 100.0,
    length: float = 50.0,
) -> np.ndarray:
    """Kernel-weighted mean of reference values around query points (no self-exclusion needed).

    Used for leave-out niche features: the caller removes the excluded cells from the
    reference before calling.

    Args:
        query_xy: Query cells x 2.
        ref_xy: Reference cells x 2.
        ref_values: Reference cells x features.
        radius: Neighbourhood radius (µm).
        length: Kernel length scale (µm).

    Returns:
        Query x features (NaN rows where no reference cell lies within ``radius``).
    """
    dist = cKDTree(query_xy).sparse_distance_matrix(
        cKDTree(ref_xy), radius, output_type="coo_matrix"
    )
    w = np.exp(-(dist.data**2) / (2 * length**2))
    wm = sparse.csr_matrix((w, (dist.row, dist.col)), shape=(query_xy.shape[0], ref_xy.shape[0]))
    tot = np.asarray(wm.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.asarray(wm @ ref_values, dtype=np.float64) / tot[:, None]
    out[tot == 0] = np.nan
    return out


def neighbour_counts(
    coords: np.ndarray, flag: np.ndarray, *, radius: float = 50.0
) -> tuple[np.ndarray, np.ndarray]:
    """Cells within ``radius`` (excluding self), and the share of them not flagged.

    Args:
        coords: Cells x 2 positions (µm).
        flag: Boolean per cell (tumour).
        radius: Radius (µm).

    Returns:
        (density, fraction of unflagged neighbours); fraction is NaN with no neighbours.
    """
    tree = cKDTree(coords)
    pairs = tree.sparse_distance_matrix(tree, radius, output_type="coo_matrix")
    rows, cols = pairs.row, pairs.col
    keep = rows != cols
    rows, cols = rows[keep], cols[keep]
    n = coords.shape[0]
    density = np.bincount(rows, minlength=n).astype(np.float64)
    other = np.bincount(rows, weights=(~flag[cols]).astype(np.float64), minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = other / density
    frac[density == 0] = np.nan
    return density, frac


def choose_domains(
    features: np.ndarray,
    *,
    ks: Sequence[int] = (3, 4, 5, 6),
    subsample: int = 20_000,
    seed: int = 0,
) -> tuple[KMeans, dict[int, float]]:
    """Fit k-means for each K and keep the best by silhouette on a subsample (protocol §3.3).

    Args:
        features: Neighbourhood features of control tumour cells (no NaN rows).
        ks: Candidate numbers of domains.
        subsample: Cells used for the silhouette.
        seed: Random seed.

    Returns:
        The chosen fitted model and the silhouette per K.
    """
    rng = np.random.default_rng(seed)
    sub = rng.choice(features.shape[0], size=min(subsample, features.shape[0]), replace=False)
    scores: dict[int, float] = {}
    models: dict[int, KMeans] = {}
    for k in ks:
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(features)
        models[k] = km
        scores[k] = float(silhouette_score(features[sub], km.labels_[sub]))
    best = max(scores, key=lambda k: scores[k])
    return models[best], scores
