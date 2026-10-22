"""Alternative tests for niche dependence (review Stage 1).

All operate on residuals of one knockout's cells, their niche labels and clone-patch ids.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from spatialpert.effects import benjamini_hochberg


def patch_means(
    resid: np.ndarray, patch: np.ndarray, patch_dom: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Patch-level pseudobulk: mean residual per patch.

    Args:
        resid: Cells x genes residuals.
        patch: Patch id per cell (0..P-1).
        patch_dom: Niche per patch.

    Returns:
        (patch x gene means, cells per patch, niche per patch).
    """
    n_p = patch_dom.size
    sums = np.zeros((n_p, resid.shape[1]))
    np.add.at(sums, patch, resid)
    sizes = np.bincount(patch, minlength=n_p).astype(np.float64)
    return sums / sizes[:, None], sizes, patch_dom


def welch_f(means: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Welch's heteroscedastic one-way ANOVA F per gene, with patches as units.

    Args:
        means: Units x genes values.
        groups: Group label per unit; every group needs >= 2 units.

    Returns:
        F per gene (NaN where a group has zero variance).
    """
    labs = np.unique(groups)
    k = labs.size
    n = np.array([np.sum(groups == g) for g in labs], dtype=np.float64)
    m = np.stack([means[groups == g].mean(0) for g in labs])
    v = np.stack([means[groups == g].var(0, ddof=1) for g in labs])
    with np.errstate(divide="ignore", invalid="ignore"):
        w = n[:, None] / v
        wsum = w.sum(0)
        mw = (w * m).sum(0) / wsum
        a = (w * (m - mw) ** 2).sum(0) / (k - 1)
        lam = ((1 - w / wsum) ** 2 / (n[:, None] - 1)).sum(0)
        b = 1 + 2 * (k - 2) / (k**2 - 1) * lam
        f: np.ndarray = a / b
    return f


def studentized_permutation_p(
    means: np.ndarray, groups: np.ndarray, *, n_perm: int, seed: int
) -> float:
    """Permutation p of sum over genes of Welch F, permuting niche labels among patches.

    Studentised statistics keep permutation tests valid when variances differ by group.

    Args:
        means: Patch x gene means.
        groups: Niche per patch.
        n_perm: Permutations.
        seed: Random seed.

    Returns:
        p-value, or NaN if a niche has fewer than 2 patches.
    """
    labs, counts = np.unique(groups, return_counts=True)
    if labs.size < 2 or counts.min() < 2:
        return float("nan")

    def stat(g: np.ndarray) -> float:
        f = welch_f(means, g)
        return float(np.nansum(f))

    obs = stat(groups)
    rng = np.random.default_rng(seed)
    null = np.array([stat(rng.permutation(groups)) for _ in range(n_perm)])
    return float((1 + np.sum(null >= obs)) / (1 + n_perm))


def pseudobulk_anova_any(means: np.ndarray, groups: np.ndarray, *, q: float = 0.05) -> bool | None:
    """Per-gene one-way ANOVA on patch means, BH across genes; True if any gene q < ``q``.

    Args:
        means: Patch x gene means.
        groups: Niche per patch.
        q: FDR level.

    Returns:
        The call, or None if a niche has fewer than 2 patches.
    """
    labs, counts = np.unique(groups, return_counts=True)
    if labs.size < 2 or counts.min() < 2:
        return None
    _, p = stats.f_oneway(*[means[groups == g] for g in labs], axis=0)
    p = np.nan_to_num(np.asarray(p, dtype=np.float64), nan=1.0)
    return bool((benjamini_hochberg(p) < q).any())


def cluster_robust_any(
    resid: np.ndarray, groups: np.ndarray, cluster: np.ndarray, *, q: float = 0.05
) -> bool | None:
    """Per-gene OLS on niche dummies with cluster-robust (CR1) SEs; BH; any gene q < ``q``.

    The Wald statistic for all niche contrasts is referred to F(k-1, G-1), with G clusters.

    Args:
        resid: Cells x genes residuals.
        groups: Niche per cell.
        cluster: Cluster (patch) per cell.
        q: FDR level.

    Returns:
        The call, or None with fewer than 3 clusters or fewer than 2 niches.
    """
    labs = np.unique(groups)
    clus = np.unique(cluster)
    n_g = clus.size
    k = labs.size
    if k < 2 or n_g < 3:
        return None
    n = resid.shape[0]
    x = np.column_stack([np.ones(n)] + [(groups == g).astype(np.float64) for g in labs[1:]])
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ resid
    e = resid - x @ beta
    p_dim = x.shape[1]
    corr = (n_g / (n_g - 1)) * ((n - 1) / (n - p_dim))
    cidx = np.searchsorted(clus, cluster)
    pvals = np.ones(resid.shape[1])
    for j in range(resid.shape[1]):
        scores = np.zeros((n_g, p_dim))
        np.add.at(scores, cidx, x * e[:, j : j + 1])
        meat = scores.T @ scores
        v = corr * xtx_inv @ meat @ xtx_inv
        b = beta[1:, j]
        vb = v[1:, 1:]
        try:
            w = float(b @ np.linalg.solve(vb, b))
        except np.linalg.LinAlgError:
            continue
        f = w / (k - 1)
        pvals[j] = float(stats.f.sf(f, k - 1, n_g - 1))
    return bool((benjamini_hochberg(pvals) < q).any())
