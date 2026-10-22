"""Domain-specific knockout effects and the location-dependence tests (protocol §4-5)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class ControlModel:
    """Per-gene control fit: domain means plus covariate slopes (protocol §4).

    Attributes:
        domain_means: Domains x genes, control mean at the covariate centre.
        beta: Covariates x genes slopes.
        centre: Covariate centring values.
        domains: Domain ids, ordered as rows of ``domain_means``.
    """

    domain_means: np.ndarray
    beta: np.ndarray
    centre: np.ndarray
    domains: np.ndarray


def _design(
    dom: np.ndarray, cov: np.ndarray, domains: np.ndarray, centre: np.ndarray
) -> np.ndarray:
    onehot = (dom[:, None] == domains[None, :]).astype(np.float64)
    return np.hstack([onehot, cov - centre])


def fit_control_model(y: sparse.csr_matrix, dom: np.ndarray, cov: np.ndarray) -> ControlModel:
    """Least-squares fit of y = mu_domain + beta' (c - centre) on control cells.

    Uses the normal equations so that ``y`` stays sparse.

    Args:
        y: Control cells x genes (log-normalised).
        dom: Domain per control cell.
        cov: Control cells x covariates (no NaN).

    Returns:
        The fitted control model.
    """
    domains = np.unique(dom)
    centre = cov.mean(axis=0)
    x = _design(dom, cov, domains, centre)
    xtx = x.T @ x
    xty = np.asarray((y.T @ x).T)  # p x genes
    coef = np.linalg.solve(xtx, xty)
    k = domains.size
    return ControlModel(domain_means=coef[:k], beta=coef[k:], centre=centre, domains=domains)


def residuals(
    model: ControlModel, y: sparse.csr_matrix | np.ndarray, dom: np.ndarray, cov: np.ndarray
) -> np.ndarray:
    """Residuals y - mu_domain - beta'(c - centre) for any cells (dense output).

    Args:
        model: Fitted control model.
        y: Cells x genes.
        dom: Domain per cell (must be in ``model.domains``).
        cov: Cells x covariates.

    Returns:
        Cells x genes residuals.
    """
    dense = np.asarray(y.toarray() if isinstance(y, sparse.csr_matrix) else y, dtype=np.float64)
    pos = np.searchsorted(model.domains, dom)
    pred = model.domain_means[pos] + (cov - model.centre) @ model.beta
    out: np.ndarray = dense - pred
    return out


def covariate_adjusted(
    model: ControlModel, y: sparse.csr_matrix | np.ndarray, cov: np.ndarray
) -> np.ndarray:
    """Covariate-adjusted expression z = y - beta'(c - centre), dense.

    Args:
        model: Fitted control model (only its slopes and centre are used).
        y: Cells x genes.
        cov: Cells x covariates.

    Returns:
        Cells x genes adjusted expression.
    """
    dense = np.asarray(y.toarray() if isinstance(y, sparse.csr_matrix) else y, dtype=np.float64)
    out: np.ndarray = dense - (cov - model.centre) @ model.beta
    return out


def local_residuals(  # noqa: PLR0913 - keyword-only options
    z_ko: np.ndarray,
    ko_coords: np.ndarray,
    z_ctrl: np.ndarray,
    ctrl_coords: np.ndarray,
    *,
    exclude: np.ndarray | None = None,
    k: int = 50,
) -> np.ndarray:
    """Local-control residuals: each cell minus the mean of its k nearest control cells.

    Removes lesion- and region-specific baselines that tissue-wide niche means leave behind
    (review AMENDMENT-005).

    Args:
        z_ko: Knockout cells x genes, covariate-adjusted.
        ko_coords: Knockout cells x 2.
        z_ctrl: Control cells x genes, covariate-adjusted.
        ctrl_coords: Control cells x 2.
        exclude: Optional boolean mask over control cells that must not serve as controls
            (for pseudo-knockouts built from controls: the pseudo-knockout cells themselves).
        k: Number of nearest controls.

    Returns:
        Knockout cells x genes residuals.
    """
    pool = np.arange(ctrl_coords.shape[0]) if exclude is None else np.where(~exclude)[0]
    tree = cKDTree(ctrl_coords[pool])
    _, nn = tree.query(ko_coords, k=min(k, pool.size))
    nn = np.atleast_2d(nn)
    base = z_ctrl[pool[nn]].mean(axis=1)
    out: np.ndarray = z_ko - base
    return out


def clone_patches(coords: np.ndarray, *, link: float = 30.0) -> np.ndarray:
    """Single-linkage patches of cells closer than ``link`` µm (protocol §4).

    Args:
        coords: Cells x 2.
        link: Linkage distance.

    Returns:
        Patch id per cell.
    """
    n = coords.shape[0]
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = cKDTree(coords).query_pairs(link, output_type="ndarray")
    graph = sparse.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, lab = connected_components(graph, directed=False)
    return np.asarray(lab)


def patch_domains(patch: np.ndarray, dom: np.ndarray) -> np.ndarray:
    """Majority domain of each patch (ties go to the smallest domain id).

    Args:
        patch: Patch id per cell (0..P-1).
        dom: Domain per cell.

    Returns:
        Domain per patch.
    """
    n_p = int(patch.max()) + 1
    present = np.bincount(patch, minlength=n_p) > 0
    if not present.all():
        raise ValueError(
            "patch ids must be contiguous 0..P-1; re-index with np.unique(ids, return_inverse=True)"
        )
    out = np.empty(n_p, dtype=dom.dtype)
    for p in range(n_p):
        vals, counts = np.unique(dom[patch == p], return_counts=True)
        out[p] = vals[np.argmax(counts)]
    return out


def drop_small_niches(
    patch: np.ndarray, patch_dom: np.ndarray, min_cells: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Keep only patches whose (majority) niche holds >= ``min_cells`` cells *after* patching.

    Applying the minimum before patch-majority assignment lets a niche shrink to a handful of
    cells, which can then dominate the n-weighted statistic (found in review: DIO3 niche 3 held
    7 cells after assignment and drove an apparent signal).

    Args:
        patch: Patch id per cell (0..P-1).
        patch_dom: Niche per patch.
        min_cells: Minimum cells per niche after assignment.

    Returns:
        (mask over the input cells, re-indexed patch ids of kept cells, niche per kept patch).
    """
    sizes = np.bincount(patch, minlength=patch_dom.size)
    keep_patch = np.ones(patch_dom.size, dtype=bool)
    while True:
        bad = [
            d
            for d in np.unique(patch_dom[keep_patch])
            if sizes[keep_patch & (patch_dom == d)].sum() < min_cells
        ]
        if not bad:
            break
        keep_patch &= ~np.isin(patch_dom, bad)
    mask: np.ndarray = keep_patch[patch]
    new_id = np.cumsum(keep_patch) - 1
    return mask, new_id[patch[mask]], patch_dom[keep_patch]


def between_domain_ss(sums: np.ndarray, sizes: np.ndarray) -> float:
    """T = sum_d n_d ||e_d - e_bar||^2 from per-domain residual sums (protocol H1).

    Args:
        sums: Domains x genes residual sums.
        sizes: Cells per domain.

    Returns:
        The between-domain sum of squares, summed over genes.
    """
    keep = sizes > 0
    s, n = sums[keep], sizes[keep]
    total = s.sum(axis=0)
    return float((s**2 / n[:, None]).sum() - (total**2).sum() / n.sum())


@dataclass(frozen=True)
class InteractionTest:
    """Result of the patch-permutation test.

    Attributes:
        statistic: Observed T.
        p_value: (1 + #perm >= obs) / (1 + n_perm).
        n_patches: Patches used.
        n_cells: Cells used.
        effects: Domains x genes mean residuals (observed, patch-majority domains).
        domains: Domain ids of ``effects`` rows.
        sizes: Cells per domain.
        null: Permutation statistics.
    """

    statistic: float
    p_value: float
    n_patches: int
    n_cells: int
    effects: np.ndarray
    domains: np.ndarray
    sizes: np.ndarray
    null: np.ndarray


def patch_permutation_test(
    resid: np.ndarray,
    patch: np.ndarray,
    patch_dom: np.ndarray,
    *,
    n_perm: int = 1000,
    seed: int = 0,
) -> InteractionTest:
    """H1 test: permute domain labels among patches, residuals fixed.

    Args:
        resid: Cells x genes residuals of one knockout's cells.
        patch: Patch id per cell (0..P-1).
        patch_dom: Domain per patch.
        n_perm: Permutations.
        seed: Random seed.

    Returns:
        The test result.
    """
    n_p = patch_dom.size
    psum = np.zeros((n_p, resid.shape[1]))
    np.add.at(psum, patch, resid)
    psize = np.bincount(patch, minlength=n_p).astype(np.float64)
    domains = np.unique(patch_dom)
    onehot = (patch_dom[:, None] == domains[None, :]).astype(np.float64)

    def stat(assign: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        sums = assign.T @ psum
        sizes = assign.T @ psize
        return between_domain_ss(sums, sizes), sums, sizes

    obs, sums, sizes = stat(onehot)
    rng = np.random.default_rng(seed)
    null = np.array([stat(onehot[rng.permutation(n_p)])[0] for _ in range(n_perm)])
    exceed = int(np.sum(null >= obs))
    with np.errstate(invalid="ignore", divide="ignore"):
        effects = sums / sizes[:, None]
    return InteractionTest(
        null=null,
        statistic=obs,
        p_value=(1 + exceed) / (1 + n_perm),
        n_patches=n_p,
        n_cells=int(psize.sum()),
        effects=effects,
        domains=domains,
        sizes=sizes,
    )


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values.

    Args:
        p: Raw p-values.

    Returns:
        Adjusted p-values (same order).
    """
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def leave_one_domain_out_r(effects: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """H2: Pearson r between each domain's effect and the pooled effect of the others.

    Args:
        effects: Domains x genes mean residuals.
        sizes: Cells per domain.

    Returns:
        One r per domain (NaN if undefined).
    """
    out = np.full(effects.shape[0], np.nan)
    for d in range(effects.shape[0]):
        others = np.arange(effects.shape[0]) != d
        w = sizes[others]
        if w.sum() == 0 or sizes[d] == 0:
            continue
        pooled = (effects[others] * w[:, None]).sum(axis=0) / w.sum()
        a, b = effects[d], pooled
        if a.std() == 0 or b.std() == 0:
            continue
        out[d] = float(np.corrcoef(a, b)[0, 1])
    return out


def matched_pseudo_patches(
    ctrl_coords: np.ndarray,
    ctrl_dom: np.ndarray,
    sizes: np.ndarray,
    doms: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw control-cell patches matching a knockout's patch sizes and domains (PC2).

    Each pseudo patch is a random control seed in the required domain plus its nearest
    control neighbours. Cells may not be reused within one pseudo-knockout.

    Args:
        ctrl_coords: Control cells x 2.
        ctrl_dom: Domain per control cell.
        sizes: Target patch sizes.
        doms: Target patch domains.
        rng: Random generator.

    Returns:
        (control-cell indices, patch id per selected cell).
    """
    tree = cKDTree(ctrl_coords)
    used = np.zeros(ctrl_coords.shape[0], dtype=bool)
    idx_out: list[np.ndarray] = []
    patch_out: list[np.ndarray] = []
    for p, (size, d) in enumerate(zip(sizes, doms, strict=True)):
        pool = np.where((ctrl_dom == d) & ~used)[0]
        if pool.size == 0:
            continue
        seed_cell = rng.choice(pool)
        k = min(int(size) * 3 + 5, ctrl_coords.shape[0])
        _, nn = tree.query(ctrl_coords[seed_cell], k=k)
        nn = np.atleast_1d(nn)
        nn = nn[~used[nn]][: int(size)]
        used[nn] = True
        idx_out.append(nn)
        patch_out.append(np.full(nn.size, p))
    if not idx_out:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    idx = np.concatenate(idx_out)
    pid = np.concatenate(patch_out)
    _, pid = np.unique(pid, return_inverse=True)
    return idx, pid
