"""Tissue-specific calibration of niche-dependence tests (review 3, Priority 1).

Lesion-level permutation assumes that lesion-mean residuals are exchangeable across niches. In
tissues where lesion baselines are large, spatially correlated and differ in spread between
niches (benchmark: mouse brain), that assumption fails and every test is inflated. Two tools:

- :func:`lesion_exchangeability` measures, from control pseudo-lesions, how strongly lesion
  baselines dominate and how unequal their spread is across niches;
- :func:`matched_null_p` and :func:`calibrate_p` replace the permutation p-value by its rank
  among control pseudo-knockouts with the knockout's own lesion geometry in the same tissue.
  This is exact if knockout lesions are placed like random control positions (stated
  assumption, not testable from the data).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from spatialpert.api import location_dependence
from spatialpert.effects import clone_patches
from spatialpert.gmn import place_geometry

PSEUDO = "__pseudo_knockout__"


def calibrate_p(p_obs: float, null_p: np.ndarray, *, seed: int = 0) -> float:
    """Randomised rank of an observed p-value among matched null p-values.

    Args:
        p_obs: Observed p-value of the knockout.
        null_p: p-values of matched pseudo-knockouts (NaN = untestable, ignored).
        seed: Seed for breaking ties uniformly (keeps the size exact for discrete p-values).

    Returns:
        Calibrated p-value in (0, 1], or NaN if no pseudo-knockout was testable.
    """
    pool = np.sort(np.asarray(null_p, dtype=np.float64))
    pool = pool[~np.isnan(pool)]
    if pool.size == 0 or np.isnan(p_obs):
        return float("nan")
    below = int(np.searchsorted(pool, p_obs, side="left"))
    ties = int(np.searchsorted(pool, p_obs, side="right")) - below
    u = np.random.default_rng(seed).uniform()
    return float((below + u * (ties + 1)) / (1 + pool.size))


def lesion_exchangeability(
    resid: np.ndarray, lesion: np.ndarray, niche: np.ndarray, *, min_cells: int = 3
) -> dict[str, float]:
    """How far lesion-mean residuals are from exchangeable across niches.

    Args:
        resid: Cells x genes residuals of (pseudo-)knockout cells, each gene in its own SD units.
        lesion: Lesion id per cell.
        niche: Niche per cell (the lesion takes its majority niche).
        min_cells: Lesions with fewer cells are ignored.

    Returns:
        ``variance_ratio``: between-lesion variance of lesion means over the variance expected
        from sampling alone (1 = no lesion baselines); ``niche_sd_ratio``: largest over smallest
        SD of lesion means across niches with at least two lesions (1 = equal spread).
    """
    labs, inv = np.unique(lesion, return_inverse=True)
    n = np.bincount(inv)
    keep = n >= min_cells
    means = np.stack([resid[inv == k].mean(0) for k in range(labs.size)])
    nan_row = np.full(resid.shape[1], np.nan)
    within = np.stack(
        [resid[inv == k].var(0, ddof=1) / n[k] if n[k] > 1 else nan_row for k in range(labs.size)]
    )
    major = np.array([np.bincount(niche[inv == k]).argmax() for k in range(labs.size)])
    means, within, major = means[keep], within[keep], major[keep]
    ratio = float("nan")
    if means.shape[0] > 1:
        ratio = float(np.nanmean(means.var(0, ddof=1)) / np.nanmean(within))
    sds = [
        means[major == d].std(0, ddof=1).mean() for d in np.unique(major) if (major == d).sum() >= 2
    ]
    sd_ratio = float(max(sds) / max(min(sds), 1e-12)) if len(sds) >= 2 else float("nan")
    return {"variance_ratio": ratio, "niche_sd_ratio": sd_ratio}


def matched_null_p(  # noqa: PLR0913 - mirrors location_dependence
    y: sparse.spmatrix | np.ndarray,
    coords: np.ndarray,
    labels: np.ndarray,
    domains: np.ndarray,
    knockout: str,
    *,
    relabel: Callable[[np.ndarray], np.ndarray] | None = None,
    control_label: str = "control",
    n_null: int = 200,
    seed: int = 0,
    **test_kwargs: object,
) -> np.ndarray:
    """p-values of control pseudo-knockouts with the knockout's lesion sizes.

    Each pseudo-knockout picks well-separated control lesions with the same cell counts as the
    knockout's lesions, is removed from the control fit, and is tested by
    :func:`location_dependence` with the same options as the knockout.

    Args:
        y: Cells x genes, as in :func:`location_dependence`.
        coords: Cells x 2 positions (µm).
        labels: Per-cell knockout label.
        domains: Per-cell niche (< 0 = unassigned).
        knockout: Knockout whose lesion geometry is matched.
        relabel: Optional function giving niche labels of pseudo-knockout cells (indices into
            the cells) computed with those cells left out of the niche definition. Without it
            the original labels are used; these are circular and make the calibration
            conservative.
        control_label: Label of control cells.
        n_null: Pseudo-knockouts.
        seed: Random seed.
        **test_kwargs: Passed to :func:`location_dependence` (unit, statistic, n_perm, ...).

    Returns:
        One p-value per pseudo-knockout (NaN if untestable).
    """
    labels = np.asarray(labels).astype(str).astype(object)  # object: no fixed-width truncation
    domains = np.asarray(domains)
    lesion_um = float(test_kwargs.get("lesion_um", 100.0))  # type: ignore[arg-type]
    ko = np.where(labels == knockout)[0]
    les = np.unique(clone_patches(coords[ko], link=lesion_um), return_inverse=True)[1]
    sizes = np.bincount(les).tolist()
    ctrl = np.where((labels == control_label) & (domains >= 0))[0]
    tree = cKDTree(coords[ctrl])
    rng = np.random.default_rng(seed)
    out = np.full(n_null, np.nan)
    for m in range(n_null):
        pick = ctrl[np.asarray(place_geometry(tree, coords[ctrl], sizes, rng))]
        lab = labels.copy()
        lab[pick] = PSEUDO
        dom = domains.copy()
        if relabel is not None:
            dom[pick] = relabel(pick)
        res = location_dependence(
            y,
            coords,
            lab,
            dom,
            control_label=control_label,
            knockouts=[PSEUDO],
            seed=int(rng.integers(1 << 31)),
            **test_kwargs,  # type: ignore[arg-type]
        )
        if len(res.table):
            out[m] = float(res.table["p"].iloc[0])
    return out
