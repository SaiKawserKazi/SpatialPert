"""Aim 6 semi-synthetic benchmark: known location-dependent effects in real control cells."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import f_oneway

from spatialpert.effects import benjamini_hochberg, patch_permutation_test


def unit_rms(v: np.ndarray) -> np.ndarray:
    """Scale a vector so that mean(v**2) == 1.

    Args:
        v: Vector.

    Returns:
        The rescaled vector.
    """
    out: np.ndarray = v / np.sqrt(np.mean(v**2))
    return out


def inject_effect(  # noqa: PLR0913 - keyword-only effect parameters
    resid: np.ndarray,
    dom: np.ndarray,
    gene_sd: np.ndarray,
    *,
    delta: float,
    lam: float,
    target_domain: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Add core + domain-specific effects to residuals.

    Every cell gets ``delta * sd * u``; cells in ``target_domain`` also get
    ``lam * delta * sd * v``. ``u`` and ``v`` are independent random unit-RMS directions.

    Args:
        resid: Cells x genes residuals of real control cells.
        dom: Domain per cell.
        gene_sd: Per-gene residual SD.
        delta: Core effect size in SD units.
        lam: Domain-specific part relative to the core.
        target_domain: Domain receiving the extra part.
        rng: Random generator.

    Returns:
        Residuals with the effect added.
    """
    g = resid.shape[1]
    u = unit_rms(np.asarray(rng.normal(size=g), dtype=np.float64))
    v = unit_rms(np.asarray(rng.normal(size=g), dtype=np.float64))
    out: np.ndarray = resid + delta * gene_sd * u
    out[dom == target_domain] += lam * delta * gene_sd * v
    return out


def cell_permutation_p(resid: np.ndarray, dom: np.ndarray, *, n_perm: int, seed: int) -> float:
    """Naive test: permute domain labels across cells (ignores clonal patches).

    Args:
        resid: Cells x genes residuals.
        dom: Domain per cell.
        n_perm: Permutations.
        seed: Random seed.

    Returns:
        Permutation p-value.
    """
    patch = np.arange(resid.shape[0])
    return patch_permutation_test(resid, patch, dom, n_perm=n_perm, seed=seed).p_value


def anova_any_gene(resid: np.ndarray, dom: np.ndarray, *, q: float = 0.05) -> bool:
    """Per-gene one-way ANOVA across domains, BH across genes; True if any gene q < ``q``.

    Args:
        resid: Cells x genes residuals.
        dom: Domain per cell.
        q: FDR level.

    Returns:
        Whether location dependence is called.
    """
    groups = [resid[dom == d] for d in np.unique(dom)]
    if len(groups) < 2:
        return False
    _, p = f_oneway(*groups, axis=0)
    p = np.nan_to_num(np.asarray(p, dtype=np.float64), nan=1.0)
    return bool((benjamini_hochberg(p) < q).any())


def scale_sizes(sizes: np.ndarray, n_total: int) -> np.ndarray:
    """Rescale patch sizes to sum to about ``n_total``, keeping every patch >= 1.

    Args:
        sizes: Template patch sizes.
        n_total: Desired total cells.

    Returns:
        New integer sizes.
    """
    scaled: np.ndarray = np.maximum(1, np.round(sizes * n_total / sizes.sum())).astype(int)
    return scaled


@dataclass(frozen=True)
class Calls:
    """Detection calls of the three tests for one replicate."""

    patch: bool
    cell: bool
    anova: bool
