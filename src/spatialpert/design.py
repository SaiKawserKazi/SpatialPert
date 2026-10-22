"""Screen design: design-effect power model and plan_screen() (review A7).

A knockout with L lesions of m cells each carries n_eff = L*m / (1 + (m - 1) * rho) effectively
independent cells, where rho is the within-lesion intraclass correlation of the tested score.
Power for a niche-specific effect is approximated by a two-sided z-test between two niches that
each receive half of the effective cells, scaled by an efficiency factor for the actual test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


def icc_oneway(values: np.ndarray, groups: np.ndarray) -> float:
    """One-way random-effects ICC(1) from a balanced or unbalanced grouping.

    Args:
        values: One score per cell.
        groups: Group (lesion) id per cell.

    Returns:
        Intraclass correlation, clipped to [0, 1].
    """
    labs, inv, counts = np.unique(groups, return_inverse=True, return_counts=True)
    k = labs.size
    n = values.size
    means = np.bincount(inv, weights=values) / counts
    grand = values.mean()
    ssb = float(np.sum(counts * (means - grand) ** 2))
    ssw = float(np.sum((values - means[inv]) ** 2))
    msb, msw = ssb / (k - 1), ssw / (n - k)
    n0 = (n - np.sum(counts**2) / n) / (k - 1)
    icc = (msb - msw) / (msb + (n0 - 1) * msw)
    return float(np.clip(icc, 0.0, 1.0))


def effective_n(lesions: int, cells_per_lesion: float, rho: float) -> float:
    """Design-effect effective sample size L*m / (1 + (m - 1) * rho).

    Args:
        lesions: Number of independent lesions.
        cells_per_lesion: Cells per lesion.
        rho: Within-lesion intraclass correlation.

    Returns:
        Effective number of independent cells.
    """
    return lesions * cells_per_lesion / (1 + (cells_per_lesion - 1) * rho)


@dataclass(frozen=True)
class DesignModel:
    """Fitted design-effect model.

    Attributes:
        rho: Within-lesion ICC of the tested score.
        sigma: Per-cell SD of the tested score (residual SD units).
        efficiency: Fitted scale (actual test vs ideal z-test), 1.0 if not fitted.
    """

    rho: float
    sigma: float
    efficiency: float = 1.0


def plan_screen_empirical(
    lesions: float | np.ndarray, scenario: str = "log_scale_lambda0.5"
) -> np.ndarray:
    """Power of the lesion-level test, interpolated from the simulated benchmark (400 cells).

    The design-effect model above did not reproduce simulated power (benchmark RMSE 0.27): power
    is limited by the number of independent lesions compared across niches, not by within-lesion
    correlation. This function therefore interpolates (in log lesions) the simulated curves.

    Args:
        lesions: Lesions per knockout.
        scenario: ``"log_scale_lambda0.5"`` (niche effect half the shared effect, log scale),
            ``"count_level_delta0.5"`` or ``"count_level_delta1.0"`` (realistic count-level
            effects, real lesion sizes).

    Returns:
        Interpolated power (0 below 3 lesions).
    """
    from spatialpert.design_table import POWER_TABLE  # noqa: PLC0415

    if scenario not in POWER_TABLE or scenario == "lesions":
        raise ValueError(f"unknown scenario {scenario!r}")
    les = np.asarray(lesions, dtype=np.float64)
    xs = np.log(np.asarray(POWER_TABLE["lesions"], dtype=np.float64))
    ys = np.asarray(POWER_TABLE[scenario], dtype=np.float64)
    out: np.ndarray = np.interp(np.log(np.clip(les, 1, None)), xs, ys)
    return np.where(les >= 3, out, 0.0)


def plan_screen(
    lesions: int | np.ndarray,
    cells_per_lesion: float | np.ndarray,
    effect_sd: float,
    model: DesignModel,
    *,
    alpha: float = 0.05,
) -> np.ndarray:
    """Predicted power to detect a niche-specific effect of ``effect_sd`` (residual SD units).

    Args:
        lesions: Lesions per knockout (scalar or array).
        cells_per_lesion: Cells per lesion (scalar or array).
        effect_sd: Niche-specific effect on the tested score, in residual SD units.
        model: Fitted :class:`DesignModel`.
        alpha: Two-sided significance level.

    Returns:
        Predicted power (array broadcast over the inputs).
    """
    les = np.asarray(lesions, dtype=np.float64)
    m = np.asarray(cells_per_lesion, dtype=np.float64)
    neff = les * m / (1 + (m - 1) * model.rho)
    se = model.sigma * np.sqrt(2.0 / (neff / 2.0))  # difference of two niche means
    z = model.efficiency * effect_sd / se
    zc = stats.norm.ppf(1 - alpha / 2)
    out: np.ndarray = stats.norm.cdf(z - zc) + stats.norm.cdf(-z - zc)
    # A lesion-level test needs at least two lesions per niche (four in total).
    return np.where(les >= 4, out, 0.0)
