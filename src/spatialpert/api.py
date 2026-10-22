"""One-call public API: test whether knockout effects depend on spatial domain.

Works on any platform's cell (or spot) table: expression, coordinates, guide labels, domain
labels and optional covariates.

Example:
    >>> from spatialpert.api import location_dependence
    >>> res = location_dependence(y, coords, labels, domains, covariates=cov)
    >>> res.table[["knockout", "T", "p", "q"]]
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import sparse

from spatialpert.baselines import patch_means, studentized_permutation_p
from spatialpert.effects import (
    benjamini_hochberg,
    clone_patches,
    drop_small_niches,
    fit_control_model,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.pipeline import tile_ids


def from_anndata(  # noqa: PLR0913 - mirrors the AnnData field names users pass
    adata: object,
    *,
    label_key: str,
    domain_key: str,
    spatial_key: str = "spatial",
    covariate_keys: list[str] | None = None,
    layer: str | None = None,
    **kwargs: object,
) -> LocationResult:
    """Run :func:`location_dependence` on an AnnData-like object.

    Duck-typed: needs ``.X`` (or ``.layers[layer]``), ``.obs`` (DataFrame), ``.obsm`` and
    ``.var_names``, so ``anndata`` is not a hard dependency. ``X`` must be log-normalised.

    Args:
        adata: AnnData (or compatible) object.
        label_key: ``obs`` column with the knockout label (controls labelled ``"control"``).
        domain_key: ``obs`` column with integer niche/domain ids.
        spatial_key: ``obsm`` key with x-y coordinates in µm.
        covariate_keys: Optional ``obs`` columns used as covariates.
        layer: Optional layer to use instead of ``X``.
        **kwargs: Passed to :func:`location_dependence`.

    Returns:
        The location-dependence result.
    """
    obs: pd.DataFrame = adata.obs  # type: ignore[attr-defined]
    y = adata.layers[layer] if layer else adata.X  # type: ignore[attr-defined]
    coords = np.asarray(adata.obsm[spatial_key], dtype=np.float64)[:, :2]  # type: ignore[attr-defined]
    cov = obs[covariate_keys].to_numpy(dtype=np.float64) if covariate_keys else None
    return location_dependence(
        y, coords, obs[label_key].to_numpy(), obs[domain_key].to_numpy(),
        covariates=cov, genes=[str(g) for g in adata.var_names], **kwargs,  # type: ignore[attr-defined, arg-type]
    )


@dataclass(frozen=True)
class LocationResult:
    """Output of :func:`location_dependence`.

    Attributes:
        table: One row per tested knockout: units, cells, T, p, q, pooled-test fields.
        effects: knockout -> DataFrame (domains x genes) of mean residual effects.
        pooled_p: p-value of the pooled statistic sum_k T_k (NaN if nothing tested).
    """

    table: pd.DataFrame
    effects: dict[str, pd.DataFrame] = field(default_factory=dict)
    pooled_p: float = float("nan")


def _units(xy: np.ndarray, unit: str, size: dict[str, float]) -> np.ndarray:
    """Permutation units: single-linkage patches or lesions, or square tiles."""
    if unit == "tile":
        return tile_ids(xy, size["tile"])
    return clone_patches(xy, link=size[unit])


def _pooled(pv: np.ndarray, statistic: str, obs: float, null: np.ndarray) -> float:
    """Pooled p: summed-statistic permutation (ss) or Cauchy combination (studentized)."""
    if pv.size == 0:
        return float("nan")
    if statistic == "ss":
        return float((1 + np.sum(null >= obs)) / (1 + null.size))
    return float(0.5 - np.arctan(np.mean(np.tan((0.5 - pv) * np.pi))) / np.pi)


def location_dependence(  # noqa: PLR0913 - public API with documented options
    y: sparse.spmatrix | np.ndarray,
    coords: np.ndarray,
    labels: np.ndarray,
    domains: np.ndarray,
    *,
    covariates: np.ndarray | None = None,
    genes: list[str] | None = None,
    control_label: str = "control",
    knockouts: list[str] | None = None,
    min_cells_per_domain: int = 30,
    unit: str = "lesion",
    statistic: str = "ss",
    link_um: float = 30.0,
    lesion_um: float = 100.0,
    tile_um: float = 250.0,
    n_perm: int = 1000,
    seed: int = 0,
) -> LocationResult:
    """Test each knockout for domain-dependent effects with a unit-level permutation.

    Residuals are taken against domain-specific control means plus covariate slopes. The
    statistic is the between-domain sum of squares of a knockout's residuals. Domain labels
    are permuted among clone patches (``unit="patch"``) or square tiles (``unit="tile"``),
    so spatially clustered clones do not inflate significance.

    Args:
        y: Cells x genes, log-normalised.
        coords: Cells x 2 positions (µm).
        labels: Per-cell knockout label; ``control_label`` marks guide-negative cells.
        domains: Per-cell integer domain (< 0 = unassigned, excluded).
        covariates: Optional cells x c covariates (no NaN).
        genes: Optional gene names for the effect tables.
        control_label: Label of control cells.
        knockouts: Knockouts to test (default: all others with enough cells).
        min_cells_per_domain: Cells needed in a domain, checked again after patch-majority
            assignment (a niche that shrinks below it is dropped).
        unit: ``"lesion"`` (default; single linkage at ``lesion_um``), ``"patch"`` (single
            linkage at ``link_um``) or ``"tile"``. Lesion units were the only unit calibrated at
            3-20 lesions in benchmarks, provided niches are defined from control cells only.
        statistic: ``"ss"`` (between-niche sum of squares) or ``"studentized"`` (sum of
            per-gene Welch F on unit means; robust to unequal variance, needs >= 2 units
            per niche).
        link_um: Patch linkage distance.
        lesion_um: Lesion linkage distance.
        tile_um: Tile edge.
        n_perm: Permutations.
        seed: Random seed.

    Returns:
        Per-knockout results, effect tables and the pooled p-value.

    Raises:
        ValueError: On an unknown ``unit`` or mismatched input lengths.
    """
    if unit not in ("patch", "lesion", "tile"):
        raise ValueError("unit must be 'patch', 'lesion' or 'tile'")
    if statistic not in ("ss", "studentized"):
        raise ValueError("statistic must be 'ss' or 'studentized'")
    n = coords.shape[0]
    if not (len(labels) == len(domains) == y.shape[0] == n):
        raise ValueError("y, coords, labels and domains must have the same number of cells")
    ycsr = sparse.csr_matrix(y)
    cov = np.zeros((n, 0)) if covariates is None else np.asarray(covariates, dtype=np.float64)
    labels = np.asarray(labels).astype(str)
    domains = np.asarray(domains)
    names = genes if genes is not None else [f"g{i}" for i in range(ycsr.shape[1])]

    ok = domains >= 0
    ctrl = ok & (labels == control_label)
    model = fit_control_model(ycsr[ctrl], domains[ctrl], cov[ctrl])
    cand = knockouts or sorted(set(labels[ok]) - {control_label})

    rows, effects = [], {}
    pooled_obs, pooled_null = 0.0, np.zeros(n_perm)
    for i, ko in enumerate(cand):
        mk = ok & (labels == ko)
        counts = pd.Series(domains[mk]).value_counts()
        use = counts.index[counts >= min_cells_per_domain].to_numpy()
        sel = np.where(mk & np.isin(domains, use))[0]
        if use.size < 2:
            continue
        uid = _units(coords[sel], unit, {"patch": link_um, "lesion": lesion_um, "tile": tile_um})
        udom = patch_domains(uid, domains[sel])
        keep, uid, udom = drop_small_niches(uid, udom, min_cells_per_domain)
        sel = sel[keep]
        if np.unique(udom).size < 2:
            continue
        r = residuals(model, ycsr[sel], domains[sel], cov[sel])
        res = patch_permutation_test(r, uid, udom, n_perm=n_perm, seed=seed + i)
        pval = res.p_value
        if statistic == "studentized":
            means, _, g = patch_means(r, uid, udom)
            pval = studentized_permutation_p(means, g, n_perm=n_perm, seed=seed + i)
            if np.isnan(pval):
                continue  # fewer than 2 units in some niche: untestable
        pooled_obs += res.statistic
        pooled_null += res.null
        rows.append(
            {"knockout": ko, "units": res.n_patches, "cells": res.n_cells,
             "domains": len(res.domains), "T": res.statistic, "p": pval}
        )
        effects[ko] = pd.DataFrame(res.effects, index=res.domains, columns=names)

    table = pd.DataFrame(rows, columns=["knockout", "units", "cells", "domains", "T", "p"])
    table["q"] = benjamini_hochberg(table["p"].to_numpy()) if len(table) else []
    pooled = _pooled(table["p"].to_numpy(), statistic, pooled_obs, pooled_null)
    return LocationResult(table=table, effects=effects, pooled_p=pooled)
