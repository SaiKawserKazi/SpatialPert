"""Xenium cell-feature matrices with custom sgRNA probes (GSE319123, Perturb-DBiT validation).

The custom probes are stored as 'Gene Expression' features whose IDs are not Ensembl IDs,
named ``<TARGET>-<guide number>`` (e.g. ``PRSS23-2``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from spatialpert.visium import CountMatrix, read_10x_h5

NEGATIVE_PROBE = "Negative Control Probe"


@dataclass(frozen=True)
class XeniumMatrix:
    """A Xenium count matrix plus feature types.

    Attributes:
        matrix: Cells x features counts with labels.
        feature_types: Feature type per column (``Gene Expression``, control codewords, ...).
    """

    matrix: CountMatrix
    feature_types: list[str]


def read_xenium_h5(path: Path) -> XeniumMatrix:
    """Read a Xenium ``cell_feature_matrix.h5`` with its feature types.

    Args:
        path: Path to the ``.h5`` file.

    Returns:
        The matrix and per-feature types.
    """
    with h5py.File(path, "r") as fh:
        types = [t.decode() for t in fh["matrix"]["features"]["feature_type"][:]]
    return XeniumMatrix(matrix=read_10x_h5(path), feature_types=types)


def guide_probe_indices(x: XeniumMatrix) -> list[int]:
    """Columns holding custom sgRNA probes (expression features without an Ensembl ID).

    Args:
        x: The Xenium matrix.

    Returns:
        Column indices of the guide probes.
    """
    m = x.matrix
    return [
        i
        for i, (fid, ftype) in enumerate(zip(m.feature_ids, x.feature_types, strict=True))
        if ftype == "Gene Expression" and not fid.startswith("ENS")
    ]


def probe_target(name: str) -> str:
    """Target gene of a guide probe name such as ``PRSS23-2``.

    Args:
        name: Probe name.

    Returns:
        The gene symbol before the final ``-<number>``.
    """
    head, _, tail = name.rpartition("-")
    return head if head and tail.isdigit() else name


def negative_probe_rate(x: XeniumMatrix) -> float:
    """Per-cell probability that one negative-control probe shows >= 1 count.

    This estimates how often a single guide probe is 'detected' in a cell by chance.

    Args:
        x: The Xenium matrix.

    Returns:
        Mean fraction of cells positive across negative-control probes.
    """
    cols = [i for i, t in enumerate(x.feature_types) if t == NEGATIVE_PROBE]
    if not cols:
        return float("nan")
    sub = x.matrix.counts[:, cols]
    return float(np.asarray((sub > 0).mean(axis=0)).mean())


def guide_cell_table(x: XeniumMatrix, *, min_umis: int = 1) -> pd.DataFrame:
    """Count guide-positive cells per probe, with the chance expectation beside it.

    Args:
        x: The Xenium matrix.
        min_umis: UMIs needed to call a cell positive.

    Returns:
        One row per probe: ``probe``, ``target``, ``cells``, ``expected_by_chance``,
        ``enrichment`` (cells / expected).
    """
    m = x.matrix
    n_cells = m.counts.shape[0]
    rate = negative_probe_rate(x)
    csc = m.counts.tocsc()
    rows = []
    for i in guide_probe_indices(x):
        cells = int((csc[:, i] >= min_umis).sum())
        expected = rate * n_cells if min_umis == 1 else float("nan")
        rows.append(
            {
                "probe": m.feature_names[i],
                "target": probe_target(m.feature_names[i]),
                "cells": cells,
                "expected_by_chance": expected,
                "enrichment": cells / expected if expected and expected > 0 else float("nan"),
            }
        )
    return pd.DataFrame(rows).sort_values("cells", ascending=False, ignore_index=True)
