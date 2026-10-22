"""Readers for 10x Visium outputs as deposited on GEO (GSE193460, Perturb-map)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

POSITION_COLUMNS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]


@dataclass(frozen=True)
class CountMatrix:
    """A spots x features UMI count matrix with its labels.

    Attributes:
        counts: CSR matrix, one row per spot barcode, one column per feature.
        barcodes: Spot barcodes, aligned with rows.
        feature_ids: Feature identifiers (Ensembl IDs, or custom names for transgenes).
        feature_names: Feature symbols, aligned with columns.
    """

    counts: sparse.csr_matrix
    barcodes: list[str]
    feature_ids: list[str]
    feature_names: list[str]


def _decode(values: np.ndarray) -> list[str]:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def read_10x_h5(path: Path) -> CountMatrix:
    """Read a Space Ranger ``filtered_feature_bc_matrix.h5`` file.

    The file stores a features x barcodes CSC matrix; it is returned transposed so that
    rows are spots.

    Args:
        path: Path to the ``.h5`` file.

    Returns:
        The count matrix with spot and feature labels.
    """
    with h5py.File(path, "r") as fh:
        grp = fh["matrix"]
        shape = tuple(int(x) for x in grp["shape"][:])
        csc = sparse.csc_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape)
        barcodes = _decode(grp["barcodes"][:])
        feature_ids = _decode(grp["features"]["id"][:])
        feature_names = _decode(grp["features"]["name"][:])
    return CountMatrix(
        counts=sparse.csr_matrix(csc.T),
        barcodes=barcodes,
        feature_ids=feature_ids,
        feature_names=feature_names,
    )


def read_positions(path: Path) -> pd.DataFrame:
    """Read a header-less Space Ranger ``tissue_positions_list.csv(.gz)``.

    Args:
        path: Path to the positions file.

    Returns:
        One row per barcode with array and pixel coordinates.
    """
    return pd.read_csv(path, header=None, names=POSITION_COLUMNS)


def read_spot_annotation(path: Path) -> pd.DataFrame:
    """Read the authors' per-spot annotation (``spot_annotation.csv.gz``).

    Missing labels are written as the string ``NA`` by R; they are kept as NaN here.

    Args:
        path: Path to the annotation file.

    Returns:
        One row per on-tissue spot with ``kmeans``, ``leiden_clusters`` and ``phenotypes``.
    """
    return pd.read_csv(path, na_values=["NA"], keep_default_na=False)
