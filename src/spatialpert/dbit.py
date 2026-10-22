"""Readers for Perturb-DBiT files as deposited on GEO (GSE319277).

Pixel identifiers are ``"<col>x<row>"`` strings. The sgRNA files prefix them with ``bc``
(``bc10x13``); expression and position files do not.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

_SG_FIRST = re.compile(r"^sg\d+_(?P<gene>[A-Za-z0-9.\-]+)$")  # sg57430_ATG4A (Brunello/Brie)
_GENE_FIRST = re.compile(r"^(?P<gene>[A-Za-z0-9.\-]+)_sg\d+$")  # Tgfbr2_sg2 (mTSG)
_CONTROL = re.compile(r"(?i)^(ntc|non[-_ ]?targeting|nt|control|safe)")


def normalise_pixel(pixel: str) -> str:
    """Strip the ``bc`` prefix so sgRNA and expression pixels share one key.

    Args:
        pixel: A pixel id such as ``bc10x13`` or ``10x13``.

    Returns:
        The id without prefix, e.g. ``10x13``.
    """
    return pixel[2:] if pixel.startswith("bc") else pixel


def guide_target(guide: str) -> str | None:
    """Return the gene a guide targets, ``"NTC"`` for controls, or None if unparseable.

    Args:
        guide: A guide name such as ``sg57430_ATG4A`` or ``Tgfbr2_sg2``.

    Returns:
        The target gene, ``"NTC"``, or None.
    """
    if _CONTROL.match(guide):
        return "NTC"
    for pattern in (_SG_FIRST, _GENE_FIRST):
        match = pattern.match(guide)
        if match is not None:
            return match.group("gene")
    return None


def read_sg_counts(path: Path) -> pd.DataFrame:
    """Read a header-less sgRNA count file (pixel, guide, count).

    Args:
        path: Path to ``sgCounts_*.txt(.gz)``.

    Returns:
        Columns ``pixel``, ``guide``, ``count``, ``gene``.
    """
    df = pd.read_csv(
        path, sep="\t", header=None, usecols=[0, 1, 2], names=["pixel", "guide", "count"]
    )
    df["pixel"] = df["pixel"].astype(str).map(normalise_pixel)
    df["gene"] = df["guide"].astype(str).map(guide_target)
    return df


def read_position_pixels(path: Path) -> list[str]:
    """Read the in-tissue pixel list (one CSV line, leading empty cell).

    Args:
        path: Path to ``position_*.txt(.gz)``.

    Returns:
        Pixel ids on tissue.
    """
    text = pd.read_csv(path, header=None, dtype=str).iloc[0].dropna()
    return [normalise_pixel(p) for p in text.tolist() if p]


_PIXEL = re.compile(r"^\d+x\d+$")


def _gene_of(feature: str) -> tuple[str, str]:
    """Split ASTRO-style ``GENE__exon`` features; plain names have an empty suffix."""
    head, sep, tail = feature.rpartition("__")
    return (head, tail) if sep else (feature, "")


def read_dbit_expression(path: Path, genes: list[str]) -> pd.DataFrame:
    """Read only the requested genes from a DBiT matrix, in either orientation.

    Pixels-as-rows files are read with ``usecols``; genes-as-rows files are streamed in
    chunks. For ASTRO feature names (``GENE__exon``, ``GENE__transcript``, ...) the ``exon``
    feature is used.

    Args:
        path: Matrix file (.tsv or .tsv.gz).
        genes: Gene symbols wanted.

    Returns:
        Pixels x genes counts (only genes found), pixel ids without prefix.
    """
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    wanted = set(genes)
    if sum(bool(_PIXEL.match(str(c).strip('"'))) for c in header[:50]) > 25:
        chunks = []
        for chunk in pd.read_csv(path, sep="\t", index_col=0, chunksize=2000):
            chunk.index = chunk.index.astype(str)
            chunks.append(chunk[chunk.index.isin(wanted)])
        df = pd.concat(chunks).T
    else:
        pick: dict[str, str] = {}
        for col in header[1:]:
            gene, suffix = _gene_of(str(col))
            if gene in wanted and suffix in ("", "exon"):
                pick.setdefault(gene, str(col))
        df = pd.read_csv(path, sep="\t", index_col=0, usecols=[header[0], *pick.values()])
        df = df.rename(columns={v: k for k, v in pick.items()})
    df.index = pd.Index([normalise_pixel(str(i)) for i in df.index])
    out: pd.DataFrame = df.astype(np.float64)
    return out


def pixel_knockouts(sg: pd.DataFrame, *, min_umis: int = 1) -> pd.DataFrame:
    """Assign each pixel the set of knockout genes detected in it.

    A pixel can hold several clones; every gene with at least ``min_umis`` guide UMIs counts.
    Unparseable guides are dropped.

    Args:
        sg: Output of :func:`read_sg_counts`.
        min_umis: Minimum guide UMIs per (pixel, gene).

    Returns:
        Columns ``pixel``, ``gene``, ``umis``: one row per detected (pixel, gene).
    """
    known = sg.dropna(subset=["gene"])
    agg: pd.DataFrame = known.groupby(["pixel", "gene"])["count"].sum().rename("umis").reset_index()
    kept: pd.DataFrame = agg[agg["umis"] >= min_umis].reset_index(drop=True)
    return kept
