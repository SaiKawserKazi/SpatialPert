"""Week-1 data gate: can GSE193460 support per-knockout spatial tests?

The proposal's rule is that at least ``min_testable`` knockouts must have at least
``min_spots`` perturbed spots in at least ``min_units`` spatial domains. Domains are not
called at this stage, so the gate uses an upper bound that cannot be beaten by any domain
calling: a knockout needs ``min_spots * min_units`` labelled spots in total to be testable at
all. A stricter, replication-based count (``min_spots`` in ``min_units`` separate sections)
is reported beside it.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

# 'Tgfbr2_1', 'Tgfbr2_4-1', 'Ifngr2_3', 'Jak2_1'. 'KP_1-2' (unassigned KP lesion) and
# 'periphery' carry no knockout identity.
_KNOCKOUT_LABEL = re.compile(r"^(?P<gene>[A-Z][A-Za-z0-9]*?)_(?P<lesion>\d+(?:-\d+)?)$")
_UNASSIGNED = {"KP"}

_ENSEMBL_PREFIX = "ENS"


def parse_knockout(label: object) -> str | None:
    """Return the knockout gene named by a ``phenotypes`` label, or None.

    Args:
        label: A value of the authors' ``phenotypes`` column.

    Returns:
        The gene symbol, or None for unassigned lesions, periphery and missing labels.
    """
    if not isinstance(label, str):
        return None
    match = _KNOCKOUT_LABEL.match(label)
    if match is None or match.group("gene") in _UNASSIGNED:
        return None
    return match.group("gene")


def knockout_spot_table(annotations: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Count labelled spots and lesions per knockout and section.

    Args:
        annotations: Section name -> the authors' spot annotation table.

    Returns:
        Columns ``section``, ``knockout``, ``n_spots``, ``n_lesions``.
    """
    rows: list[dict[str, object]] = []
    for section, ann in annotations.items():
        labels = ann["phenotypes"]
        genes = labels.map(parse_knockout)
        frame = pd.DataFrame({"gene": genes, "label": labels}).dropna(subset=["gene"])
        for gene, grp in frame.groupby("gene"):
            rows.append(
                {
                    "section": section,
                    "knockout": gene,
                    "n_spots": len(grp),
                    "n_lesions": grp["label"].nunique(),
                }
            )
    return pd.DataFrame(rows, columns=["section", "knockout", "n_spots", "n_lesions"])


def spot_category_counts(annotations: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Count spots per section as normal, tumour with knockout label, unassigned tumour, periphery.

    Args:
        annotations: Section name -> the authors' spot annotation table.

    Returns:
        One row per section with the four counts and the total.
    """
    rows: list[dict[str, object]] = []
    for section, ann in annotations.items():
        labels = ann["phenotypes"]
        is_ko = labels.map(parse_knockout).notna()
        is_periphery = labels.eq("periphery")
        is_tumour_cluster = labels.notna()
        rows.append(
            {
                "section": section,
                "n_spots": len(ann),
                "kmeans_tumour": int(ann["kmeans"].eq("tumor").sum()),
                "knockout_labelled": int(is_ko.sum()),
                "unassigned_tumour": int((is_tumour_cluster & ~is_ko & ~is_periphery).sum()),
                "periphery": int(is_periphery.sum()),
            }
        )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class GateResult:
    """Outcome of the week-1 data gate.

    Attributes:
        library_size: Knockouts in the screen library.
        labelled: Knockouts with at least one labelled spot.
        testable_upper_bound: Knockouts with >= min_spots * min_units labelled spots.
        testable_replicated: Knockouts with >= min_spots spots in >= min_units sections.
        min_testable: Knockouts the proposal requires.
        passed: Whether the upper bound reaches ``min_testable``.
    """

    library_size: int
    labelled: tuple[str, ...]
    testable_upper_bound: tuple[str, ...]
    testable_replicated: tuple[str, ...]
    min_testable: int
    passed: bool


def evaluate_gate(
    table: pd.DataFrame,
    *,
    library_size: int,
    min_spots: int = 30,
    min_units: int = 2,
    min_testable: int = 10,
) -> GateResult:
    """Apply the pre-set week-1 rule to a knockout spot table.

    Args:
        table: Output of :func:`knockout_spot_table`.
        library_size: Number of knockouts in the screen library.
        min_spots: Spots needed per domain (or section).
        min_units: Domains (or sections) needed per knockout.
        min_testable: Testable knockouts needed to pass.

    Returns:
        The gate outcome.
    """
    totals = table.groupby("knockout")["n_spots"].sum()
    upper = sorted(str(g) for g in totals[totals >= min_spots * min_units].index)
    big_sections = table[table["n_spots"] >= min_spots].groupby("knockout")["section"].nunique()
    replicated = sorted(str(g) for g in big_sections[big_sections >= min_units].index)
    return GateResult(
        library_size=library_size,
        labelled=tuple(sorted(str(g) for g in totals.index)),
        testable_upper_bound=tuple(upper),
        testable_replicated=tuple(replicated),
        min_testable=min_testable,
        passed=len(upper) >= min_testable,
    )


def guaranteed_by_median(
    n_knockouts: int,
    median_units: float,
    *,
    min_spots: int = 30,
    min_units: int = 2,
) -> int:
    """Knockouts that must clear the gate's upper bound, given only a published median.

    When per-knockout counts are not public, a reported median still fixes a floor: at least
    ``ceil(n / 2)`` knockouts have ``median_units`` or more cells. If the median reaches the
    upper bound (``min_spots * min_units``), that many knockouts are guaranteed to pass it.
    This is a necessary condition only; domain balance still needs the cell-level data.

    Args:
        n_knockouts: Knockouts in the library.
        median_units: Published median cells (or spots) per knockout.
        min_spots: Cells needed per domain.
        min_units: Domains needed per knockout.

    Returns:
        The guaranteed number of knockouts at or above the upper bound (0 if the median is
        below it).
    """
    if median_units < min_spots * min_units:
        return 0
    return (n_knockouts + 1) // 2


def non_ensembl_features(feature_ids: Sequence[str], feature_names: Sequence[str]) -> list[str]:
    """List features that are not Ensembl genes (transgenes such as WPRE or barcode tags).

    Args:
        feature_ids: Feature identifiers from the count matrix.
        feature_names: Feature symbols, aligned with ``feature_ids``.

    Returns:
        Names of the non-Ensembl features, in file order.
    """
    return [
        name
        for fid, name in zip(feature_ids, feature_names, strict=True)
        if not fid.startswith(_ENSEMBL_PREFIX)
    ]
