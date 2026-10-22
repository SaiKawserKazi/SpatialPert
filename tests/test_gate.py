import numpy as np
import pandas as pd
import pytest

from spatialpert.gate import (
    evaluate_gate,
    knockout_spot_table,
    non_ensembl_features,
    parse_knockout,
    spot_category_counts,
)


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        ("Tgfbr2_1", "Tgfbr2"),
        ("Tgfbr2_4-2", "Tgfbr2"),
        ("Ifngr2_3", "Ifngr2"),
        ("Jak2_1", "Jak2"),
        ("KP_1-1", None),
        ("periphery", None),
        (np.nan, None),
        (None, None),
    ],
)
def test_parse_knockout(label, expected):
    assert parse_knockout(label) == expected


def _ann(labels, kmeans=None):
    return pd.DataFrame(
        {
            "phenotypes": labels,
            "kmeans": kmeans or ["tumor" if isinstance(x, str) else "normal" for x in labels],
        }
    )


def test_knockout_spot_table_counts_spots_and_lesions():
    anns = {
        "s1": _ann(["Tgfbr2_4-1"] * 3 + ["Tgfbr2_4-2"] * 2 + ["KP_4-1"] * 4 + [np.nan]),
        "s2": _ann(["Jak2_1"] * 2 + ["periphery"]),
    }
    table = knockout_spot_table(anns).set_index(["section", "knockout"])
    assert table.loc[("s1", "Tgfbr2"), "n_spots"] == 5
    assert table.loc[("s1", "Tgfbr2"), "n_lesions"] == 2
    assert table.loc[("s2", "Jak2"), "n_spots"] == 2
    assert len(table) == 2


def test_spot_category_counts_partitions_tumour_spots():
    ann = _ann(["Tgfbr2_1", "KP_1-1", "KP_1-1", "periphery", np.nan, np.nan])
    row = spot_category_counts({"s": ann}).iloc[0]
    assert row["knockout_labelled"] == 1
    assert row["unassigned_tumour"] == 2
    assert row["periphery"] == 1
    assert row["n_spots"] == 6


def _table(rows):
    return pd.DataFrame(rows, columns=["section", "knockout", "n_spots", "n_lesions"])


def test_gate_upper_bound_and_replication_differ():
    table = _table(
        [
            ("s1", "A", 70, 1),  # 70 >= 60 total, but one section only
            ("s1", "B", 30, 1),
            ("s2", "B", 30, 1),  # replicated in two sections
            ("s1", "C", 59, 1),  # just under the bound
        ]
    )
    res = evaluate_gate(table, library_size=35, min_testable=2)
    assert res.testable_upper_bound == ("A", "B")
    assert res.testable_replicated == ("B",)
    assert res.labelled == ("A", "B", "C")
    assert res.passed


def test_gate_fails_below_min_testable():
    table = _table([("s1", "A", 100, 2)])
    res = evaluate_gate(table, library_size=35, min_testable=10)
    assert not res.passed


@pytest.mark.parametrize(
    ("n", "median", "expected"),
    [(35, 125, 18), (34, 125, 17), (35, 60, 18), (35, 59, 0)],
)
def test_guaranteed_by_median(n, median, expected):
    from spatialpert.gate import guaranteed_by_median

    assert guaranteed_by_median(n, median) == expected


def test_non_ensembl_features():
    ids = ["ENSMUSG000001", "WPRE", "ENSMUSG000002", "HA-tag"]
    names = ["Actb", "WPRE", "Gapdh", "HA"]
    assert non_ensembl_features(ids, names) == ["WPRE", "HA"]
