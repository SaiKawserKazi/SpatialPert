import gzip

import pandas as pd
import pytest

from spatialpert.dbit import (
    guide_target,
    normalise_pixel,
    pixel_knockouts,
    read_position_pixels,
    read_sg_counts,
)


@pytest.mark.parametrize(
    ("guide", "gene"),
    [
        ("sg57430_ATG4A", "ATG4A"),
        ("sg16727_S100A4", "S100A4"),
        ("Tgfbr2_sg2", "Tgfbr2"),
        ("Rps19_sg1", "Rps19"),
        ("NTC_3", "NTC"),
        ("Non-Targeting_12", "NTC"),
        ("garbage", None),
    ],
)
def test_guide_target(guide, gene):
    assert guide_target(guide) == gene


def test_normalise_pixel():
    assert normalise_pixel("bc10x13") == "10x13"
    assert normalise_pixel("10x13") == "10x13"


def test_read_sg_counts_and_pixel_knockouts(tmp_path):
    p = tmp_path / "sg.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(
            "bc1x1\tTgfbr2_sg1\t1\t\nbc1x1\tTgfbr2_sg2\t2\t\nbc1x2\tNTC_1\t1\t\nbc2x2\tjunk\t5\t\n"
        )
    sg = read_sg_counts(p)
    assert sg["pixel"].tolist() == ["1x1", "1x1", "1x2", "2x2"]
    ko = pixel_knockouts(sg)
    assert ko.set_index(["pixel", "gene"])["umis"].to_dict() == {
        ("1x1", "Tgfbr2"): 3,
        ("1x2", "NTC"): 1,
    }
    assert len(pixel_knockouts(sg, min_umis=2)) == 1


def test_read_position_pixels(tmp_path):
    p = tmp_path / "pos.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(",1x1,2x1,3x1\n")
    assert read_position_pixels(p) == ["1x1", "2x1", "3x1"]


def test_pixel_knockouts_empty():
    sg = pd.DataFrame({"pixel": ["1x1"], "guide": ["junk"], "count": [1], "gene": [None]})
    assert pixel_knockouts(sg).empty


def test_read_dbit_expression_pixels_as_rows_astro(tmp_path):
    from spatialpert.dbit import read_dbit_expression

    p = tmp_path / "a.tsv"
    p.write_text(
        "\tACTB__exon\tACTB__transcript\tKRT20__exon\tXYZ__exon\n1x1\t3\t9\t1\t5\n2x1\t0\t9\t2\t5\n"
    )
    df = read_dbit_expression(p, ["ACTB", "KRT20"])
    assert sorted(df.columns) == ["ACTB", "KRT20"]
    assert df.loc["1x1", "ACTB"] == 3.0
    assert df.loc["2x1", "KRT20"] == 2.0


def test_read_dbit_expression_genes_as_rows(tmp_path):
    from spatialpert.dbit import read_dbit_expression

    pix = [f"{i}x1" for i in range(1, 31)]
    lines = ["\t".join(pix)]
    lines.append("ACTB\t" + "\t".join(["1"] * 30))
    lines.append("OTHER\t" + "\t".join(["7"] * 30))
    p = tmp_path / "b.tsv"
    p.write_text("\n".join(lines) + "\n")
    df = read_dbit_expression(p, ["ACTB"])
    assert df.shape == (30, 1)
    assert df.loc["5x1", "ACTB"] == 1.0
