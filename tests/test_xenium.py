import h5py
import numpy as np
import pytest
from scipy import sparse

from spatialpert.xenium import (
    guide_cell_table,
    guide_probe_indices,
    negative_probe_rate,
    probe_target,
    read_xenium_h5,
)


@pytest.mark.parametrize(
    ("name", "target"),
    [("PRSS23-2", "PRSS23"), ("C10orf53-3", "C10orf53"), ("NKX2-1", "NKX2"), ("ACTB", "ACTB")],
)
def test_probe_target(name, target):
    assert probe_target(name) == target


@pytest.fixture
def xen(tmp_path):
    # 4 cells; features: gene, guide probe, negative probe x2
    dense = np.array(
        [
            [5, 1, 0, 0],
            [3, 2, 1, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
        ],
        dtype=np.int32,
    ).T  # features x cells
    csc = sparse.csc_matrix(dense)
    p = tmp_path / "x.h5"
    with h5py.File(p, "w") as fh:
        g = fh.create_group("matrix")
        g["data"], g["indices"], g["indptr"] = csc.data, csc.indices, csc.indptr
        g["shape"] = np.array(dense.shape)
        g["barcodes"] = np.array([b"c1", b"c2", b"c3", b"c4"])
        f = g.create_group("features")
        f["id"] = np.array([b"ENSG1", b"PRSS23-2", b"NegControlProbe_1", b"NegControlProbe_2"])
        f["name"] = np.array([b"ACTB", b"PRSS23-2", b"NegControlProbe_1", b"NegControlProbe_2"])
        f["feature_type"] = np.array(
            [
                b"Gene Expression",
                b"Gene Expression",
                b"Negative Control Probe",
                b"Negative Control Probe",
            ]
        )
    return read_xenium_h5(p)


def test_guide_probe_indices(xen):
    assert guide_probe_indices(xen) == [1]


def test_negative_probe_rate(xen):
    # probe 1 positive in 1/4 cells, probe 2 in 0/4 -> mean 0.125
    assert negative_probe_rate(xen) == pytest.approx(0.125)


def test_guide_cell_table(xen):
    t = guide_cell_table(xen)
    row = t.iloc[0]
    assert row["target"] == "PRSS23"
    assert row["cells"] == 2
    assert row["expected_by_chance"] == pytest.approx(0.5)
    assert row["enrichment"] == pytest.approx(4.0)
    assert guide_cell_table(xen, min_umis=2).iloc[0]["cells"] == 1
