import gzip

import h5py
import numpy as np
from scipy import sparse

from spatialpert.visium import read_10x_h5, read_positions, read_spot_annotation


def test_read_10x_h5_transposes_to_spots_by_features(tmp_path):
    # features x barcodes, as Space Ranger stores it
    dense = np.array([[1, 0, 2], [0, 3, 0]], dtype=np.int32)  # 2 features, 3 spots
    csc = sparse.csc_matrix(dense)
    path = tmp_path / "m.h5"
    with h5py.File(path, "w") as fh:
        g = fh.create_group("matrix")
        g["data"], g["indices"], g["indptr"] = csc.data, csc.indices, csc.indptr
        g["shape"] = np.array(dense.shape)
        g["barcodes"] = np.array([b"AAA-1", b"CCC-1", b"GGG-1"])
        f = g.create_group("features")
        f["id"] = np.array([b"ENSMUSG1", b"WPRE"])
        f["name"] = np.array([b"Actb", b"WPRE"])
    m = read_10x_h5(path)
    assert m.counts.shape == (3, 2)
    np.testing.assert_array_equal(m.counts.toarray(), dense.T)
    assert m.barcodes == ["AAA-1", "CCC-1", "GGG-1"]
    assert m.feature_names == ["Actb", "WPRE"]


def test_read_positions_and_annotation(tmp_path):
    pos = tmp_path / "p.csv.gz"
    with gzip.open(pos, "wt") as fh:
        fh.write("AAA-1,1,0,0,100,200\nCCC-1,0,0,2,100,220\n")
    df = read_positions(pos)
    assert list(df["in_tissue"]) == [1, 0]
    assert df.loc[1, "pxl_col"] == 220

    ann = tmp_path / "a.csv.gz"
    with gzip.open(ann, "wt") as fh:
        fh.write("barcode,kmeans,phenotypes\nAAA-1,tumor,Tgfbr2_1\nCCC-1,normal,NA\n")
    a = read_spot_annotation(ann)
    assert a["phenotypes"].isna().tolist() == [False, True]
