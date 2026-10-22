import numpy as np
import pytest
from scipy.spatial import cKDTree

from spatialpert.gmn import gmn_pvalue, niche_ss, place_geometry


def test_niche_ss_zero_when_niches_equal():
    r = np.array([[1.0], [1.0], [1.0], [1.0]])
    assert niche_ss(r, np.array([0, 0, 1, 1])) == pytest.approx(0.0)
    r2 = np.array([[0.0], [0.0], [2.0], [2.0]])
    assert niche_ss(r2, np.array([0, 0, 1, 1])) == pytest.approx(4.0)


def test_place_geometry_sizes_and_separation():
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 5000, size=(20000, 2))
    tree = cKDTree(xy)
    idx = place_geometry(tree, xy, [30, 10, 5], rng, min_sep=800)
    assert idx.size == 45
    assert len(set(idx.tolist())) == 45
    c = [xy[idx[:30]].mean(0), xy[idx[30:40]].mean(0), xy[idx[40:]].mean(0)]
    assert min(np.hypot(*(c[i] - c[j])) for i in range(3) for j in range(i + 1, 3)) > 600


def test_gmn_pvalue():
    assert gmn_pvalue(5.0, np.array([1.0, 2.0, 6.0])) == pytest.approx(0.5)
    assert gmn_pvalue(10.0, np.zeros(99)) == pytest.approx(0.01)


def test_place_geometry_returns_true_lesion_ids():
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 5000, size=(20000, 2))
    idx, les = place_geometry(cKDTree(xy), xy, [12, 5, 8], rng, min_sep=600, return_lesions=True)
    assert np.bincount(les).tolist() == [12, 5, 8]
    assert idx.size == les.size
