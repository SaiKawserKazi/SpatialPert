import numpy as np
import pytest
from scipy import sparse

from spatialpert.domains import (
    call_tumour,
    choose_domains,
    knockout_calls,
    marker_score,
    neighbour_counts,
    neighbourhood_mean,
    normalise_log1p,
)


def test_normalise_log1p_scales_to_median():
    counts = sparse.csr_matrix(np.array([[2.0, 2.0], [8.0, 0.0], [0.0, 0.0]]))
    out = normalise_log1p(counts).toarray()
    # median of nonzero totals (4, 8) = 6
    np.testing.assert_allclose(out[0], np.log1p([3.0, 3.0]))
    np.testing.assert_allclose(out[1], np.log1p([6.0, 0.0]))
    assert out[2].sum() == 0


def test_marker_score_and_tumour_call():
    rng = np.random.default_rng(0)
    hi = rng.normal(3, 0.3, size=(200, 2))
    lo = rng.normal(0, 0.3, size=(200, 2))
    x = sparse.csr_matrix(np.clip(np.vstack([hi, lo]), 0, None))
    s = marker_score(x, [0, 1])
    flag = call_tumour(s)
    assert flag[:200].mean() > 0.95
    assert flag[200:].mean() < 0.05


def test_knockout_calls():
    g = sparse.csr_matrix(np.array([[0, 0], [1, 0], [0, 2], [1, 1]]))
    labels = knockout_calls(g, ["A", "B"])
    assert labels.tolist() == ["control", "A", "B", "multi"]


def test_neighbourhood_mean_excludes_self():
    coords = np.array([[0.0, 0.0], [10.0, 0.0], [1000.0, 0.0]])
    ref = np.array([True, True, True])
    vals = np.array([[1.0], [3.0], [5.0]])
    feats, total = neighbourhood_mean(coords, vals, ref, radius=50, length=50)
    assert feats[0, 0] == pytest.approx(3.0)  # only cell 1 is a neighbour of cell 0
    assert feats[1, 0] == pytest.approx(1.0)
    assert np.isnan(feats[2, 0])  # isolated
    assert total[2] == 0


def test_neighbourhood_mean_reference_only():
    coords = np.array([[0.0, 0.0], [5.0, 0.0], [0.0, 5.0]])
    ref = np.array([False, True, False])
    vals = np.array([[7.0]])
    feats, _ = neighbourhood_mean(coords, vals, ref, radius=50, length=50)
    assert feats[0, 0] == pytest.approx(7.0)
    assert feats[2, 0] == pytest.approx(7.0)
    assert np.isnan(feats[1, 0])  # its only reference neighbour is itself


def test_neighbour_counts():
    coords = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0], [500.0, 0.0]])
    flag = np.array([True, False, True, True])
    dens, frac = neighbour_counts(coords, flag, radius=15)
    assert dens.tolist() == [1, 2, 1, 0]
    assert frac[0] == pytest.approx(1.0)
    assert frac[1] == pytest.approx(0.0)
    assert np.isnan(frac[3])


def test_choose_domains_finds_three_blobs():
    rng = np.random.default_rng(0)
    centres = np.array([[0, 0], [10, 0], [0, 10]])
    x = np.vstack([c + rng.normal(0, 0.5, size=(100, 2)) for c in centres])
    km, scores = choose_domains(x, ks=(2, 3, 4), subsample=300)
    assert km.n_clusters == 3
    assert set(scores) == {2, 3, 4}


def test_tile_ids():
    from spatialpert.pipeline import tile_ids

    c = np.array([[10.0, 10.0], [240.0, 20.0], [260.0, 10.0], [10.0, 300.0]])
    t = tile_ids(c, 250.0)
    assert t[0] == t[1]
    assert len(set(t.tolist())) == 3


def test_neighbourhood_mean_query_matches_full_version():
    from spatialpert.domains import neighbourhood_mean_query

    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 300, size=(200, 2))
    vals = rng.normal(size=(200, 3))
    ref = np.ones(200, dtype=bool)
    ref[:20] = False  # first 20 cells are queries, not references
    full, _ = neighbourhood_mean(xy, vals[ref], ref, radius=100, length=50)
    q = neighbourhood_mean_query(xy[:20], xy[ref], vals[ref], radius=100, length=50)
    np.testing.assert_allclose(q, full[:20], rtol=1e-10)
