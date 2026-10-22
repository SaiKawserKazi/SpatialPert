import numpy as np
import pytest
from scipy import sparse

from spatialpert.effects import (
    benjamini_hochberg,
    between_domain_ss,
    clone_patches,
    fit_control_model,
    leave_one_domain_out_r,
    matched_pseudo_patches,
    patch_domains,
    patch_permutation_test,
    residuals,
)


def test_control_model_recovers_domain_means_and_slope():
    rng = np.random.default_rng(0)
    n = 2000
    dom = rng.integers(0, 3, n)
    cov = rng.normal(size=(n, 1))
    mu = np.array([[1.0, 0.0], [2.0, 1.0], [3.0, -1.0]])
    beta = np.array([[0.5, -0.2]])
    y = mu[dom] + cov @ beta + rng.normal(0, 0.01, size=(n, 2))
    m = fit_control_model(sparse.csr_matrix(y), dom, cov)
    np.testing.assert_allclose(m.beta, beta, atol=0.01)
    r = residuals(m, y, dom, cov)
    assert np.abs(r).max() < 0.05


def test_between_domain_ss_matches_definition():
    rng = np.random.default_rng(1)
    r = rng.normal(size=(30, 4))
    d = np.repeat([0, 1, 2], 10)
    sums = np.stack([r[d == k].sum(0) for k in range(3)])
    sizes = np.array([10.0, 10.0, 10.0])
    grand = r.mean(0)
    expected = sum(10 * ((r[d == k].mean(0) - grand) ** 2).sum() for k in range(3))
    assert between_domain_ss(sums, sizes) == pytest.approx(expected)


def test_clone_patches_and_majority():
    coords = np.array([[0, 0], [10, 0], [100, 0], [110, 0], [500, 0]], dtype=float)
    p = clone_patches(coords, link=30)
    assert p[0] == p[1]
    assert p[2] == p[3]
    assert len(set(p)) == 3
    dom = np.array([1, 1, 2, 0, 0])
    pd_ = patch_domains(p, dom)
    assert pd_[p[0]] == 1
    assert pd_[p[2]] == 0  # tie 2 vs 0 -> smallest id
    assert pd_[p[4]] == 0


def _ko(n_patches=20, cells=10, genes=5, shift=0.0, seed=0):
    rng = np.random.default_rng(seed)
    patch = np.repeat(np.arange(n_patches), cells)
    pdom = np.arange(n_patches) % 2
    r = rng.normal(size=(n_patches * cells, genes))
    r[pdom[patch] == 1] += shift
    return r, patch, pdom


def test_permutation_detects_domain_shift():
    r, patch, pdom = _ko(shift=1.0)
    res = patch_permutation_test(r, patch, pdom, n_perm=200)
    assert res.p_value < 0.01
    assert res.n_cells == 200


def test_permutation_null_is_calibrated():
    ps = [patch_permutation_test(*_ko(seed=s), n_perm=99, seed=s).p_value for s in range(60)]
    assert np.mean(np.array(ps) < 0.05) <= 0.15


def test_bh_matches_scipy():
    from scipy.stats import false_discovery_control

    rng = np.random.default_rng(3)
    for p in (np.array([0.01, 0.04, 0.03, 0.5]), rng.uniform(size=50)):
        np.testing.assert_allclose(benjamini_hochberg(p), false_discovery_control(p))


def test_leave_one_domain_out_r():
    base = np.array([1.0, -1.0, 2.0, 0.0])
    eff = np.stack([base, base * 2, base + 0.01])
    r = leave_one_domain_out_r(eff, np.array([10.0, 10.0, 10.0]))
    assert np.all(r > 0.99)


def test_matched_pseudo_patches():
    rng = np.random.default_rng(0)
    coords = rng.uniform(0, 1000, size=(500, 2))
    dom = (coords[:, 0] > 500).astype(int)
    idx, pid = matched_pseudo_patches(coords, dom, np.array([5, 3]), np.array([0, 1]), rng)
    assert len(idx) == 8
    assert len(set(idx)) == 8
    assert np.all(dom[idx[pid == 0]] == 0)  # seed domain; neighbours mostly same side
    assert np.bincount(pid).tolist() == [5, 3]


def test_drop_small_niches_applies_minimum_after_patching():
    from spatialpert.effects import drop_small_niches

    # patches: 0 (niche A, 40 cells), 1 (niche B, 5 cells), 2 (niche C, 35 cells)
    patch = np.repeat([0, 1, 2], [40, 5, 35])
    pdom = np.array([0, 1, 2])
    mask, p2, d2 = drop_small_niches(patch, pdom, 30)
    assert mask.sum() == 75
    assert not mask[40:45].any()
    assert d2.tolist() == [0, 2]
    assert sorted(set(p2.tolist())) == [0, 1]


def test_drop_small_niches_cascades_to_empty():
    from spatialpert.effects import drop_small_niches

    mask, _, d2 = drop_small_niches(np.repeat([0, 1], [10, 12]), np.array([0, 1]), 30)
    assert mask.sum() == 0
    assert d2.size == 0


def test_local_residuals_remove_regional_offset():
    from spatialpert.effects import local_residuals

    rng = np.random.default_rng(0)
    cc = rng.uniform(0, 1000, size=(2000, 2))
    offset = (cc[:, 0] > 500) * 5.0  # region-specific baseline
    zc = offset[:, None] + rng.normal(0, 0.1, size=(2000, 2))
    ko_xy = np.array([[100.0, 500.0], [900.0, 500.0]])
    zk = np.array([[0.0, 0.0], [5.0, 5.0]]) + 1.0  # same true effect (+1) in both regions
    r = local_residuals(zk, ko_xy, zc, cc, k=30)
    np.testing.assert_allclose(r, 1.0, atol=0.1)


def test_local_residuals_exclude_mask():
    from spatialpert.effects import local_residuals

    cc = np.array([[0.0, 0.0], [1.0, 0.0], [100.0, 0.0]])
    zc = np.array([[10.0], [10.0], [0.0]])
    r = local_residuals(
        np.array([[0.0]]),
        np.array([[0.0, 0.0]]),
        zc,
        cc,
        exclude=np.array([True, True, False]),
        k=1,
    )
    assert r[0, 0] == 0.0  # only the far, non-excluded control is used
