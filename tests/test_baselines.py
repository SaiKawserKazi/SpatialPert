import numpy as np
import pytest
from scipy import stats

from spatialpert.baselines import (
    cluster_robust_any,
    patch_means,
    pseudobulk_anova_any,
    studentized_permutation_p,
    welch_f,
)


def test_patch_means():
    r = np.array([[1.0], [3.0], [10.0]])
    m, s, d = patch_means(r, np.array([0, 0, 1]), np.array([5, 6]))
    np.testing.assert_allclose(m.ravel(), [2.0, 10.0])
    assert s.tolist() == [2.0, 1.0]
    assert d.tolist() == [5, 6]


def test_welch_f_matches_reference_formula():
    # two groups: Welch F equals the squared Welch t
    rng = np.random.default_rng(0)
    a, b = rng.normal(0, 1, 12), rng.normal(0.5, 3, 7)
    f = welch_f(np.concatenate([a, b])[:, None], np.array([0] * 12 + [1] * 7))[0]
    t = stats.ttest_ind(a, b, equal_var=False).statistic
    assert f == pytest.approx(t**2, rel=1e-9)


def _units(shift, seed=0, n_per=8, genes=4, sd=(1.0, 1.0)):
    rng = np.random.default_rng(seed)
    g = np.repeat([0, 1], n_per)
    m = np.vstack([rng.normal(0, sd[0], (n_per, genes)), rng.normal(shift, sd[1], (n_per, genes))])
    return m, g


def test_studentized_detects_shift_and_needs_two_per_group():
    m, g = _units(2.0)
    assert studentized_permutation_p(m, g, n_perm=199, seed=0) < 0.01
    assert np.isnan(studentized_permutation_p(m[:9], np.array([0] * 8 + [1]), n_perm=9, seed=0))


def test_studentized_calibrated_under_unequal_variance():
    # Valid with enough units per group (15: 0.05 over 200 seeds). With 6 per group it
    # inflates (0.165 vs 0.07 for the plain statistic): small-sample Welch variances are
    # unstable. That regime is judged by the benchmark, not asserted here.
    ps = [
        studentized_permutation_p(*_units(0.0, seed=s, n_per=15, sd=(0.3, 3.0)), n_perm=99, seed=s)
        for s in range(80)
    ]
    assert np.mean(np.array(ps) < 0.05) <= 0.12


def test_pseudobulk_and_cluster_robust():
    m, g = _units(3.0)
    assert pseudobulk_anova_any(m, g) is True
    assert pseudobulk_anova_any(m[:9], np.array([0] * 8 + [1])) is None
    rng = np.random.default_rng(2)
    cl = np.repeat(np.arange(10), 20)
    grp = (cl >= 5).astype(int)
    r = rng.normal(size=(200, 3)) + rng.normal(size=10)[cl][:, None]
    r[grp == 1] += 3.0
    assert cluster_robust_any(r, grp, cl) is True
    assert cluster_robust_any(r[:40], grp[:40], cl[:40]) is None  # one niche only
