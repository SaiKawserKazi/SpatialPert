import numpy as np
import pytest

from spatialpert.benchmark import (
    anova_any_gene,
    cell_permutation_p,
    inject_effect,
    scale_sizes,
    unit_rms,
)


def test_unit_rms():
    v = unit_rms(np.array([1.0, -3.0, 2.0]))
    assert np.mean(v**2) == pytest.approx(1.0)


def test_inject_effect_only_target_gets_extra():
    rng = np.random.default_rng(0)
    r = np.zeros((4, 3))
    dom = np.array([0, 0, 1, 1])
    out = inject_effect(r, dom, np.ones(3), delta=1.0, lam=1.0, target_domain=1, rng=rng)
    np.testing.assert_allclose(out[0], out[1])
    assert not np.allclose(out[0], out[2])
    out0 = inject_effect(r, dom, np.ones(3), delta=1.0, lam=0.0, target_domain=1, rng=rng)
    np.testing.assert_allclose(out0[0], out0[3])


def test_scale_sizes():
    s = scale_sizes(np.array([10, 30, 1]), 82)
    assert s.sum() == pytest.approx(82, abs=2)
    assert s.min() >= 1


def test_anova_and_cell_perm_detect_strong_shift():
    rng = np.random.default_rng(1)
    r = rng.normal(size=(200, 5))
    dom = np.repeat([0, 1], 100)
    r[dom == 1] += 1.0
    assert anova_any_gene(r, dom)
    assert cell_permutation_p(r, dom, n_perm=99, seed=0) < 0.05
    null = rng.normal(size=(200, 5))
    assert not anova_any_gene(null, dom)
