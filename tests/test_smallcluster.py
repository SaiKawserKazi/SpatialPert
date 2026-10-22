import numpy as np
import pytest

from spatialpert.smallcluster import (
    cr2_two_niche_min_q,
    cr2_vcov,
    mixedlm_min_q,
    wild_cluster_bootstrap_p,
    within_lesion_signflip_p,
)


def test_cr2_with_singleton_clusters_equals_hc2():
    import statsmodels.api as sm

    rng = np.random.default_rng(0)
    n = 60
    x = np.column_stack([np.ones(n), rng.integers(0, 2, n).astype(float)])
    y = x @ np.array([1.0, 0.5]) + rng.normal(size=n) * (1 + x[:, 1])
    fit = sm.OLS(y, x).fit(cov_type="HC2")
    e = y - x @ np.linalg.lstsq(x, y, rcond=None)[0]
    v = cr2_vcov(x, e, np.arange(n))
    np.testing.assert_allclose(v, fit.cov_params(), rtol=1e-8)


def _lesions(shift, seed=0, *, n_les=8, per=40, within=True):
    genes = 5
    rng = np.random.default_rng(seed)
    les = np.repeat(np.arange(n_les), per)
    niche = np.tile(np.repeat([0, 1], per // 2), n_les) if within else (les % 2)
    r = rng.normal(size=(n_les * per, genes)) + rng.normal(0, 1, size=(n_les, genes))[les]
    r[niche == 1] += shift
    return r, niche, les


@pytest.mark.parametrize("fn", [wild_cluster_bootstrap_p, within_lesion_signflip_p])
def test_permutation_style_methods_detect_strong_effect(fn):
    r, n, l_ = _lesions(1.5)
    kw = {"n_boot": 199} if fn is wild_cluster_bootstrap_p else {}
    assert fn(r, n, l_, **kw) < 0.05


def test_q_style_methods_detect_strong_effect():
    r, n, l_ = _lesions(1.5)
    assert cr2_two_niche_min_q(r, n, l_) < 0.05
    assert mixedlm_min_q(r, n, l_) < 0.05


def test_untestable_designs_return_nan():
    r, n, l_ = _lesions(0.0)
    one = np.zeros_like(n)
    assert np.isnan(wild_cluster_bootstrap_p(r, one, l_))
    assert np.isnan(within_lesion_signflip_p(r[:80], n[:80], l_[:80]))  # 2 lesions only


def test_signflip_exact_small_case():
    # 3 lesions, contrasts +1 each in one gene: obs = 1; null over 8 flips gives
    # mean in {+-1, +-1/3}: squares {1, 1/9}; P(>=1) = 2/8
    r = np.zeros((12, 1))
    niche = np.tile([0, 0, 1, 1], 3)
    les = np.repeat([0, 1, 2], 4)
    r[niche == 1] = 1.0
    assert within_lesion_signflip_p(r, niche, les, min_cells=2) == pytest.approx(0.25)


def test_null_calibration_signflip_and_wcb():
    ps_sf = [within_lesion_signflip_p(*_lesions(0.0, seed=s), n_flip=999) for s in range(60)]
    ps_wb = [
        wild_cluster_bootstrap_p(*_lesions(0.0, seed=s, within=False), n_boot=99, seed=s)
        for s in range(40)
    ]
    assert np.mean(np.array(ps_sf) < 0.05) <= 0.15
    assert np.mean(np.array(ps_wb) < 0.05) <= 0.15


def test_vectorised_summed_wald_equals_loop():
    from spatialpert.smallcluster import _design, _summed_wald_loop, summed_wald

    rng = np.random.default_rng(4)
    niche = rng.integers(0, 3, 300)
    cl = rng.integers(0, 12, 300)
    r = rng.normal(size=(300, 20))
    x, _ = _design(niche)
    _, cidx = np.unique(cl, return_inverse=True)
    assert summed_wald(r, x, cidx, 12) == pytest.approx(_summed_wald_loop(r, x, cidx, 12), rel=1e-8)


def test_niche_contrast_recovers_effect_and_covers():
    from spatialpert.smallcluster import niche_contrast

    rng = np.random.default_rng(5)
    cover, est = 0, []
    for _ in range(200):
        les = np.repeat(np.arange(12), 15)
        niche = np.where(np.arange(les.size) % 15 < 8, 0, 1)  # every lesion straddles both niches
        base = rng.normal(0, 0.5, 12)[les]
        r = (base + 0.3 * (niche == 1) + rng.normal(size=les.size))[:, None]
        fit = niche_contrast(r, niche, les, 0, 1)
        assert fit is not None
        est.append(fit["estimate"][0])
        cover += fit["lo"][0] <= 0.3 <= fit["hi"][0]
    assert np.mean(est) == pytest.approx(0.3, abs=0.05)
    assert 0.9 <= cover / 200 <= 0.99
