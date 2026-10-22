import numpy as np
import pytest

from spatialpert.calibration import calibrate_p, lesion_exchangeability, matched_null_p


def test_calibrate_p_is_uniform_under_exchangeability():
    rng = np.random.default_rng(0)
    pool = rng.choice([0.004, 0.02, 0.3, 0.7, 1.0], size=200)  # discrete, with ties
    obs = rng.choice([0.004, 0.02, 0.3, 0.7, 1.0], size=2000)
    pc = np.array([calibrate_p(p, pool, seed=i) for i, p in enumerate(obs)])
    assert np.mean(pc <= 0.05) == pytest.approx(0.05, abs=0.015)
    assert np.isnan(calibrate_p(0.1, np.array([np.nan])))


def test_lesion_exchangeability_detects_baselines_and_unequal_spread():
    rng = np.random.default_rng(1)
    les = np.repeat(np.arange(40), 20)
    niche = np.repeat(np.arange(40) % 2, 20)
    noise = rng.normal(size=(les.size, 5))
    flat = lesion_exchangeability(noise, les, niche)
    assert flat["variance_ratio"] == pytest.approx(1.0, abs=0.4)
    base = rng.normal(size=(40, 5)) * np.where(np.arange(40) % 2, 3.0, 1.0)[:, None]
    strong = lesion_exchangeability(noise + base[les], les, niche)
    assert strong["variance_ratio"] > 10
    assert strong["niche_sd_ratio"] > 2


def test_matched_null_p_runs_the_same_test_on_pseudo_knockouts():
    rng = np.random.default_rng(2)
    n, g = 4000, 6
    xy = rng.uniform(0, 4000, size=(n, 2))
    dom = (xy[:, 0] > 2000).astype(int)
    y = rng.normal(size=(n, g)) + dom[:, None]
    lab = np.array(["control"] * n, dtype=object)
    centres = np.array([[500, 500], [500, 3000], [3000, 800], [3200, 3200]])
    for c in centres:
        lab[np.argsort(np.hypot(*(xy - c).T))[:40]] = "KO"
    p = matched_null_p(y, xy, lab, dom, "KO", n_null=5, n_perm=99, seed=3)
    assert p.shape == (5,)
    assert np.sum(~np.isnan(p)) >= 3  # pseudo-knockouts must actually be tested
    assert np.all((p[~np.isnan(p)] > 0) & (p[~np.isnan(p)] <= 1))
