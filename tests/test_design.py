import numpy as np
import pytest

from spatialpert.design import DesignModel, effective_n, icc_oneway, plan_screen


def test_icc_recovers_simulated_value():
    rng = np.random.default_rng(0)
    g = np.repeat(np.arange(300), 20)
    # between-group var 0.25, within 0.75 -> ICC 0.25
    v = rng.normal(0, 0.5, 300)[g] + rng.normal(0, np.sqrt(0.75), g.size)
    assert icc_oneway(v, g) == pytest.approx(0.25, abs=0.04)


def test_effective_n_limits():
    assert effective_n(10, 40, 0.0) == pytest.approx(400)
    assert effective_n(10, 40, 1.0) == pytest.approx(10)


def test_plan_screen_monotone_in_lesions_and_floor():
    m = DesignModel(rho=0.2, sigma=1.0)
    p = plan_screen(np.array([2, 4, 10, 20]), 40, 0.5, m)
    assert p[0] == 0.0
    assert np.all(np.diff(p[1:]) > 0)
    # same total cells, more lesions -> more power
    assert plan_screen(20, 20, 0.5, m) > plan_screen(5, 80, 0.5, m)


def test_plan_screen_empirical_interpolates_table():
    from spatialpert.design import plan_screen_empirical
    from spatialpert.design_table import POWER_TABLE

    t = POWER_TABLE
    i = t["lesions"].index(10)
    assert plan_screen_empirical(10) == pytest.approx(t["log_scale_lambda0.5"][i])
    assert plan_screen_empirical(2) == 0.0
    p = plan_screen_empirical(np.array([5, 10, 20]))
    assert np.all(np.diff(p) >= 0)
    with pytest.raises(ValueError, match="scenario"):
        plan_screen_empirical(10, "nope")
