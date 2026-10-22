import numpy as np
import pytest

from spatialpert.api import location_dependence as run_test


def _tissue(shift: float, seed: int = 0):
    rng = np.random.default_rng(seed)
    n_ctrl, genes = 2000, 6
    coords = rng.uniform(0, 2000, size=(n_ctrl, 2))
    dom = (coords[:, 0] > 1000).astype(int)
    labels = np.array(["control"] * n_ctrl, dtype=object)
    y = rng.normal(size=(n_ctrl, genes)) + dom[:, None] * 0.5
    # knockout clones: 12 small patches per domain
    kc, kd = [], []
    for d in (0, 1):
        for c in rng.uniform([100 + 1000 * d, 100], [900 + 1000 * d, 1900], size=(12, 2)):
            kc.append(c + rng.normal(0, 5, size=(8, 2)))
            kd.append(np.full(8, d))
    kc, kd = np.vstack(kc), np.concatenate(kd)
    ky = rng.normal(size=(kd.size, genes)) + kd[:, None] * 0.5 + 1.0
    ky[kd == 1, :3] += shift
    return (
        np.vstack([y, ky]),
        np.vstack([coords, kc]),
        np.concatenate([labels, np.array(["KO"] * kd.size, dtype=object)]),
        np.concatenate([dom, kd]),
    )


def test_detects_domain_dependent_effect():
    y, c, lab, dom = _tissue(shift=1.5)
    res = run_test(y, c, lab, dom, unit="patch", n_perm=200)
    row = res.table.set_index("knockout").loc["KO"]
    assert row["p"] < 0.01
    assert row["units"] == 24
    assert res.effects["KO"].shape == (2, 6)


def test_null_rejection_rate_is_controlled():
    # 400-seed check gave 25/400 = 0.0625 (binomial p = 0.25 vs 0.05); here a fast version.
    ps = np.array(
        [
            run_test(*_tissue(shift=0.0, seed=s), unit="patch", n_perm=99).table.iloc[0]["p"]
            for s in range(40)
        ]
    )
    assert np.mean(ps < 0.05) <= 0.15
    assert 0.3 < np.median(ps) < 0.7


def test_tile_unit_and_validation():
    y, c, lab, dom = _tissue(shift=1.5)
    res = run_test(y, c, lab, dom, unit="tile", tile_um=300, n_perm=200)
    assert len(res.table) == 1
    with pytest.raises(ValueError, match="unit"):
        run_test(y, c, lab, dom, unit="bad")
    with pytest.raises(ValueError, match="same number"):
        run_test(y[:-1], c, lab, dom)


def test_tutorial_runs_and_separates_knockouts():
    import runpy
    from pathlib import Path

    ns = runpy.run_path(str(Path(__file__).resolve().parents[1] / "docs" / "tutorial.py"))
    t = ns["res"].table.set_index("knockout")
    assert t.loc["KO_A", "p"] < 0.01
    assert t.loc["KO_B", "p"] > 0.05


def test_from_anndata_duck_typed():
    import types

    import pandas as pd

    from spatialpert.api import from_anndata

    y, c, lab, dom = _tissue(shift=1.5)
    ad = types.SimpleNamespace(
        X=y,
        obs=pd.DataFrame({"ko": lab, "niche": dom}),
        obsm={"spatial": c},
        var_names=[f"g{i}" for i in range(y.shape[1])],
        layers={},
    )
    res = from_anndata(ad, label_key="ko", domain_key="niche", unit="patch", n_perm=200)
    assert res.table.iloc[0]["p"] < 0.01


def test_lesion_unit_and_studentized_statistic():
    y, c, lab, dom = _tissue(shift=1.5)
    res = run_test(y, c, lab, dom, unit="lesion", lesion_um=20, statistic="studentized", n_perm=200)
    assert res.table.iloc[0]["p"] < 0.05
    with pytest.raises(ValueError, match="statistic"):
        run_test(y, c, lab, dom, statistic="bad")
