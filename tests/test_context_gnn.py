import numpy as np
import pytest

torch = pytest.importorskip("torch")

from spatialpert.context_gnn import ContextEncoder, ring_means  # noqa: E402


def test_ring_means_exclude_self_and_split_rings():
    xy = np.array([[0.0, 0.0], [10.0, 0.0], [100.0, 0.0]])
    x = np.array([[1.0], [2.0], [4.0]])
    out = ring_means(xy, xy, x, (50.0, 150.0), same=True)
    assert out[0, 0] == pytest.approx(2.0)  # inner ring of cell 0: cell 1 only (self excluded)
    assert out[0, 1] == pytest.approx(4.0)  # outer ring: cell 2


def test_encoder_learns_context_and_leaves_out_cells():
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 2000, size=(3000, 2))
    region = (xy[:, 0] > 1000).astype(float)
    x = np.column_stack([region + rng.normal(0, 0.3, 3000), rng.normal(size=3000)])
    enc = ContextEncoder.fit(xy, x, epochs=150, dim=4)
    z = enc.embed_reference()
    # the embedding separates the two regions although it never sees a cell's own features
    from sklearn.linear_model import LogisticRegression

    assert LogisticRegression().fit(z, region).score(z, region) > 0.9
    excl = np.zeros(3000, dtype=bool)
    excl[:10] = True
    assert enc.embed(xy[:10], exclude=excl).shape == (10, 4)
