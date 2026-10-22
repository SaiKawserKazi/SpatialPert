"""Learned graph context encoder for control-defined niches (review 3, Priority 2).

Each cell's context is the mean expression (PCs) of control neighbours in two rings
(0-r1 and r1-r2 µm, the cell itself always excluded), mapped by a small neural network to an
embedding. The network is trained on control cells to predict a cell's own PCs from its
context only. One-hop aggregation is deliberate: with two message-passing hops a cell's own
features return to it through its neighbours (i -> j -> i), so the encoder could learn the
cell rather than its surroundings. Knockout (or pseudo-knockout) cells are removed from the
reference before their context is computed, so they never contribute to the context that
labels them.

Requires the optional dependency ``torch`` (``pip install spatialpert[gnn]``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree


def ring_means(
    query_xy: np.ndarray,
    ref_xy: np.ndarray,
    ref_x: np.ndarray,
    rings: tuple[float, float],
    *,
    same: bool,
) -> np.ndarray:
    """Mean reference features in the rings [0, r1) and [r1, r2) around each query point.

    Args:
        query_xy: Query positions x 2.
        ref_xy: Reference (control) positions x 2.
        ref_x: Reference features.
        rings: (r1, r2) in µm.
        same: Query and reference are the same cells (self-pairs are dropped).

    Returns:
        Query x (2 * features); rows are zero where a ring is empty.
    """
    r1, r2 = rings
    d = cKDTree(query_xy).sparse_distance_matrix(cKDTree(ref_xy), r2, output_type="coo_matrix")
    keep = d.row != d.col if same else np.ones(d.row.size, dtype=bool)
    rows, cols, dist = d.row[keep], d.col[keep], d.data[keep]
    parts = []
    for inner in (dist < r1, dist >= r1):
        w = sparse.csr_matrix(
            (np.ones(int(inner.sum())), (rows[inner], cols[inner])),
            shape=(query_xy.shape[0], ref_xy.shape[0]),
        )
        deg = np.asarray(w.sum(1)).ravel()
        inv = np.divide(1.0, deg, out=np.zeros_like(deg), where=deg > 0)
        parts.append(np.asarray(sparse.diags(inv) @ w @ ref_x))
    return np.hstack(parts)


@dataclass
class ContextEncoder:
    """Trained encoder plus the control reference it was trained on.

    Attributes:
        model: The torch module (context features -> embedding -> predicted own PCs).
        ref_xy: Control-cell coordinates.
        ref_x: Control-cell input features (expression PCs).
        rings: Inner and outer ring radii (µm).
    """

    model: Any
    ref_xy: np.ndarray
    ref_x: np.ndarray
    rings: tuple[float, float]

    @classmethod
    def fit(  # noqa: PLR0913 - training options
        cls,
        ref_xy: np.ndarray,
        ref_x: np.ndarray,
        *,
        rings: tuple[float, float] = (50.0, 150.0),
        dim: int = 16,
        epochs: int = 400,
        seed: int = 0,
    ) -> ContextEncoder:
        """Train on control cells to predict each cell's features from its context only.

        Args:
            ref_xy: Control cells x 2 (µm).
            ref_x: Control cells x features (e.g. 30 expression PCs).
            rings: Ring radii (µm).
            dim: Embedding size.
            epochs: Full-batch Adam epochs.
            seed: Torch seed.

        Returns:
            The trained encoder.
        """
        import torch  # noqa: PLC0415

        torch.manual_seed(seed)
        ctx = torch.tensor(ring_means(ref_xy, ref_xy, ref_x, rings, same=True), dtype=torch.float32)
        y = torch.tensor(ref_x, dtype=torch.float32)
        model = _mlp(ctx.shape[1], dim, ref_x.shape[1])
        opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-5)
        for _ in range(epochs):
            opt.zero_grad()
            loss = ((model.head(model.encode(ctx)) - y) ** 2).mean()
            loss.backward()
            opt.step()
        model.eval()
        return cls(model=model, ref_xy=ref_xy, ref_x=ref_x, rings=rings)

    def embed(self, query_xy: np.ndarray, exclude: np.ndarray | None = None) -> np.ndarray:
        """Context embedding of query positions, optionally leaving some control cells out.

        Args:
            query_xy: Query positions x 2.
            exclude: Boolean mask over control cells that must not serve as context (e.g. the
                pseudo-knockout's own cells).

        Returns:
            Query x dim embedding.
        """
        import torch  # noqa: PLC0415

        keep = np.ones(self.ref_xy.shape[0], dtype=bool) if exclude is None else ~exclude
        ctx = ring_means(query_xy, self.ref_xy[keep], self.ref_x[keep], self.rings, same=False)
        with torch.no_grad():
            z = self.model.encode(torch.tensor(ctx, dtype=torch.float32))
        out: np.ndarray = z.numpy()
        return out

    def embed_reference(self) -> np.ndarray:
        """Embedding of every control cell (self excluded from its own context)."""
        import torch  # noqa: PLC0415

        ctx = ring_means(self.ref_xy, self.ref_xy, self.ref_x, self.rings, same=True)
        with torch.no_grad():
            z = self.model.encode(torch.tensor(ctx, dtype=torch.float32))
        out: np.ndarray = z.numpy()
        return out


def _mlp(n_ctx: int, dim: int, n_out: int) -> Any:
    import torch  # noqa: PLC0415
    from torch import nn  # noqa: PLC0415

    class Encoder(nn.Module):  # type: ignore[misc,unused-ignore]
        def __init__(self) -> None:
            super().__init__()
            self.enc = nn.Sequential(nn.Linear(n_ctx, 64), nn.ReLU(), nn.Linear(64, dim))
            self.dec = nn.Linear(dim, n_out)

        def encode(self, ctx: torch.Tensor) -> torch.Tensor:
            out: torch.Tensor = self.enc(ctx)
            return out

        def head(self, z: torch.Tensor) -> torch.Tensor:
            out: torch.Tensor = self.dec(z)
            return out

    return Encoder()
