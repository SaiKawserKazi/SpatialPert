"""Small-cluster inference for niche dependence with the lesion as cluster (AMENDMENT-008, A1).

All functions take one knockout's cell residuals (cells x genes), each cell's niche and each
cell's lesion id, and return a p-value (NaN when the design cannot be tested).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import stats

from spatialpert.effects import benjamini_hochberg

WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])


def _design(niche: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    labs = np.unique(niche)
    x = np.column_stack([np.ones(niche.size)] + [(niche == g).astype(np.float64) for g in labs[1:]])
    return x, labs


def _testable(niche: np.ndarray, lesion: np.ndarray, min_cells: int) -> bool:
    """Every niche must appear (>= min_cells cells) in at least two lesions."""
    labs = np.unique(niche)
    if labs.size < 2:
        return False
    for g in labs:
        les = np.unique(lesion[niche == g])
        if sum(np.sum((niche == g) & (lesion == c)) >= min_cells for c in les) < 2:
            return False
    return True


def summed_wald(resid: np.ndarray, x: np.ndarray, cluster_idx: np.ndarray, n_c: int) -> float:
    """Sum over genes of the CR1 Wald statistic for all niche coefficients (vectorised)."""
    n, p = x.shape
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ resid  # p x G
    e = resid - x @ beta
    scores = np.zeros((n_c, p, resid.shape[1]))
    np.add.at(scores, cluster_idx, x[:, :, None] * e[:, None, :])
    corr = (n_c / (n_c - 1)) * ((n - 1) / max(n - p, 1))
    meat = np.einsum("cpg,cqg->gpq", scores, scores)
    v = corr * np.einsum("ip,gpq,qj->gij", xtx_inv, meat, xtx_inv)[:, 1:, 1:]
    b = beta[1:].T  # G x (p-1)
    vinv = np.linalg.pinv(v)
    return float(np.einsum("gi,gij,gj->", b, vinv, b))


def _summed_wald_loop(resid: np.ndarray, x: np.ndarray, cluster_idx: np.ndarray, n_c: int) -> float:
    """Reference (per-gene loop) implementation, kept for testing the vectorised version."""
    n, p = x.shape
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ resid
    e = resid - x @ beta
    scores = np.zeros((n_c, p, resid.shape[1]))
    np.add.at(scores, cluster_idx, x[:, :, None] * e[:, None, :])
    corr = (n_c / (n_c - 1)) * ((n - 1) / max(n - p, 1))
    total = 0.0
    for g in range(resid.shape[1]):
        s = scores[:, :, g]
        v = corr * xtx_inv @ (s.T @ s) @ xtx_inv
        b = beta[1:, g]
        total += float(b @ np.linalg.pinv(v[1:, 1:]) @ b)
    return total


def wild_cluster_bootstrap_p(  # noqa: PLR0913 - keyword-only options
    resid: np.ndarray,
    niche: np.ndarray,
    lesion: np.ndarray,
    *,
    n_boot: int = 499,
    seed: int = 0,
    min_cells: int = 3,
) -> float:
    """Restricted wild cluster bootstrap (Webb 6-point weights) of the summed Wald statistic.

    Args:
        resid: Cells x genes residuals.
        niche: Niche per cell.
        lesion: Lesion (cluster) id per cell.
        n_boot: Bootstrap draws.
        seed: Random seed.
        min_cells: Cells a lesion needs in a niche to count as visiting it.

    Returns:
        Bootstrap p-value, or NaN if some niche is seen in fewer than two lesions.
    """
    if not _testable(niche, lesion, min_cells):
        return float("nan")
    x, _ = _design(niche)
    clus, cidx = np.unique(lesion, return_inverse=True)
    obs = summed_wald(resid, x, cidx, clus.size)
    mu = resid.mean(0)
    u = resid - mu  # residuals of the restricted (no-niche) model
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_boot):
        w = rng.choice(WEBB, size=clus.size)[cidx]
        if summed_wald(mu + u * w[:, None], x, cidx, clus.size) >= obs:
            exceed += 1
    return (1 + exceed) / (1 + n_boot)


def cr2_vcov(x: np.ndarray, e: np.ndarray, cluster: np.ndarray) -> np.ndarray:
    """CR2 (Bell-McCaffrey) cluster-robust covariance of OLS coefficients for one outcome.

    Args:
        x: Cells x p design.
        e: OLS residuals for one gene.
        cluster: Cluster id per cell.

    Returns:
        p x p covariance matrix.
    """
    xtx_inv = np.linalg.inv(x.T @ x)
    h = x @ xtx_inv @ x.T
    meat = np.zeros((x.shape[1], x.shape[1]))
    for c in np.unique(cluster):
        ix = np.where(cluster == c)[0]
        w, vec = np.linalg.eigh(np.eye(ix.size) - h[np.ix_(ix, ix)])
        ac = vec @ np.diag(np.clip(w, 1e-10, None) ** -0.5) @ vec.T
        s = x[ix].T @ ac @ e[ix]
        meat += np.outer(s, s)
    out: np.ndarray = xtx_inv @ meat @ xtx_inv
    return out


def cr2_two_niche_min_q(
    resid: np.ndarray, niche: np.ndarray, lesion: np.ndarray, *, min_cells: int = 3
) -> float:
    """CR2 (Bell-McCaffrey) t-test of the two largest niches per gene, Imbens-Kolesar df.

    Returns the minimum BH q over genes (so 'any gene q < 0.05' is 'result < 0.05'), or NaN.
    """
    labs, counts = np.unique(niche, return_counts=True)
    if labs.size < 2:
        return float("nan")
    a, b = labs[np.argsort(-counts)[:2]]
    fit = niche_contrast(resid, niche, lesion, a, b, min_cells=min_cells)
    if fit is None:
        return float("nan")
    pv = 2 * stats.t.sf(np.abs(fit["estimate"] / fit["se"]), fit["df"])
    return float(benjamini_hochberg(pv).min())


def niche_contrast(  # noqa: PLR0913 - keyword-only options
    resid: np.ndarray,
    niche: np.ndarray,
    lesion: np.ndarray,
    a: int,
    b: int,
    *,
    min_cells: int = 3,
    level: float = 0.95,
) -> dict[str, np.ndarray] | None:
    """Knockout-by-niche contrast Psi_g = Delta_b,g - Delta_a,g with a CR2 interval per gene.

    Residuals are knockout minus niche-matched controls, so the mean residual in a niche
    estimates the knockout effect there (Delta) and the difference between two niches estimates
    the niche dependence Psi. Standard errors are CR2 (Bell-McCaffrey) clustered on lesions with
    Imbens-Kolesar degrees of freedom.

    Args:
        resid: Cells x genes residuals of one knockout.
        niche: Niche per cell.
        lesion: Lesion id per cell.
        a: Reference niche.
        b: Comparison niche.
        min_cells: Cells a lesion needs in a niche to count as visiting it.
        level: Confidence level.

    Returns:
        Arrays (one per gene) ``estimate``, ``se``, ``lo``, ``hi`` and the scalar ``df``, or
        None if either niche is seen in fewer than two lesions.
    """
    keep = np.isin(niche, [a, b])
    r, nic, les = resid[keep], niche[keep], lesion[keep]
    if np.unique(nic).size < 2 or not _testable(nic, les, min_cells):
        return None
    n = r.shape[0]
    x = np.column_stack([np.ones(n), (nic == b).astype(np.float64)])
    xtx_inv = np.linalg.inv(x.T @ x)
    h = x @ xtx_inv @ x.T
    clus = np.unique(les)
    ell = np.array([0.0, 1.0])
    adj, gcols = {}, []
    eye = np.eye(n)
    for c in clus:
        ix = np.where(les == c)[0]
        m = np.eye(ix.size) - h[np.ix_(ix, ix)]
        w, vec = np.linalg.eigh(m)
        w = np.clip(w, 1e-10, None)
        adj[c] = (ix, vec @ np.diag(w**-0.5) @ vec.T)
        # Imbens-Kolesar G column: (I - H)[:, ix] A_c X_c (X'X)^-1 ell
        gcols.append((eye - h)[:, ix] @ adj[c][1] @ x[ix] @ xtx_inv @ ell)
    gm = np.column_stack(gcols)
    lam = np.linalg.eigvalsh(gm.T @ gm)
    df = float(lam.sum() ** 2 / np.sum(lam**2))
    beta = xtx_inv @ x.T @ r
    e = r - x @ beta
    se = np.ones(r.shape[1])
    for g in range(r.shape[1]):
        meat = np.zeros((2, 2))
        for ix, ac in adj.values():
            s = x[ix].T @ ac @ e[ix, g]
            meat += np.outer(s, s)
        v = xtx_inv @ meat @ xtx_inv
        se[g] = np.sqrt(max(v[1, 1], 1e-300))
    half = stats.t.ppf(0.5 + level / 2, df) * se
    return {
        "estimate": beta[1],
        "se": se,
        "lo": beta[1] - half,
        "hi": beta[1] + half,
        "df": np.array(df),
    }


def within_lesion_signflip_p(  # noqa: PLR0913 - keyword-only options
    resid: np.ndarray,
    niche: np.ndarray,
    lesion: np.ndarray,
    *,
    min_cells: int = 3,
    n_flip: int = 2000,
    seed: int = 0,
) -> float:
    """Paired within-lesion test: niche contrast inside lesions, sign-flip null.

    The two niches with most cells are compared. Each lesion with >= ``min_cells`` cells in both
    gives a contrast vector (mean residual in niche B minus niche A); the statistic is the
    squared mean contrast summed over genes. Each lesion's baseline cancels in its own contrast.

    Returns:
        Sign-flip p-value (exact when 2^L <= n_flip), or NaN with fewer than 3 paired lesions.
    """
    labs, counts = np.unique(niche, return_counts=True)
    if labs.size < 2:
        return float("nan")
    a, b = labs[np.argsort(-counts)[:2]]
    cons = []
    for c in np.unique(lesion):
        ma, mb = (lesion == c) & (niche == a), (lesion == c) & (niche == b)
        if ma.sum() >= min_cells and mb.sum() >= min_cells:
            cons.append(resid[mb].mean(0) - resid[ma].mean(0))
    if len(cons) < 3:
        return float("nan")
    cm = np.stack(cons)
    obs = float(np.sum(cm.mean(0) ** 2))
    n_l = cm.shape[0]
    if 2**n_l <= n_flip:
        signs = np.array(
            [[1 - 2 * ((k >> i) & 1) for i in range(n_l)] for k in range(2**n_l)], dtype=np.float64
        )
    else:
        rng = np.random.default_rng(seed)
        signs = np.asarray(rng.choice([-1.0, 1.0], size=(n_flip, n_l)), dtype=np.float64)
    null = np.sum((signs @ cm / n_l) ** 2, axis=1)
    return float(np.mean(null >= obs - 1e-12))


def mixedlm_min_q(
    resid: np.ndarray, niche: np.ndarray, lesion: np.ndarray, *, min_cells: int = 3
) -> float:
    """Statsmodels MixedLM on lesion x niche pseudobulk, lesion random intercept; min BH q.

    No Satterthwaite / Kenward-Roger df are available in statsmodels; a Wald chi-square is used.
    """
    import statsmodels.formula.api as smf  # noqa: PLC0415 - optional dependency

    if not _testable(niche, lesion, min_cells):
        return float("nan")
    rows = []
    for c in np.unique(lesion):
        for g in np.unique(niche):
            m = (lesion == c) & (niche == g)
            if m.sum() >= min_cells:
                rows.append((c, g, resid[m].mean(0)))
    les = np.array([r[0] for r in rows])
    nic = np.array([r[1] for r in rows])
    y = np.stack([r[2] for r in rows])
    import pandas as pd  # noqa: PLC0415

    pv = np.ones(y.shape[1])
    for j in range(y.shape[1]):
        df = pd.DataFrame({"y": y[:, j], "niche": nic.astype(str), "lesion": les.astype(str)})
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = smf.mixedlm("y ~ C(niche)", df, groups=df["lesion"]).fit(reml=True)
            names = [k for k in fit.params.index if k.startswith("C(niche)")]
            b = fit.params[names].to_numpy()
            v = fit.cov_params().loc[names, names].to_numpy()
            w = float(b @ np.linalg.pinv(v) @ b)
            pv[j] = float(stats.chi2.sf(w, len(names)))
        except (ValueError, np.linalg.LinAlgError):
            pv[j] = 1.0
    return float(benjamini_hochberg(pv).min())
