"""V5 P1.2-P1.4 and P2.1: null-scenario benchmark of lesion-level tests and their calibration.

All scenarios have no niche-dependent knockout effect (Psi = 0). Per test draw it records the
cell-level permutation p, the lesion permutation p, the CR2 minimum q, the lesion p calibrated
against a random-placement pool, and the lesion p calibrated against a niche-matched pool
(pseudo-lesions seeded in the same niches as the draw's lesions; pools cached by occupancy).
Scenarios (see v5/SPATIALPERT_V5_PLAN.md section 2): A, B1-B3 (niche-enriched seeds, odds
2/5/20), B4 (boundary-enriched seeds), C1 (lesion size doubled in one niche), C2 (survival
0.5 in one niche), E (heavy-tailed sizes), G1 (adjacent lesions merged by 100 um linkage),
G2 (lesions defined by 50 um linkage), H1 (knockout main effect), H2 (lesion random effects).

Usage: step40_scenarios.py TISSUE L SCENARIO[,SCENARIO...] N_DRAWS
Writes reports/v5/step40_<tissue>_L<L>_<scenario>.csv (one row per draw; appended per draw).
"""

import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step22_brain import load_brain  # noqa: E402

from spatialpert.calibration import calibrate_p
from spatialpert.effects import clone_patches, covariate_adjusted, fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import cr2_two_niche_min_q

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "v5"
N_CELLS, N_PERM, MIN_SEP = 400, 250, 400.0
POOL_M = int(os.environ.get("SP_POOL_M", 1000))
MATCH_M = int(os.environ.get("SP_MATCH_M", 100))


class Tissue:
    """Precomputed control reference: coordinates, niches, covariate-adjusted expression."""

    def __init__(self, name: str) -> None:
        a = load_brain() if name == "brain" else load_xenium_analysis(ROOT)
        self.nm = NicheModel.fit(a, n_pcs=20 if name == "brain" else 30)
        ok = self.nm.labels >= 0
        self.ci = np.where((a.is_control & ok) if name == "tumour" else ok)[0]
        self.xy = self.nm.coords[self.ci]
        self.niche = self.nm.labels[self.ci]
        model = fit_control_model(a.y[self.ci], self.niche, a.cov[self.ci])
        self.dom = model.domains
        self.mu = model.domain_means
        sd = residuals(model, a.y[self.ci[:20000]], self.niche[:20000], a.cov[self.ci[:20000]]).std(0)
        self.z = covariate_adjusted(model, a.y[self.ci], a.cov[self.ci]) / sd
        self.mu = self.mu / sd
        self.tree = cKDTree(self.xy)
        nb = self.tree.query_ball_point(self.xy, r=50.0)
        self.boundary = np.array([np.any(self.niche[np.asarray(n)] != self.niche[i]) for i, n in enumerate(nb)])
        self.k = self.dom.size

    def leave_out(self, idx: np.ndarray) -> np.ndarray:
        g = self.ci[idx]
        return self.nm.leave_out_labels(g, g)

    def resid(self, idx: np.ndarray, d: np.ndarray) -> np.ndarray:
        out: np.ndarray = self.z[idx] - self.mu[np.searchsorted(self.dom, d)]
        return out


def grow(t: Tissue, sizes: list[int], rng: np.random.Generator, *, weights: np.ndarray | None = None,
         min_sep: float = MIN_SEP, seed_niches: list[int] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Seed lesions (weighted or niche-constrained), then take each seed's nearest unused cells."""
    n = t.xy.shape[0]
    used = np.zeros(n, dtype=bool)
    seeds: list[int] = []
    for j in range(len(sizes)):
        if seed_niches is not None:
            cand = np.where(t.niche == seed_niches[j])[0]
            p = None
        else:
            cand = np.arange(n)
            p = None if weights is None else weights / weights.sum()
        for _ in range(4000):
            s = int(rng.choice(cand, p=p)) if p is not None else int(cand[rng.integers(cand.size)])
            if not used[s] and all(np.hypot(*(t.xy[s] - t.xy[q])) > min_sep for q in seeds):
                seeds.append(s)
                break
    idx, les = [], []
    for j, s in enumerate(seeds):
        _, nn = t.tree.query(t.xy[s], k=min(sizes[j] * 3 + 10, n))
        nn = np.atleast_1d(nn)
        nn = nn[~used[nn]][: sizes[j]]
        used[nn] = True
        idx.append(nn)
        les.append(np.full(nn.size, j))
    return np.concatenate(idx), np.concatenate(les)


def ss_perm_p(r: np.ndarray, unit: np.ndarray, niche_of_unit: np.ndarray, rng: np.random.Generator) -> float:
    """Between-niche SS of residuals with niche labels permuted among units (vectorised)."""
    labs = np.unique(niche_of_unit)
    if labs.size < 2:
        return float("nan")
    u_n = np.bincount(unit).astype(np.float64)
    u_s = np.zeros((u_n.size, r.shape[1]))
    np.add.at(u_s, unit, r)
    tot = u_s.sum(0)
    big_n = u_n.sum()

    def stat(lab: np.ndarray) -> np.ndarray:  # lab: (perms, units)
        out = np.zeros(lab.shape[0])
        for d in labs:
            m = (lab == d).astype(np.float64)
            n_d = m @ u_n
            s_d = m @ u_s
            with np.errstate(invalid="ignore", divide="ignore"):
                out += np.where(n_d > 0, ((s_d - np.outer(n_d / big_n, tot)) ** 2).sum(1) / n_d, 0.0)
        return out

    obs = stat(niche_of_unit[None, :])[0]
    perms = np.stack([rng.permutation(niche_of_unit) for _ in range(N_PERM)])
    null = stat(perms)
    return float((1 + np.sum(null >= obs - 1e-12)) / (1 + N_PERM))


def lesion_majority(les: np.ndarray, d: np.ndarray) -> np.ndarray:
    return np.array([np.bincount(d[les == k]).argmax() for k in range(les.max() + 1)])


def evaluate(t: Tissue, idx: np.ndarray, les: np.ndarray, rng: np.random.Generator, *, link: float | None = None,
             main: float = 0.0, lesion_re: float = 0.0) -> dict[str, float]:
    d = t.leave_out(idx)
    keep = d >= 0
    idx, d, les = idx[keep], d[keep], np.unique(les[keep], return_inverse=True)[1]
    r = t.resid(idx, d)
    if main:
        r = r + main * rng.normal(size=r.shape[1]) / np.sqrt(1.0)
    if lesion_re:
        r = r + lesion_re * rng.normal(size=(les.max() + 1, r.shape[1]))[les]
    unit = les if link is None else np.unique(clone_patches(t.xy[idx], link=link), return_inverse=True)[1]
    um = lesion_majority(unit, d)
    return {"n_units": int(unit.max() + 1), "occupancy": ",".join(str(int((um == k).sum())) for k in t.dom),
            "p_cell": ss_perm_p(r, np.arange(r.shape[0]), d, rng),
            "p_lesion": ss_perm_p(r, unit, um, rng),
            "q_cr2": cr2_two_niche_min_q(r, d, unit)}


def sizes_equal(L: int) -> list[int]:
    return [N_CELLS // L] * L


def sizes_heavy(L: int, rng: np.random.Generator) -> list[int]:
    w = rng.lognormal(0, 1.0, L)
    return [int(max(3, round(x))) for x in w / w.sum() * N_CELLS]


def draw(t: Tissue, scen: str, L: int, rng: np.random.Generator) -> dict[str, float]:  # noqa: PLR0911
    target = int(t.dom[rng.integers(t.k)])
    if scen in ("A", "H1", "H2"):
        idx, les = grow(t, sizes_equal(L), rng)
        return evaluate(t, idx, les, rng, main=0.5 if scen == "H1" else 0.0, lesion_re=0.5 if scen == "H2" else 0.0)
    if scen in ("B1", "B2", "B3"):
        odds = {"B1": 2.0, "B2": 5.0, "B3": 20.0}[scen]
        w = np.where(t.niche == target, odds, 1.0)
        idx, les = grow(t, sizes_equal(L), rng, weights=w)
        return evaluate(t, idx, les, rng)
    if scen == "B4":
        idx, les = grow(t, sizes_equal(L), rng, weights=np.where(t.boundary, 10.0, 1.0))
        return evaluate(t, idx, les, rng)
    if scen == "C1":
        # seed first, then double the size of lesions whose seed lies in the target niche
        base = sizes_equal(L)
        idx0, les0 = grow(t, base, rng)
        seeds_niche = np.array([t.niche[idx0[les0 == k][0]] for k in range(L)])
        f = np.where(seeds_niche == target, 2.0, 1.0)
        sz = [int(max(3, round(x))) for x in f / f.sum() * N_CELLS]
        idx, les = grow(t, sz, rng, seed_niches=list(seeds_niche))
        return evaluate(t, idx, les, rng)
    if scen == "C2":
        idx0, les0 = grow(t, sizes_equal(L) + sizes_equal(L)[: L // 2], rng)
        keepl = [k for k in range(les0.max() + 1) if not (t.niche[idx0[les0 == k][0]] == target and rng.random() < 0.5)][:L]
        m = np.isin(les0, keepl)
        return evaluate(t, idx0[m], np.unique(les0[m], return_inverse=True)[1], rng)
    if scen == "E":
        idx, les = grow(t, sizes_heavy(L, rng), rng)
        return evaluate(t, idx, les, rng)
    if scen == "G1":
        idx, les = grow(t, sizes_equal(L), rng, min_sep=150.0)
        return evaluate(t, idx, les, rng, link=100.0)
    if scen == "G2":
        idx, les = grow(t, sizes_equal(L), rng)
        return evaluate(t, idx, les, rng, link=50.0)
    raise ValueError(scen)


class Pools:
    """Random-placement pools (built once) and niche-matched pools (cached by occupancy)."""

    def __init__(self, t: Tissue, L: int, rng: np.random.Generator) -> None:
        self.t, self.L, self.rng = t, L, rng
        self.random: dict[str, np.ndarray] = {}
        self.matched: dict[tuple[str, str], np.ndarray] = {}

    def random_pool(self, kind: str) -> np.ndarray:
        if kind not in self.random:
            ps = []
            for _ in range(POOL_M):
                sz = sizes_heavy(self.L, self.rng) if kind == "heavy" else sizes_equal(self.L)
                link = {"link100": 100.0, "link50": 50.0}.get(kind)
                ms = MIN_SEP
                idx, les = grow(self.t, sz, self.rng, min_sep=ms)
                ps.append(evaluate(self.t, idx, les, self.rng, link=link)["p_lesion"])
            self.random[kind] = np.array(ps)
        return self.random[kind]

    def matched_pool(self, kind: str, occupancy: str, sizes: list[int]) -> np.ndarray:
        key = (kind, occupancy)
        if key not in self.matched:
            counts = [int(c) for c in occupancy.split(",")]
            niches = [int(d) for d, c in zip(self.t.dom, counts, strict=True) for _ in range(c)]
            ps = []
            for _ in range(MATCH_M):
                sz = sizes[: len(niches)] if len(sizes) >= len(niches) else sizes_equal(len(niches))
                idx, les = grow(self.t, sz, self.rng, seed_niches=niches)
                ps.append(evaluate(self.t, idx, les, self.rng)["p_lesion"])
            self.matched[key] = np.array(ps)
        return self.matched[key]


def main() -> None:
    tissue, L, scens, n_draws = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(","), int(sys.argv[4])
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    t = Tissue(tissue)
    pools = Pools(t, L, np.random.default_rng(9000 + L + (0 if tissue == "tumour" else 500)))
    print(f"{tissue} loaded {time.time() - t0:.0f}s", flush=True)
    for scen in scens:
        kind = {"E": "heavy", "G1": "link100", "G2": "link50"}.get(scen, "equal")
        f = OUT / f"step40_{tissue}_L{L}_{scen}.csv"
        done = len(pd.read_csv(f)) if f.exists() else 0
        srng = np.random.default_rng(zlib.crc32(f"{tissue}-{L}-{scen}".encode()) + done)
        for i in range(done, n_draws):
            rec = draw(t, scen, L, srng)
            rec["p_cal_random"] = calibrate_p(rec["p_lesion"], pools.random_pool(kind), seed=i)
            if kind == "equal":
                rec["p_cal_matched"] = calibrate_p(rec["p_lesion"], pools.matched_pool(kind, rec["occupancy"], sizes_equal(L)), seed=i)
            else:
                rec["p_cal_matched"] = float("nan")
            rec |= {"tissue": tissue, "L": L, "scenario": scen, "draw": i}
            pd.DataFrame([rec]).to_csv(f, mode="a", header=not f.exists(), index=False)
            if i % 50 == 0:
                print(f"{scen} {i}/{n_draws} {time.time() - t0:.0f}s pools(random={len(pools.random)}, matched={len(pools.matched)})", flush=True)
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
