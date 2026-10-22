"""AMENDMENT-010: niche dependence of immune infiltration in Perturb-map imaging (S-BIAD267).

Implements the registered plan exactly: lesions (60 px single linkage, >= 30 labelled cells,
>= 80% majority Gene), per-lobe marker thresholds (95th percentile), outcome = logit fraction
of marker-positive unlabelled cells within 40 px of the lesion, niche = k-means (K = 2) on the
100-300 µm ring composition (cells within 200 px of any labelled cell excluded), residual
against same-mouse same-niche reference lesions, within-mouse permutation, matched-null
calibration and CR2 (mouse clusters) intervals.

Usage:
  step36_perturbmap_niche.py --design-only   lesion/niche counts per knockout (no outcomes)
  step36_perturbmap_niche.py                 full analysis -> reports/step36_perturbmap.json
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib.path import Path as Poly
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

from spatialpert.calibration import calibrate_p
from spatialpert.effects import benjamini_hochberg
from spatialpert.smallcluster import niche_contrast

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data/interim/perturbmap"
ANN = ROOT / "data/raw/S-BIAD267"
MARKERS = ("CD8", "CD4", "B220")
IFN = ["Jak2", "Stat1", "Stat2", "Ifnar", "Ifngr2", "Irf1", "Irf9", "Ptpn2", "Socs1"]
HALO, RING_EXCL, RING_IN, RING_OUT = 40.0, 200.0, 200.0, 600.0


def polygons(mouse: str, lobe: str) -> list[tuple[str, np.ndarray]]:
    """Curated lesion polygons (Pro-Code name, vertices) of one lobe (AMENDMENT-010a)."""
    f = ANN / f"{mouse}_Lymphoid_{lobe}_curated.fix.json"
    if not f.exists():
        return []
    d = json.loads(f.read_text())
    out = []
    for x in d if isinstance(d, list) else d.get("features", []):
        name = (x.get("properties") or {}).get("name")
        geom = x.get("geometry") or {}
        if not name or geom.get("type") not in ("Polygon", "MultiPolygon"):
            continue
        ring = geom["coordinates"][0] if geom["type"] == "Polygon" else max(geom["coordinates"], key=lambda q: len(q[0]))[0]
        out.append((re.sub(r"_\d+$", "", name), np.asarray(ring, dtype=np.float64)))
    return out


def lesions_of_lobe(df: pd.DataFrame, code2gene: dict[str, str]) -> list[dict[str, object]]:
    lab = df[df.Gene.notna()]
    unl = df[df.Gene.isna()]
    thr = {m: float(np.quantile(df[m], 0.95)) for m in MARKERS}
    pos = {m: (unl[m].to_numpy() > thr[m]) for m in MARKERS}
    all_pos = {m: (df[m].to_numpy() > thr[m]) for m in MARKERS}
    uxy, axy, lxy = unl[["X", "Y"]].to_numpy(), df[["X", "Y"]].to_numpy(), lab[["X", "Y"]].to_numpy()
    if len(lxy) == 0:
        return []
    utree, atree = cKDTree(uxy), cKDTree(axy)
    near_lab = cKDTree(lxy).query(axy, k=1, distance_upper_bound=RING_EXCL)[0] < np.inf
    out = []
    for code, poly in polygons(str(df.mouse.iloc[0]), str(df.lobe.iloc[0])):
        gene = code2gene.get(code)
        if gene is None:
            continue
        inside = Poly(poly).contains_points(lxy)
        if inside.sum() < 30:
            continue
        pts = lxy[inside]
        in_unl = np.where(Poly(poly).contains_points(uxy))[0]
        rim = np.concatenate([np.asarray(h, dtype=int) for h in utree.query_ball_point(pts, r=HALO)])
        halo = np.unique(np.concatenate([in_unl, rim]))
        cen, rad = pts.mean(0), np.max(np.hypot(*(pts - pts.mean(0)).T))
        ring = np.asarray(atree.query_ball_point(cen, r=rad + RING_OUT), dtype=int)
        dist = np.hypot(*(axy[ring] - cen).T)
        ring = ring[(dist >= rad + RING_IN) & ~near_lab[ring]]
        rec: dict[str, object] = {"gene": gene, "procode": code, "cells": int(inside.sum()), "x": float(cen[0]), "y": float(cen[1]),
                                  "halo_n": int(halo.size), "ring_n": int(ring.size)}
        for mk in MARKERS:
            k_pos = int(pos[mk][halo].sum()) if halo.size else 0
            rec[f"y_{mk}"] = float(np.log((k_pos + 0.5) / (halo.size - k_pos + 0.5)))
            rec[f"ring_{mk}"] = float(all_pos[mk][ring].mean()) if ring.size else np.nan
        rec["ring_density"] = float(ring.size / max(np.pi * ((rad + RING_OUT) ** 2 - (rad + RING_IN) ** 2), 1.0))
        out.append(rec)
    return out


def build() -> pd.DataFrame:
    files = sorted(D.glob("*_Lymphoid.parquet"))
    pairs = pd.concat([pd.read_parquet(f, columns=["procode", "Gene"]).dropna() for f in files])
    code2gene = pairs.groupby("procode").Gene.agg(lambda v: v.value_counts().index[0]).to_dict()
    rows = []
    for f in files:
        df = pd.read_parquet(f)
        for r in lesions_of_lobe(df, code2gene):
            r |= {"mouse": str(df.mouse.iloc[0]), "lobe": f"{df.mouse.iloc[0]}_{df.lobe.iloc[0]}"}
            rows.append(r)
    les = pd.DataFrame(rows)
    feats = les[[f"ring_{m}" for m in MARKERS] + ["ring_density"]]
    ok = feats.notna().all(1) & (les.ring_n >= 20)
    les = les[ok].reset_index(drop=True)
    z = (feats[ok] - feats[ok].mean()) / feats[ok].std()
    km = KMeans(2, n_init=20, random_state=0).fit(z)
    immune = int(np.argmax(km.cluster_centers_[:, :3].sum(1)))
    les["niche"] = (km.labels_ == immune).astype(int)  # 1 = immune-rich
    return les


def residual(les: pd.DataFrame, gene: str, y: str) -> np.ndarray:
    ref = les[les.gene != gene].groupby(["mouse", "niche"])[y].mean()
    k = les[les.gene == gene]
    return np.array([r[y] - ref.get((r.mouse, r.niche), np.nan) for _, r in k.iterrows()])


def psi(e: np.ndarray, niche: np.ndarray) -> float:
    return float(e[niche == 1].mean() - e[niche == 0].mean())


def within_mouse_perm(e: np.ndarray, niche: np.ndarray, mouse: np.ndarray, rng: np.random.Generator, n: int = 9999) -> float:
    obs = abs(psi(e, niche))
    groups = [np.where(mouse == m)[0] for m in np.unique(mouse)]
    hit = 0
    for _ in range(n):
        nn = niche.copy()
        for g in groups:
            nn[g] = rng.permutation(niche[g])
        if (nn == 1).any() and (nn == 0).any() and abs(psi(e, nn)) >= obs:
            hit += 1
    return (1 + hit) / (1 + n)


def analyse(les: pd.DataFrame, gene: str, y: str, rng: np.random.Generator) -> dict[str, float]:
    k = les[les.gene == gene].reset_index(drop=True)
    e = residual(les, gene, y)
    ok = ~np.isnan(e)
    k, e = k[ok].reset_index(drop=True), e[ok]
    nic, mou = k.niche.to_numpy(), k.mouse.to_numpy()
    p = within_mouse_perm(e, nic, mou, rng)
    # matched null: pseudo-knockouts from reference lesions with the same per-mouse counts
    ref = les[les.gene != gene]
    null = []
    for _ in range(1000):
        pick = pd.concat([ref[ref.mouse == m].sample(n=min(c, int((ref.mouse == m).sum())), random_state=int(rng.integers(1 << 31)))
                          for m, c in k.mouse.value_counts().items()])
        sub = les.drop(index=pick.index).copy()
        tag = "__pseudo__"
        sub = pd.concat([sub, pick.assign(gene=tag)])
        ep = residual(sub, tag, y)
        pm = sub[sub.gene == tag]
        okp = ~np.isnan(ep)
        nn = pm.niche.to_numpy()[okp]
        null.append(within_mouse_perm(ep[okp], nn, pm.mouse.to_numpy()[okp], rng, n=199) if len(np.unique(nn)) == 2 else np.nan)
    fit = niche_contrast(e[:, None], nic, mou_codes(mou), 0, 1, min_cells=1)
    lo, hi = (float(fit["lo"][0]), float(fit["hi"][0])) if fit is not None else (np.nan, np.nan)
    return {"lesions": int(len(k)), "immune_rich": int(nic.sum()), "mice": int(len(np.unique(mou))),
            "delta_main": float(e.mean()), "psi": psi(e, nic), "ci_lo": lo, "ci_hi": hi,
            "p_perm": p, "p_calibrated": calibrate_p(p, np.array(null), seed=1)}


def mou_codes(mouse: np.ndarray) -> np.ndarray:
    return np.unique(mouse, return_inverse=True)[1]


def main() -> None:
    les = build()
    les.to_csv(ROOT / "reports/step36_lesions_design.csv", index=False, columns=[c for c in les.columns if not c.startswith("y_")])
    cnt = les.groupby("gene").agg(lesions=("niche", "size"), immune_rich=("niche", "sum"), mice=("mouse", "nunique"))
    cnt["immune_poor"] = cnt.lesions - cnt.immune_rich
    elig = cnt[(cnt.lesions >= 20) & (cnt.immune_rich >= 5) & (cnt.immune_poor >= 5) & (cnt.mice >= 3)]
    design = {"lobes": int(les.lobe.nunique()), "mice": int(les.mouse.nunique()), "lesions": int(len(les)),
              "immune_rich_fraction": float(les.niche.mean()), "eligible": elig.index.tolist(),
              "per_gene": cnt.sort_values("lesions", ascending=False).to_dict(orient="index")}
    if "--design-only" in sys.argv:
        print(json.dumps({k: v for k, v in design.items() if k != "per_gene"}, indent=2))
        print(cnt.sort_values("lesions", ascending=False).head(40).to_string())
        return
    rng = np.random.default_rng(36)
    res: dict[str, object] = {"design": design, "results": {}}
    for y in ("y_CD8", "y_CD4", "y_B220"):
        tab = {g: analyse(les, g, y, rng) for g in elig.index}
        q = benjamini_hochberg(np.array([v["p_calibrated"] for v in tab.values()]))
        for (g, v), qq in zip(tab.items(), q, strict=True):
            v["q_calibrated"] = float(qq)
            v["ifn_pathway"] = g in IFN
            v["niche_dependent"] = bool(qq < 0.10 and (v["ci_lo"] > 0 or v["ci_hi"] < 0))
        res["results"][y] = tab  # type: ignore[index]
        print(y, {g: (round(v["psi"], 2), round(v["p_calibrated"], 3), round(v["q_calibrated"], 3)) for g, v in tab.items()}, flush=True)
    (ROOT / "reports/step36_perturbmap.json").write_text(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    _ = stats  # scipy kept for downstream extensions
    main()
