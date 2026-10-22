"""Step 6: confirmatory H1-H5 (protocol §5, pooled branch per AMENDMENT-001).

Writes reports/step6_hypotheses.json and reports/step6_effects.csv.
"""

import io
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from spatialpert.dbit import pixel_knockouts, read_dbit_expression, read_sg_counts
from spatialpert.domains import normalise_log1p
from spatialpert.effects import (
    clone_patches,
    fit_control_model,
    leave_one_domain_out_r,
    matched_pseudo_patches,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.pipeline import load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
TESTABLE = ["CIITA", "DIO3", "ESAM", "PRSS23"]
N_PERM = 1000
N_PSEUDO = 200
HALLMARK_URL = (
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
)
HALLMARK_SETS = [
    "Interferon Gamma Response",
    "Interferon Alpha Response",
    "TNF-alpha Signaling via NF-kB",
    "Hypoxia",
    "Epithelial Mesenchymal Transition",
]
DBIT = ROOT / "data" / "raw" / "GSE319277"
DBIT_PAIRED = {
    "Bru.LC.6": ("GSM9514752_sgCounts_Bru.LC.6_grnalib.txt.gz", "GSM9514752_Bru.LC.6_allRNA.tsv.gz"),
    "Bru.LC.8": ("GSM9514756_sgCounts_Bru.LC.8_grnalib.txt.gz", "GSM9514756_Bru.LC.8.tsv.gz"),
    "Bru.LC.1.33.3": ("GSM9514758_sgCounts_Bru.LC.1.33.3.txt.gz", "GSM9514758_Bru.LC.1.33.3.tsv.gz"),
}


def hallmark_genes() -> set[str]:
    path = ROOT / "data" / "external" / "MSigDB_Hallmark_2020.txt"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(HALLMARK_URL, timeout=120) as resp:  # noqa: S310
            path.write_bytes(resp.read())
    genes: set[str] = set()
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts and parts[0] in HALLMARK_SETS:
            genes.update(g for g in parts[2:] if g)
    return genes


def main() -> None:  # noqa: C901, PLR0912 - one linear analysis script
    a = load_xenium_analysis(ROOT)
    ctrl = a.is_control
    dom = a.cells["domain"].to_numpy()
    coords = a.cells[["x", "y"]].to_numpy()
    model = fit_control_model(a.y[ctrl], dom[ctrl], a.cov[ctrl])
    ctrl_idx = np.where(ctrl)[0]

    # ---------------- H1 + H2 per knockout ----------------
    per_ko, eff_rows, templates = {}, [], {}
    pooled_obs, pooled_null = 0.0, np.zeros(N_PERM)
    lodo: list[float] = []
    for i, ko in enumerate(TESTABLE):
        mk = a.knockout_mask(ko)
        counts = pd.Series(dom[mk]).value_counts()
        ok_d = np.sort(counts.index[counts >= 30].to_numpy())
        sel = np.where(mk & np.isin(dom, ok_d))[0]
        patch = clone_patches(coords[sel])
        pdom = patch_domains(patch, dom[sel])
        templates[ko] = (np.bincount(patch), pdom)
        r = residuals(model, a.y[sel], dom[sel], a.cov[sel])
        res = patch_permutation_test(r, patch, pdom, n_perm=N_PERM, seed=i)
        pooled_obs += res.statistic
        pooled_null += res.null
        lodo_k = leave_one_domain_out_r(res.effects, res.sizes)
        lodo.extend(x for x in lodo_k if np.isfinite(x))
        per_ko[ko] = {
            "domains": [int(d) for d in res.domains],
            "cells_per_domain": [int(s) for s in res.sizes],
            "patches": res.n_patches,
            "T": res.statistic,
            "p_exploratory": res.p_value,
            "lodo_r": [float(x) for x in lodo_k],
        }
        for d, s, e in zip(res.domains, res.sizes, res.effects, strict=True):
            eff_rows.append({"knockout": ko, "domain": int(d), "cells": int(s), **dict(zip(a.genes, e, strict=True))})
    p_h1 = float((1 + np.sum(pooled_null >= pooled_obs)) / (1 + N_PERM))

    # H2 null: pseudo-KO sets (one per template) -> pooled median LODO r
    rng = np.random.default_rng(0)
    null_h2 = []
    for _ in range(N_PSEUDO):
        rs: list[float] = []
        for ko in TESTABLE:
            sizes, pdom = templates[ko]
            idx_l, pid = matched_pseudo_patches(coords[ctrl_idx], dom[ctrl_idx], sizes, pdom, rng)
            idx = ctrl_idx[idx_l]
            pd_ = patch_domains(pid, dom[idx])
            if np.unique(pd_).size < 2:
                continue
            r = residuals(model, a.y[idx], dom[idx], a.cov[idx])
            res = patch_permutation_test(r, pid, pd_, n_perm=0)
            rs.extend(x for x in leave_one_domain_out_r(res.effects, res.sizes) if np.isfinite(x))
        null_h2.append(float(np.median(rs)) if rs else np.nan)
    null_h2_a = np.array(null_h2)
    med_lodo = float(np.median(lodo))
    p_h2 = float((1 + np.sum(null_h2_a >= med_lodo)) / (1 + np.isfinite(null_h2_a).sum()))

    # ---------------- H3 ----------------
    eff = pd.DataFrame(eff_rows)
    energy = np.zeros(len(a.genes))
    for ko in TESTABLE:
        sub = eff[eff["knockout"] == ko]
        e = sub[a.genes].to_numpy()
        n = sub["cells"].to_numpy(dtype=float)
        ebar = (e * n[:, None]).sum(0) / n.sum()
        energy += (n[:, None] * (e - ebar) ** 2).sum(0)
    hall = hallmark_genes()
    in_set = np.array([g in hall for g in a.genes])
    if in_set.sum() >= 3 and (~in_set).sum() >= 3:
        p_h3 = float(mannwhitneyu(energy[in_set], energy[~in_set], alternative="greater").pvalue)
    else:
        p_h3 = float("nan")

    # ---------------- H4 ----------------
    tum_xy = coords
    xb = np.minimum((4 * (tum_xy[:, 0] - tum_xy[:, 0].min()) / np.ptp(tum_xy[:, 0])).astype(int), 3)
    yb = np.minimum((4 * (tum_xy[:, 1] - tum_xy[:, 1].min()) / np.ptp(tum_xy[:, 1])).astype(int), 3)
    block = xb * 4 + yb
    h4_units = []
    for ko in TESTABLE:
        mk = np.where(a.knockout_mask(ko))[0]
        r_all = residuals(model, a.y[mk], dom[mk], a.cov[mk])
        for b in np.unique(block[mk]):
            inb = block[mk] == b
            if inb.sum() < 15:
                continue
            obs = r_all[inb].mean(0)
            bdom = pd.Series(dom[mk][inb]).mode().iloc[0]
            pooled = r_all[~inb].mean(0)
            same = (~inb) & (dom[mk] == bdom)
            if same.sum() < 15:
                continue
            dom_pred = r_all[same].mean(0)
            ra = float(np.corrcoef(obs, pooled)[0, 1])
            rb = float(np.corrcoef(obs, dom_pred)[0, 1])
            h4_units.append({"knockout": ko, "block": int(b), "n": int(inb.sum()), "r_pooled": ra, "r_domain": rb})
    h4 = pd.DataFrame(h4_units)
    if len(h4) >= 3:
        diff = h4["r_domain"] - h4["r_pooled"]
        p_h4 = float(wilcoxon(diff, alternative="greater").pvalue)
        med_h4 = float(diff.median())
    else:
        p_h4, med_h4 = float("nan"), float("nan")

    # ---------------- H5 ----------------
    xen_counts = {ko: int(a.knockout_mask(ko).sum()) for ko in a.cells["label"].unique()}
    pooled_xen = {}
    for ko, n in xen_counts.items():
        if ko in ("control", "multi") or n < 30:
            continue
        mk = np.where(a.knockout_mask(ko))[0]
        pooled_xen[ko] = pd.Series(residuals(model, a.y[mk], dom[mk], a.cov[mk]).mean(0), index=a.genes)
    h5_units = []
    for sec, (sgf, exf) in DBIT_PAIRED.items():
        ko_px = pixel_knockouts(read_sg_counts(DBIT / sgf))
        mat = read_dbit_expression(DBIT / exf, a.genes)
        genes = list(mat.columns)
        if len(genes) < 10:
            h5_units.append({"section": sec, "note": f"only {len(genes)} shared genes"})
            continue
        # normalise on the analysis genes only (pixel totals of the full transcriptome are
        # not loaded); this is a within-gene-set scaling, applied identically to all pixels
        from scipy import sparse

        y = normalise_log1p(sparse.csr_matrix(mat.to_numpy()))
        yg = y.toarray()
        guided = set(ko_px["pixel"])
        is_ctrl = ~mat.index.isin(guided)
        base = yg[is_ctrl].mean(0)
        for ko, xv in pooled_xen.items():
            px = set(ko_px.loc[ko_px["gene"] == ko, "pixel"])
            rows = mat.index.isin(px)
            if rows.sum() < 30:
                continue
            dv = yg[rows].mean(0) - base
            xs = xv[genes].to_numpy()
            r = float(np.corrcoef(xs, dv)[0, 1])
            null = np.array([np.corrcoef(xs, rng.permutation(dv))[0, 1] for _ in range(N_PERM)])
            h5_units.append({"section": sec, "knockout": ko, "pixels": int(rows.sum()), "genes": len(genes), "r": r, "p_perm": float((1 + np.sum(null >= r)) / (1 + N_PERM))})

    h5_rs = [u["r"] for u in h5_units if "r" in u]
    report = {
        "branch": "pooled (AMENDMENT-001)",
        "H1": {"T_pooled": pooled_obs, "p": p_h1, "supported": p_h1 < 0.05, "per_knockout_exploratory": per_ko},
        "H2": {"median_lodo_r": med_lodo, "n_units": len(lodo), "null_median": float(np.nanmedian(null_h2_a)), "p": p_h2, "supported": bool(med_lodo >= 0.3 and p_h2 < 0.05)},
        "H3": {"n_hallmark_genes": int(in_set.sum()), "n_other": int((~in_set).sum()), "hallmark_genes": [g for g, s in zip(a.genes, in_set, strict=True) if s], "p": p_h3, "supported": bool(p_h3 < 0.05)},
        "H4": {"n_units": len(h4), "median_diff": med_h4, "p": p_h4, "supported": bool(med_h4 > 0 and p_h4 < 0.05), "units": h4_units},
        "H5": {"units": h5_units, "median_r": float(np.median(h5_rs)) if h5_rs else None, "supported": bool(h5_rs) and all(u.get("p_perm", 1) < 0.05 for u in h5_units if "r" in u) and float(np.median(h5_rs)) > 0},
    }
    eff.to_csv(REPORTS / "step6_effects.csv", index=False)
    (REPORTS / "step6_hypotheses.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("per_knockout_exploratory", "units", "hallmark_genes")} if isinstance(v, dict) else v for k, v in report.items()}, indent=2, default=float))


if __name__ == "__main__":
    main()
