"""Review Stage 1 (AMENDMENT-003): seven tests on identical semi-synthetic draws.

lambda = 0: 500 replicates per size; lambda > 0: 200. 250 permutations. Records patch counts.
Writes reports/step10_stage1.csv and reports/step10_stage1.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from spatialpert.baselines import (
    cluster_robust_any,
    patch_means,
    pseudobulk_anova_any,
    studentized_permutation_p,
)
from spatialpert.benchmark import anova_any_gene, cell_permutation_p, inject_effect, scale_sizes
from spatialpert.effects import (
    clone_patches,
    fit_control_model,
    matched_pseudo_patches,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.pipeline import load_xenium_analysis, tile_ids

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
TESTABLE = ["CIITA", "DIO3", "ESAM", "PRSS23"]
NS = [50, 100, 200, 400, 800]
LAMBDAS = [0.0, 0.25, 0.5, 1.0]
REPS = {0.0: 500, 0.25: 200, 0.5: 200, 1.0: 200}
N_PERM = 250
DELTA = 0.5
A = 0.05
METHODS = ["patch", "studentized", "pseudobulk", "sandwich", "tile", "cell", "anova"]


def main() -> None:
    a = load_xenium_analysis(ROOT)
    ctrl = a.is_control
    dom = a.cells["domain"].to_numpy()
    coords = a.cells[["x", "y"]].to_numpy()
    model = fit_control_model(a.y[ctrl], dom[ctrl], a.cov[ctrl])
    ctrl_idx = np.where(ctrl)[0]
    sub = ctrl_idx[:20000]
    gene_sd = residuals(model, a.y[sub], dom[sub], a.cov[sub]).std(axis=0)
    templates = []
    for ko in TESTABLE:
        mk = a.knockout_mask(ko)
        c = pd.Series(dom[mk]).value_counts()
        sel = mk & np.isin(dom, c.index[c >= 30].to_numpy())
        p = clone_patches(coords[sel])
        templates.append((np.bincount(p), patch_domains(p, dom[sel])))

    rows = []
    rng = np.random.default_rng(10)
    for n in NS:
        for lam in LAMBDAS:
            for rep in range(REPS[lam]):
                sizes0, pdom0 = templates[rep % len(templates)]
                idx_l, pid = matched_pseudo_patches(coords[ctrl_idx], dom[ctrl_idx], scale_sizes(sizes0, n), pdom0, rng)
                idx = ctrl_idx[idx_l]
                pd_ = patch_domains(pid, dom[idx])
                if np.unique(pd_).size < 2:
                    continue
                cdom = pd_[pid]
                r = residuals(model, a.y[idx], dom[idx], a.cov[idx])
                r = inject_effect(r, cdom, gene_sd, delta=DELTA, lam=lam, target_domain=int(np.unique(pd_)[0]), rng=rng)
                seed = int(rng.integers(1 << 31))
                means, _, g = patch_means(r, pid, pd_)
                tid = tile_ids(coords[idx], 250.0)
                tdom = patch_domains(tid, cdom)
                sp = studentized_permutation_p(means, g, n_perm=N_PERM, seed=seed)
                pb = pseudobulk_anova_any(means, g)
                sw = cluster_robust_any(r, cdom, pid)
                _, per_niche = np.unique(pd_, return_counts=True)
                rows.append({
                    "n_target": n, "lambda": lam, "rep": rep, "n_cells": int(r.shape[0]),
                    "n_patches": int(pd_.size), "min_patches_per_niche": int(per_niche.min()),
                    "patch": patch_permutation_test(r, pid, pd_, n_perm=N_PERM, seed=seed).p_value < A,
                    "studentized": None if np.isnan(sp) else sp < A,
                    "pseudobulk": pb,
                    "sandwich": sw,
                    "tile": patch_permutation_test(r, tid, tdom, n_perm=N_PERM, seed=seed).p_value < A if np.unique(tdom).size >= 2 else None,
                    "cell": cell_permutation_p(r, cdom, n_perm=N_PERM, seed=seed) < A,
                    "anova": anova_any_gene(r, cdom),
                })
            print("n", n, "lambda", lam, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(R / "step10_stage1.csv", index=False)

    def rate(s: pd.Series) -> dict[str, float]:
        v = s.dropna().astype(bool)
        k, m = int(v.sum()), int(len(v))
        return {"rate": k / m if m else float("nan"), "k": k, "m": m, "na": int(s.isna().sum())}

    summary = []
    for (n, lam), s in df.groupby(["n_target", "lambda"]):
        summary.append({"n": int(n), "lambda": float(lam), **{m: rate(s[m]) for m in METHODS}})
    t1 = {m: max(r[m]["rate"] for r in summary if r["lambda"] == 0) for m in METHODS}
    pw = {m: float(np.mean([r[m]["rate"] for r in summary if r["lambda"] > 0])) for m in METHODS}
    qual = [m for m in METHODS if t1[m] <= 0.06]
    adopted = max(qual, key=lambda m: pw[m]) if qual else None
    by_patches = df[df["lambda"] == 0].groupby(pd.cut(df.loc[df["lambda"] == 0, "n_patches"], [0, 10, 20, 30, 40, 80]), observed=True)["patch"].agg(["mean", "size"]).reset_index()
    by_patches["n_patches"] = by_patches["n_patches"].astype(str)
    out = {
        "settings": {"ns": NS, "lambdas": LAMBDAS, "reps": {str(k): v for k, v in REPS.items()}, "n_perm": N_PERM, "delta": DELTA},
        "summary": summary,
        "max_type1": t1,
        "mean_power": pw,
        "qualifying": qual,
        "adopted": adopted,
        "patch_type1_by_n_patches": by_patches.to_dict(orient="records"),
    }
    (R / "step10_stage1.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ("max_type1", "mean_power", "qualifying", "adopted")}, indent=2))


if __name__ == "__main__":
    main()
