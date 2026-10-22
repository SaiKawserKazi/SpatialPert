"""A8: count-level effects, effect size and sparsity, real lesion-size distribution.

Effects are applied to raw counts of pseudo-knockout cells (binomial thinning for decreases,
Poisson boosting for increases), then cells are renormalised with the same median target as
the real data. Lesion sizes are resampled from the real knockouts' 100 µm lesion sizes.
Tests: lesion-unit permutation and within-lesion sign-flip. Writes reports/step26_realism.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from spatialpert.effects import clone_patches, fit_control_model, patch_domains, patch_permutation_test, residuals
from spatialpert.gmn import place_geometry
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import within_lesion_signflip_p
from spatialpert.xenium import read_xenium_h5

ROOT = Path(__file__).resolve().parents[1]
N = 400


def raw_counts(a) -> sparse.csr_matrix:  # noqa: ANN001
    x = read_xenium_h5(ROOT / "data/raw/GSE319123/GSM9511110_Xenium_lung_cell_feature_matrix.h5")
    m = x.matrix
    pos = {b: i for i, b in enumerate(m.barcodes)}
    col = {g: i for i, (g, fid, t) in enumerate(zip(m.feature_names, m.feature_ids, x.feature_types, strict=True)) if t == "Gene Expression" and fid.startswith("ENS")}
    rows = np.array([pos[c] for c in a.cells["cell_id"]])
    return m.counts[rows][:, [col[g] for g in a.genes]].tocsr()


def main() -> None:  # noqa: PLR0915 - linear benchmark script
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    ok = a.is_control & (nm.labels >= 0)
    model = fit_control_model(a.y[ok], nm.labels[ok], a.cov[ok])
    ci = np.where(ok)[0]
    cxy = nm.coords[ci]
    tree = cKDTree(cxy)
    counts = raw_counts(a).astype(np.float64)
    gsd = residuals(model, a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    # real lesion sizes (all four testable knockouts, 100 µm lesions)
    pool = []
    for ko in ("CIITA", "DIO3", "ESAM", "PRSS23"):
        mk = np.where(a.knockout_mask(ko))[0]
        pool += np.bincount(clone_patches(nm.coords[mk], link=100.0)).tolist()
    pool_a = np.array(pool, dtype=float)
    # Step 4 normalisation: each cell's total over ALL panel genes, target = median over all kept cells
    step4 = json.loads((ROOT / "reports/step4_domains.json").read_text())
    target = float(step4["median_total_counts_kept"])
    tot_full = np.exp(a.cells["log_total"].to_numpy())  # all-gene totals of the analysis cells
    analysis_tot = np.asarray(counts.sum(1)).ravel()
    rng = np.random.default_rng(26)

    def draw(n_les: int) -> tuple[np.ndarray, np.ndarray]:
        s = rng.choice(pool_a, size=n_les, replace=True)
        sizes = np.maximum(3, np.round(s / s.sum() * N)).astype(int)
        il, les = place_geometry(tree, cxy, sizes.tolist(), rng, return_lesions=True)
        return ci[il], les

    def tests(idx: np.ndarray, les: np.ndarray, y: sparse.csr_matrix) -> dict[str, float]:
        d = nm.leave_out_labels(idx, idx)
        keep = d >= 0
        lk = np.unique(les[keep], return_inverse=True)[1]  # contiguous lesion ids after filtering
        r = residuals(model, y[keep], d[keep], a.cov[idx][keep])
        ud = patch_domains(lk, d[keep])
        lp = patch_permutation_test(r, lk, ud, n_perm=250, seed=int(rng.integers(1 << 31))).p_value if np.unique(ud).size > 1 else np.nan
        return {"lesion": lp, "signflip": within_lesion_signflip_p(r, d[keep], lk, seed=int(rng.integers(1 << 31)))}

    def renorm(c: sparse.csr_matrix, idx: np.ndarray) -> sparse.csr_matrix:
        tot = tot_full[idx] + (np.asarray(c.sum(1)).ravel() - analysis_tot[idx])
        sc_ = np.divide(target, tot, out=np.zeros_like(tot), where=tot > 0)
        out = sparse.csr_matrix(sparse.diags(sc_) @ c)
        out.data = np.log1p(out.data)
        return out

    rows = []
    for n_les in (5, 10, 20):
        for rep in range(200):
            idx, les = draw(n_les)
            rows.append({"lesions": n_les, "delta": 0.0, "sparsity": 0.0, **tests(idx, les, renorm(counts[idx], idx))})
        for delta in (0.25, 0.5, 1.0):
            for spars in (0.05, 0.2, 1.0):
                for rep in range(100):
                    idx, les = draw(n_les)
                    c = counts[idx].toarray()
                    d = nm.leave_out_labels(idx, idx)
                    tgt = int(rng.choice(np.unique(d[d >= 0])))
                    genes = rng.random(c.shape[1]) < spars
                    if not genes.any():
                        genes[rng.integers(c.shape[1])] = True
                    lfc = np.zeros(c.shape[1])
                    lfc[genes] = delta * gsd[genes] * rng.choice([-1.0, 1.0], size=genes.sum()) * np.sqrt(1 / max(spars, 1e-9))
                    lfc = np.clip(lfc, -3, 3)
                    fold = np.exp(lfc)
                    mask = d == tgt
                    for g in np.where(genes)[0]:
                        cg = c[mask, g]
                        if fold[g] < 1:
                            c[mask, g] = rng.binomial(cg.astype(int), fold[g])
                        else:
                            c[mask, g] = cg + rng.poisson(cg * (fold[g] - 1))
                    rows.append({"lesions": n_les, "delta": delta, "sparsity": spars, **tests(idx, les, renorm(sparse.csr_matrix(c), idx))})
        print("L", n_les, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "reports/step26_realism.csv", index=False)
    summ = []
    for L in (5, 10, 20):
        nul = df[(df.lesions == L) & (df.delta == 0)]
        rec = {"lesions": L}
        for t in ("lesion", "signflip"):
            pn = nul[t].dropna()
            crit = float(np.quantile(pn, 0.05)) if len(pn) else float("nan")
            rec[t] = {"type1": float(np.mean(pn < 0.05)) if len(pn) else None, "applicable": round(len(pn) / len(nul), 3), "power": {}}
            for (dl, sp), s in df[(df.lesions == L) & (df.delta > 0)].groupby(["delta", "sparsity"]):
                pa = s[t].dropna()
                rec[t]["power"][f"d{dl}_s{sp}"] = {"raw": float(np.mean(pa < 0.05)) if len(pa) else None,
                                                    "size_adjusted": float(np.mean(pa <= crit)) if len(pa) and len(pn) else None}
        summ.append(rec)
    (ROOT / "reports/step26_realism.json").write_text(json.dumps({"lesion_size_pool": sorted(pool, reverse=True), "summary": summ}, indent=2))
    print("done")


if __name__ == "__main__":
    main()
