"""Step 5 (protocol §6): PC1 estimator check on Perturb-map, PC2 calibration on Xenium.

Writes reports/step5_controls.json. No H1-H5 statistic is computed here.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from spatialpert.domains import normalise_log1p
from spatialpert.effects import (
    clone_patches,
    fit_control_model,
    matched_pseudo_patches,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.gate import parse_knockout
from spatialpert.pipeline import load_xenium_analysis
from spatialpert.visium import read_10x_h5, read_spot_annotation

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "step5_controls.json"
PMAP = ROOT / "data" / "raw" / "GSE193460"
SECTIONS = {"KP_1": "GSM5808054", "KP_2": "GSM5808055", "KP_3": "GSM5808056", "KP_4": "GSM5808057"}
TESTABLE = ["CIITA", "DIO3", "ESAM", "PRSS23"]
N_PSEUDO = 200
N_PERM = 1000
ALPHA = 0.05
ARRAYS: dict[str, object] = {}


def _median_pairwise_r(vecs: list[np.ndarray]) -> float:
    rs = [
        float(np.corrcoef(vecs[i], vecs[j])[0, 1])
        for i in range(len(vecs))
        for j in range(i + 1, len(vecs))
    ]
    return float(np.median(rs))


def pc1() -> dict[str, object]:
    """Tgfbr2 (and Ifngr2) vs unassigned KP tumour spots; agreement across sections."""
    mats, anns = {}, {}
    for sec, gsm in SECTIONS.items():
        m = read_10x_h5(PMAP / f"{gsm}_{sec}_filtered_feature_bc_matrix.h5")
        ann = read_spot_annotation(PMAP / f"{gsm}_{sec}_spot_annotation.csv.gz")
        pos = {b: i for i, b in enumerate(m.barcodes)}
        ann = ann[ann["barcode"].isin(pos)].reset_index(drop=True)
        rows = np.array([pos[b] for b in ann["barcode"]])
        mats[sec] = (m.counts[rows].tocsr(), m.feature_names)
        anns[sec] = ann

    def label(s: object) -> str | None:
        if not isinstance(s, str):
            return None
        g = parse_knockout(s)
        if g is not None:
            return g
        return "KP" if s.startswith("KP_") else None

    # Analysis genes: mean normalised expression >= 0.1 in KP control spots, all sections.
    names = mats["KP_1"][1]
    logs, labs, secs = [], [], []
    for sec in SECTIONS:
        counts, fn = mats[sec]
        assert fn == names
        lab = anns[sec]["phenotypes"].map(label).to_numpy()
        keep = pd.notna(lab)
        logs.append(normalise_log1p(counts[keep]))
        labs.append(lab[keep])
        secs.append(np.full(int(keep.sum()), sec))
    y = sparse.vstack(logs).tocsr()
    lab = np.concatenate(labs).astype(str)
    sec = np.concatenate(secs)
    ctrl = lab == "KP"
    norm = y[ctrl].copy()
    norm.data = np.expm1(norm.data)
    genes = np.where(np.asarray(norm.mean(axis=0)).ravel() >= 0.1)[0]
    y = y[:, genes].tocsr()
    totals = []  # log total counts of the kept spots, in the same order as y
    for s in SECTIONS:
        counts, _ = mats[s]
        k = pd.notna(anns[s]["phenotypes"].map(label).to_numpy())
        totals.append(np.log(np.asarray(counts[k].sum(axis=1)).ravel()))
    log_total = np.concatenate(totals)
    cov = log_total[:, None]

    model = fit_control_model(y[ctrl], sec[ctrl], cov[ctrl])
    rng = np.random.default_rng(0)
    out: dict[str, object] = {"n_genes": int(genes.size), "n_kp_spots": int(ctrl.sum())}
    for ko in ("Tgfbr2", "Ifngr2"):
        ko_secs = [s for s in SECTIONS if ((sec == s) & (lab == ko)).sum() >= 30]
        mask = lab == ko
        r = residuals(model, y[mask], sec[mask], cov[mask])
        vecs = [r[sec[mask] == s].mean(axis=0) for s in ko_secs]
        obs = _median_pairwise_r(vecs) if len(vecs) >= 2 else float("nan")
        # null: shuffle KO/KP labels within each section, same KO counts
        null = []
        tum = (lab == ko) | ctrl
        for _ in range(N_PERM):
            vv = []
            for s in ko_secs:
                idx = np.where(tum & (sec == s))[0]
                n_ko = int(((lab == ko) & (sec == s)).sum())
                pick = rng.choice(idx, size=n_ko, replace=False)
                vv.append(residuals(model, y[pick], sec[pick], cov[pick]).mean(axis=0))
            null.append(_median_pairwise_r(vv))
        null_a = np.array(null)
        ARRAYS[f'PC1_null_{ko}'] = null_a.tolist()
        p = float((1 + np.sum(null_a >= obs)) / (1 + N_PERM))
        out[ko] = {
            "sections": ko_secs,
            "spots_per_section": {s: int(((lab == ko) & (sec == s)).sum()) for s in ko_secs},
            "median_pairwise_r": obs,
            "null_median": float(np.median(null_a)),
            "p_value": p,
            "pass": bool(obs >= 0.3 and p < 0.05),
        }
    return out


def pc2() -> dict[str, object]:
    """Pseudo-knockout false-positive rate of the H1 patch-permutation test."""
    a = load_xenium_analysis(ROOT)
    ctrl = a.is_control
    model = fit_control_model(a.y[ctrl], a.cells["domain"].to_numpy()[ctrl], a.cov[ctrl])
    coords = a.cells[["x", "y"]].to_numpy()
    dom = a.cells["domain"].to_numpy()
    ctrl_idx = np.where(ctrl)[0]
    rng = np.random.default_rng(0)

    templates = []
    for ko in TESTABLE:
        mk = a.knockout_mask(ko)
        counts = pd.Series(dom[mk]).value_counts()
        ok_d = counts.index[counts >= 30].to_numpy()
        sel = mk & np.isin(dom, ok_d)
        patch = clone_patches(coords[sel])
        pdom = patch_domains(patch, dom[sel])
        sizes = np.bincount(patch)
        templates.append((ko, sizes, pdom))

    pvals, used = [], []
    for i in range(N_PSEUDO):
        ko, sizes, pdom = templates[i % len(templates)]
        idx_local, pid = matched_pseudo_patches(coords[ctrl_idx], dom[ctrl_idx], sizes, pdom, rng)
        idx = ctrl_idx[idx_local]
        r = residuals(model, a.y[idx], dom[idx], a.cov[idx])
        pd_ = patch_domains(pid, dom[idx])
        if np.unique(pd_).size < 2:
            continue
        res = patch_permutation_test(r, pid, pd_, n_perm=N_PERM, seed=i)
        pvals.append(res.p_value)
        used.append(ko)
    p = np.array(pvals)
    ARRAYS['PC2_pvalues'] = p.tolist()
    fpr = float(np.mean(p < ALPHA))
    return {
        "n_pseudo": int(p.size),
        "templates": {k: int(sum(u == k for u in used)) for k in TESTABLE},
        "fpr_at_0.05": fpr,
        "threshold": 0.08,
        "pass": bool(fpr <= 0.08),
        "p_quantiles": {q: float(np.quantile(p, q)) for q in (0.05, 0.25, 0.5, 0.75)},
    }


def main() -> None:
    report = {"PC1": pc1(), "PC2": pc2()}
    REPORT.write_text(json.dumps(report, indent=2))
    (ROOT / 'reports' / 'step5_arrays.json').write_text(json.dumps(ARRAYS))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
