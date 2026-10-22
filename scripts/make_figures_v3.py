"""Nature Methods figures for manuscript v3 (review 2, figure items).

Six display items + three Extended Data figures, 89/183 mm, Arial 5-7 pt, Okabe-Ito, scale bars,
one naming scheme (nm_style.METHOD_NAME). Every plotted number is also written as source data.
Writes reports/figures/v3/{fig*,ed*}.{pdf,png} and reports/figures/v3/source_data/*.csv.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import patches
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nm_style as S  # noqa: E402

S.setup()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
F = R / "figures" / "v3"
SD = F / "source_data"
SD.mkdir(parents=True, exist_ok=True)
LS = [1, 3, 5, 10, 20]
KOS = ["CIITA", "DIO3", "ESAM", "PRSS23"]
ADOPT = json.loads((R / "step19_adoption.json").read_text())
REAL = json.loads((R / "step28_real_smallcluster.json").read_text())

# test key in step19 -> nm_style method key
T19 = {"loo_cell": "cell", "loo_patch30": "patch30", "loo_lesion": "lesion", "loo_stud_lesion": "stud_lesion",
       "loo_cr2": "cr2", "loo_wcb": "wcb", "loo_signflip": "signflip", "loo_mixedlm": "mixedlm"}


def src(name: str, df: pd.DataFrame) -> None:
    df.to_csv(SD / f"{name}.csv", index=False)


def fig2(a, nm) -> None:  # noqa: ANN001
    """Tissue, benchmark niches, knockout lesions, marker dot plot, lesions per niche."""
    lab = a.cells["label"].to_numpy().astype(object)
    xy = nm.coords
    ok = nm.labels >= 0
    fig = plt.figure(figsize=(S.DOUBLE, 120 * S.MM))
    ax = fig.add_axes([0.02, 0.45, 0.44, 0.48])
    for d in range(4):
        s = ok & (nm.labels == d)
        ax.scatter(xy[s, 0], xy[s, 1], s=0.05, c=S.NICHE_COLOR[d], rasterized=True, label=f"N{d}", lw=0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(markerscale=12, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.02), handletextpad=0.1)
    S.scale_bar(ax, 500, "500 µm")
    ax.set_title("Control-defined niches (k-means, K = 4)", loc="left", x=0.05)
    S.letter(ax, "a", x=0.0, y=0.985)
    bx = fig.add_axes([0.52, 0.45, 0.44, 0.48])
    bx.scatter(xy[ok, 0], xy[ok, 1], s=0.05, c="#DDDDDD", rasterized=True, lw=0)
    for k in KOS:
        s = lab == k
        n = REAL[k]["lesions"]
        bx.scatter(xy[s, 0], xy[s, 1], s=0.8, c=S.KO_COLOR[k], lw=0, label=f"{k} ({n} lesion{'s' if n > 1 else ''})")
    bx.set_aspect("equal")
    bx.axis("off")
    bx.legend(markerscale=5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 0.02), handletextpad=0.1)
    S.scale_bar(bx, 500, "500 µm")
    bx.set_title("Knockout cells (100 µm lesions)", loc="left")
    S.letter(bx, "b", x=0.0, y=0.985)
    # c: dot plot of niche marker genes in control cells
    ctrl = a.is_control & ok
    yc = a.y[ctrl]
    dc = nm.labels[ctrl]
    mean = np.stack([np.asarray(yc[dc == d].mean(0)).ravel() for d in range(4)])
    frac = np.stack([np.asarray((yc[dc == d] > 0).mean(0)).ravel() for d in range(4)])
    z = (mean - mean.mean(0)) / (mean.std(0) + 1e-9)
    genes = list(dict.fromkeys(int(g) for d in range(4) for g in np.argsort(-z[d])[:3]))
    cx = fig.add_axes([0.08, 0.08, 0.40, 0.24])
    rows = []
    for i, d in enumerate(range(4)):
        for j, g in enumerate(genes):
            rows.append({"niche": f"N{d}", "gene": a.genes[g], "z_mean": z[d, g], "frac_expressing": frac[d, g]})
    dfc = pd.DataFrame(rows)
    sc = cx.scatter([genes.index(next(g for g in genes if a.genes[g] == r)) for r in dfc.gene],
                    [int(n[1]) for n in dfc.niche], s=dfc.frac_expressing * 30, c=dfc.z_mean, cmap="RdBu_r",
                    vmin=-1.5, vmax=1.5, edgecolor="black", linewidth=0.2)
    cx.set_xticks(range(len(genes)), [a.genes[g] for g in genes], rotation=60, ha="right", style="italic")
    cx.set_yticks(range(4), [f"N{d} (n = {int((dc == d).sum()):,})" for d in range(4)])
    cx.set_ylim(3.6, -0.6)
    cb = fig.colorbar(sc, ax=cx, fraction=0.035, pad=0.02)
    cb.set_label("mean expression (z)")
    cb.outline.set_linewidth(0.3)
    for f_ in (0.25, 0.5, 1.0):
        cx.scatter([], [], s=f_ * 30, c="white", edgecolor="black", linewidth=0.3, label=f"{int(f_ * 100)}%")
    cx.legend(title="cells expressing", loc="upper left", bbox_to_anchor=(1.18, 1.0), title_fontsize=5.5)
    S.letter(cx, "c", x=-0.2, y=1.05)
    src("fig2c_dotplot", dfc)
    # d: lesions per niche, stacked
    dx = fig.add_axes([0.70, 0.08, 0.27, 0.27])
    rows = []
    bottom = np.zeros(len(KOS))
    for d in range(4):
        v = np.array([REAL[k]["lesions_per_niche"].get(str(d), 0) for k in KOS], dtype=float)
        dx.bar(range(len(KOS)), v, bottom=bottom, color=S.NICHE_COLOR[d], width=0.65, lw=0, label=f"N{d}")
        bottom += v
        rows += [{"knockout": k, "niche": f"N{d}", "lesions_visiting": int(x)} for k, x in zip(KOS, v, strict=True)]
    for i, k in enumerate(KOS):
        dx.text(i, bottom[i] + 0.4, f"{REAL[k]['lesions']}", ha="center", fontsize=5.5)
    dx.set_xticks(range(len(KOS)), KOS, rotation=30, ha="right")
    dx.set_ylabel("lesion-niche visits (stacked)")
    dx.legend(ncol=2, loc="upper right", handlelength=0.8, columnspacing=0.6)
    S.letter(dx, "d", x=-0.35, y=1.05)
    src("fig2d_lesions_per_niche", pd.DataFrame(rows))
    S.save(fig, F / "fig2_tissue")


def fig3() -> None:
    """Calibration-power map: every test x lesion number, type-I and size-adjusted power."""
    per = ADOPT["per_test"]
    order = ["loo_cell", "loo_patch30", "loo_wcb", "loo_mixedlm", "loo_lesion", "loo_stud_lesion", "loo_cr2", "loo_signflip"]
    t1 = np.full((len(order), len(LS)), np.nan)
    pw = np.full_like(t1, np.nan)
    rows = []
    for i, t in enumerate(order):
        for r in per[t]:
            j = LS.index(r["lesions"])
            if "type1" in r:
                t1[i, j] = r["type1"]
                pw[i, j] = r["power_size_adjusted"] if r["power_size_adjusted"] is not None else np.nan
                rows.append({"test": S.METHOD_NAME[T19[t]], "lesions": r["lesions"], "type1": r["type1"],
                             "type1_lo": r["type1_ci"][0], "type1_hi": r["type1_ci"][1], "power_size_adjusted": r["power_size_adjusted"],
                             "m_null": r["m_null"], "m_alt": r["m_alt"]})
    src("fig3_calibration_power", pd.DataFrame(rows))
    fig, axes = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    names = [S.METHOD_NAME[T19[t]] for t in order]
    for ax, mat, title, cmap, vmax in ((axes[0], t1, "Type-I error (300 null draws)", "Reds", 0.25),
                                       (axes[1], pw, "Size-adjusted power (150 draws)", "Blues", 1.0)):
        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isnan(v):
                    ax.add_patch(patches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", ec="#AAAAAA", lw=0))
                    continue
                bad = ax is axes[0] and v > 0.06
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5,
                        color="white" if v > 0.6 * vmax else "black", fontweight="bold" if bad else "normal")
        ax.set_xticks(range(len(LS)), LS)
        ax.set_xlabel("lesions per knockout (400 cells)")
        ax.set_yticks(range(len(order)), names if ax is axes[0] else [""] * len(order))
        ax.set_title(title, loc="left")
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.outline.set_linewidth(0.3)
        ax.axhline(3.5, color="black", lw=0.6)
    axes[0].text(4.55, 3.3, "inflated", ha="right", va="bottom", fontsize=5, style="italic")
    axes[0].text(4.55, 3.7, "calibrated (≤ 0.06)", ha="right", va="top", fontsize=5, style="italic")
    S.letter(axes[0], "a", x=-0.75)
    S.letter(axes[1], "b", x=-0.08)
    # c: circular vs leave-out niche labels (A2)
    ax = axes[2]
    for t, key, mk in (("circ_lesion", "lesion", "o"), ("circ_patch30", "patch30", "s")):
        c = [next((r["type1"] for r in per[t] if r["lesions"] == L and "type1" in r), np.nan) for L in LS]
        lo = [next((r["type1"] for r in per["loo_" + t[5:]] if r["lesions"] == L and "type1" in r), np.nan) for L in LS]
        ax.plot(LS, c, marker=mk, ls="--", c=S.METHOD_COLOR[key], mfc="white", label=f"{S.METHOD_NAME[key]}, circular")
        ax.plot(LS, lo, marker=mk, ls="-", c=S.METHOD_COLOR[key], label=f"{S.METHOD_NAME[key]}, leave-out")
    ax.axhline(0.05, c="black", lw=0.4, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(LS, LS)
    ax.minorticks_off()
    ax.set_xlabel("lesions per knockout")
    ax.set_ylabel("type-I error")
    ax.set_title("Knockout cells define their niche", loc="left")
    ax.legend(loc="upper left", fontsize=5)
    S.letter(ax, "c", x=-0.3)
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.2, top=0.88, wspace=0.45)
    S.save(fig, F / "fig3_calibration_power")


def fig4() -> None:
    """Design: power curve with CIs, formula failure, count-level realism, Perturb-map lesion counts."""
    per = ADOPT["per_test"]
    d27 = json.loads((R / "step27_design.json").read_text())
    d26 = json.loads((R / "step26_realism.json").read_text())
    d25 = json.loads((R / "step25_perturbmap_lesions.json").read_text())
    fig, axes = plt.subplots(1, 4, figsize=(S.DOUBLE, 50 * S.MM))
    ax = axes[0]
    rows = []
    for t in ("loo_lesion", "loo_stud_lesion", "loo_cr2", "loo_signflip"):
        rr = [r for r in per[t] if r.get("power_size_adjusted") is not None]
        x = [r["lesions"] for r in rr]
        y = [r["power_size_adjusted"] for r in rr]
        ax.fill_between(x, [r["power_ci"][0] for r in rr], [r["power_ci"][1] for r in rr], color=S.METHOD_COLOR[T19[t]], alpha=0.15, lw=0)
        ax.plot(x, y, marker="o", c=S.METHOD_COLOR[T19[t]], label=S.METHOD_NAME[T19[t]])
        rows += [{"test": S.METHOD_NAME[T19[t]], "lesions": r["lesions"], "power": r["power_size_adjusted"],
                  "lo": r["power_ci"][0], "hi": r["power_ci"][1]} for r in rr]
    src("fig4a_power_curve", pd.DataFrame(rows))
    ax.axvspan(2.8, 10, color=S.OI["grey"], alpha=0.12, lw=0)
    ax.text(4.1, 0.3, "single-section\nscreens", ha="center", fontsize=5, color="#555555")
    ax.axhline(0.8, c="black", lw=0.4, ls=":")
    ax.set_title("Calibrated tests", loc="left")
    ax.set_ylabel("size-adjusted power")
    ax.legend(loc="upper left", fontsize=5)
    ax = axes[1]
    ax.plot(d27["lesions"], d27["simulated_size_adjusted_power"], marker="o", c=S.METHOD_COLOR["lesion"], label="simulated")
    ax.plot(d27["lesions"], d27["formula_power_efficiency1"], marker="^", ls="--", c=S.OI["grey"], label="design effect, eff. 1")
    ax.plot(d27["lesions"], d27["formula_power_fitted"], marker="v", ls="--", c=S.OI["black"],
            label=f"design effect, eff. {d27['fitted_efficiency']:.2f}")
    ax.set_title(f"Design effect (RMSE {d27['rmse_fitted']:.2f})", loc="left")
    ax.legend(loc="upper left", fontsize=5)
    src("fig4b_formula", pd.DataFrame({k: d27[k] for k in ("lesions", "simulated_size_adjusted_power", "formula_power_efficiency1", "formula_power_fitted")}))
    ax = axes[2]
    rows = []
    for dl, mk in (("d0.5", "o"), ("d1.0", "s")):
        x = [r["lesions"] for r in d26["summary"]]
        y = [r["lesion"]["power"][f"{dl}_s1.0"]["size_adjusted"] for r in d26["summary"]]
        t = [r["lesion"]["type1"] for r in d26["summary"]]
        ax.plot(x, y, marker=mk, c=S.METHOD_COLOR["lesion"], label=f"power, δ = {dl[1:]}")
        rows += [{"lesions": a_, "delta": dl[1:], "power_size_adjusted": b_, "type1": c_} for a_, b_, c_ in zip(x, y, t, strict=True)]
    ax.plot(x, t, marker="x", c=S.OI["vermillion"], ls=":", label="type-I error")
    ax.set_title("Count-level realism", loc="left")
    ax.legend(loc="upper left", fontsize=5)
    src("fig4c_realism", pd.DataFrame(rows))
    for ax in axes[:3]:
        ax.set_xscale("log")
        ax.set_xticks([3, 5, 10, 20], [3, 5, 10, 20])
        ax.minorticks_off()
        ax.set_xlabel("lesions per knockout")
        ax.set_ylim(-0.02, 1.02)
    ax = axes[3]
    n = np.array(sorted((v["lesions"] for v in d25["per_code"].values()), reverse=True))
    ax.bar(np.arange(n.size), n, color=S.OI["grey"], width=0.8, lw=0)
    for thr, c in ((10, S.OI["orange"]), (20, S.OI["blue"])):
        ax.axhline(thr, c=c, lw=0.6, ls="--")
        ax.text(n.size - 1, thr * 1.08, f"{int((n >= thr).sum())} codes ≥ {thr}", ha="right", fontsize=5, color=c, bbox={"fc": "white", "ec": "none", "pad": 0.5})
    ax.set_yscale("log")
    ax.set_xlabel("Pro-Code (ranked)")
    ax.set_ylabel("lesions (8 mice)")
    ax.set_title("Perturb-map lesion counts", loc="left")
    src("fig4d_perturbmap", pd.DataFrame({"procode": list(d25["per_code"]), "lesions": [v["lesions"] for v in d25["per_code"].values()],
                                          "mice": [v["mice"] for v in d25["per_code"].values()]}))
    for i, ax in enumerate(axes):
        S.letter(ax, "abcd"[i], x=-0.3)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.2, top=0.86, wspace=0.42)
    S.save(fig, F / "fig4_design")


def fig5() -> None:
    """Generalisation: second tissue (brain) and four niche callers."""
    br = json.loads((R / "step22_brain.json").read_text())
    ca = json.loads((R / "step24_callers.json").read_text()) | json.loads((R / "step24_callers_cellcharter.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(S.DOUBLE, 58 * S.MM), gridspec_kw={"width_ratios": [1, 1.5]})
    ax = axes[0]
    rows = []
    for key in ("cell", "patch30", "lesion", "wcb"):
        x = [r["lesions"] for r in br["summary"] if r[key]["type1"] is not None and r[key]["applicable"] >= 0.1]
        y = [r[key]["type1"] for r in br["summary"] if r[key]["type1"] is not None and r[key]["applicable"] >= 0.1]
        ax.plot(x, y, marker="o", c=S.METHOD_COLOR[key], label=S.METHOD_NAME[key])
        rows += [{"test": S.METHOD_NAME[key], "lesions": r["lesions"], "type1": r[key]["type1"], "applicable": r[key]["applicable"]}
                 for r in br["summary"]]
    src("fig5a_brain", pd.DataFrame(rows))
    ax.axhline(0.05, c="black", lw=0.4, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(LS, LS)
    ax.minorticks_off()
    ax.set_xlabel("lesions per knockout")
    ax.set_ylabel("type-I error")
    ax.set_title("Mouse brain (10 µm bins)", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=5)
    ax = axes[1]
    callers = [("kmeans", "k-means"), ("banksy", "BANKSY"), ("spagcn", "SpaGCN"), ("cellcharter", "CellCharter")]
    tests = ["lesion", "wcb", "signflip"]
    Lc = [3, 5, 10, 20]
    mat = np.full((len(callers) * len(tests), len(Lc)), np.nan)
    ylabs, rows = [], []
    for i, (ck, cn) in enumerate(callers):
        for j, t in enumerate(tests):
            ylabs.append(f"{cn} · {S.METHOD_NAME[t]}")
            for r in ca[ck]["summary"]:
                v = r[t]["type1"]
                mat[i * len(tests) + j, Lc.index(r["lesions"])] = np.nan if v is None else v
                rows.append({"caller": cn, "test": S.METHOD_NAME[t], "lesions": r["lesions"], "type1": v, "applicable": r[t]["applicable"]})
    src("fig5b_callers", pd.DataFrame(rows))
    im = ax.imshow(mat, cmap="Reds", vmin=0, vmax=0.25, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                ax.add_patch(patches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", ec="#AAAAAA", lw=0))
            else:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5, color="white" if v > 0.15 else "black",
                        fontweight="bold" if v > 0.06 else "normal")
    for i in range(1, len(callers)):
        ax.axhline(i * len(tests) - 0.5, c="black", lw=0.5)
    ax.set_yticks(range(len(ylabs)), ylabs)
    ax.set_xticks(range(len(Lc)), Lc)
    ax.set_xlabel("lesions per knockout")
    ax.set_title("Type-I error by niche caller (tumour)", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.outline.set_linewidth(0.3)
    S.letter(axes[0], "a", x=-0.25)
    S.letter(axes[1], "b", x=-0.75)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.18, top=0.88, wspace=0.75)
    S.save(fig, F / "fig5_generalisation")


def fig6() -> None:
    """Real data: p-value by test, untestable hatched, and lesions per niche."""
    cols = [("cell", "cell_p"), ("patch30", "patch30_p"), ("lesion", "lesion_p"), ("stud_lesion", "stud_lesion_p"),
            ("cr2", "cr2_min_q"), ("wcb", "wcb_p"), ("signflip", "signflip_p"), ("mixedlm", "mixedlm_min_q")]
    mat = np.full((len(KOS), len(cols)), np.nan)
    rows = []
    for i, k in enumerate(KOS):
        for j, (m, key) in enumerate(cols):
            v = REAL[k].get(key)
            v = None if v is None or v != v else float(v)
            mat[i, j] = np.nan if v is None else -np.log10(max(v, 1e-3))
            rows.append({"knockout": k, "test": S.METHOD_NAME[m], "p_or_min_q": v, "testable": v is not None})
    src("fig6a_real_p", pd.DataFrame(rows))
    fig, axes = plt.subplots(1, 2, figsize=(S.DOUBLE, 55 * S.MM), gridspec_kw={"width_ratios": [2.2, 1]})
    ax = axes[0]
    im = ax.imshow(mat, cmap="Reds", vmin=0, vmax=3, aspect="auto")
    for i in range(len(KOS)):
        for j in range(len(cols)):
            v = mat[i, j]
            if np.isnan(v):
                ax.add_patch(patches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", ec="#AAAAAA", lw=0))
                ax.text(j, i, "untestable", ha="center", va="center", fontsize=4.5, color="#555555",
                        bbox={"fc": "white", "ec": "none", "pad": 0.3})
                continue
            p = 10 ** -v
            txt = "p ≤ 0.001" if p <= 0.001 + 1e-12 else f"{p:.3f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5, color="white" if v > 2 else "black")
    ax.set_xticks(range(len(cols)), [S.METHOD_NAME[m] for m, _ in cols], rotation=30, ha="right")
    ax.set_yticks(range(len(KOS)), KOS)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("−log10 p (CR2, mixed model: min q)")
    cb.outline.set_linewidth(0.3)
    ax.axvline(0.5, c="black", lw=0.6)
    ax.set_title("Only the cell-level test rejects; lesion-aware tests are mostly untestable", loc="left")
    ax = axes[1]
    bottom = np.zeros(len(KOS))
    for d in range(4):
        v = np.array([REAL[k]["cells_per_niche"].get(str(d), 0) for k in KOS], dtype=float)
        ax.barh(range(len(KOS)), v, left=bottom, color=S.NICHE_COLOR[d], height=0.65, lw=0, label=f"N{d}")
        bottom += v
    for i, k in enumerate(KOS):
        per = REAL[k]["lesions_per_niche"]
        ax.text(bottom[i] + 8, i, "/".join(str(per.get(str(d), 0)) for d in range(4)), va="center", fontsize=5)
    ax.set_yticks(range(len(KOS)), KOS)
    ax.invert_yaxis()
    ax.set_xlabel("cells per niche (lesions N0/N1/N2/N3)")
    ax.set_xlim(0, bottom.max() * 1.45)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=0.8, columnspacing=0.6)
    S.letter(axes[0], "a", x=-0.12)
    S.letter(axes[1], "b", x=-0.25, y=1.12)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.3, top=0.85, wspace=0.3)
    S.save(fig, F / "fig6_real_data")


def ed1() -> None:
    """QQ plots of null p-values with 95% pointwise beta bands."""
    df = pd.read_csv(R / "step19_stage1b.csv")
    tests = ["loo_cell", "loo_patch30", "loo_lesion", "loo_cr2", "loo_wcb", "loo_signflip"]
    fig, axes = plt.subplots(2, 3, figsize=(S.DOUBLE, 100 * S.MM))
    rows = []
    for ax, t in zip(axes.ravel(), tests, strict=True):
        for L, c in ((5, S.OI["sky"]), (10, S.OI["blue"]), (20, S.OI["black"])):
            p = np.sort(df[(df.lesions == L) & (df["lambda"] == 0)][t].dropna().to_numpy())
            if p.size < 10:
                continue
            n = p.size
            i = np.arange(1, n + 1)
            e = i / (n + 1)
            lo, hi = stats.beta.ppf(0.025, i, n - i + 1), stats.beta.ppf(0.975, i, n - i + 1)
            if L == 20:
                ax.fill_between(-np.log10(e), -np.log10(lo), -np.log10(hi), color=S.OI["grey"], alpha=0.25, lw=0)
            ax.plot(-np.log10(e), -np.log10(np.maximum(p, 1e-4)), ".", ms=1.5, c=c, label=f"{L} lesions (n = {n})")
            rows += [{"test": S.METHOD_NAME[T19[t]], "lesions": L, "expected": a_, "observed": b_} for a_, b_ in zip(e, p, strict=True)]
        m = 2.6
        ax.plot([0, m], [0, m], c="black", lw=0.4, ls=":")
        ax.set_xlim(0, m)
        ax.set_ylim(0, 4.1)
        ax.set_title(S.METHOD_NAME[T19[t]], loc="left")
        ax.set_xlabel("expected −log10 p")
        ax.set_ylabel("observed −log10 p")
    axes[0, 0].legend(loc="upper left", fontsize=5, markerscale=3)
    for k, ax in enumerate(axes.ravel()):
        S.letter(ax, "abcdef"[k], x=-0.3)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.09, top=0.94, wspace=0.4, hspace=0.6)
    src("ed1_qq", pd.DataFrame(rows))
    S.save(fig, F / "ed1_qq")


def ed2() -> None:
    """ICC by lesion size and runtime."""
    icc = json.loads((R / "step23_icc.json").read_text())
    rt = json.loads((R / "step20_runtime.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(S.DOUBLE, 50 * S.MM))
    ax = axes[0]
    m = [int(k) for k in icc]
    ax.errorbar(m, [icc[k]["icc_median"] for k in icc], yerr=np.array([[icc[k]["icc_median"] - icc[k]["icc_iqr"][0],
                icc[k]["icc_iqr"][1] - icc[k]["icc_median"]] for k in icc]).T, marker="o", c=S.OI["blue"], capsize=1.5, label="summed score")
    ax.plot(m, [icc[k]["per_gene_icc_median"] for k in icc], marker="s", c=S.OI["orange"], label="per gene (median)")
    ax.set_xscale("log")
    ax.set_xticks(m, m)
    ax.minorticks_off()
    ax.set_xlabel("cells per lesion")
    ax.set_ylabel("within-lesion ICC")
    ax.set_title("Within-lesion correlation is small", loc="left")
    ax.legend(loc="upper right", fontsize=5)
    src("ed2a_icc", pd.DataFrame([{"cells_per_lesion": int(k), **{a: b for a, b in v.items() if a != "icc_iqr"}} for k, v in icc.items()]))
    ax = axes[1]
    items = [("load + normalise\n(352,913 cells)", rt["load_and_normalise_full_section"]["seconds"]),
             ("fit niches", rt["niche_model_fit"]["seconds"]),
             ("patch test /KO", rt["location_dependence_patch_1000perm_4KO"]["seconds_per_knockout"]),
             ("lesion test /KO", rt["location_dependence_lesion_1000perm_4KO"]["seconds_per_knockout"]),
             ("wild bootstrap\n(999 draws)", rt["wild_cluster_bootstrap_999_pseudoKO_10lesions_400cells"]["seconds"]),
             ("sign-flip\n(2,000 flips)", rt["signflip_2000_pseudoKO_10lesions_400cells"]["seconds"])]
    ax.barh(range(len(items)), [v for _, v in items], color=S.OI["grey"], height=0.6)
    for i, (_, v) in enumerate(items):
        ax.text(v * 1.15, i, f"{v:.3g} s", va="center", fontsize=5)
    ax.set_xscale("log")
    ax.set_yticks(range(len(items)), [n for n, _ in items])
    ax.invert_yaxis()
    ax.set_xlabel("seconds (4-core laptop, peak RSS 0.63 GB)")
    ax.set_title("Runtime", loc="left")
    src("ed2b_runtime", pd.DataFrame(items, columns=["step", "seconds"]))
    S.letter(axes[0], "a", x=-0.15)
    S.letter(axes[1], "b", x=-0.45)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.2, top=0.86, wspace=0.6)
    S.save(fig, F / "ed2_icc_runtime")


def ed3(a) -> None:  # noqa: ANN001
    """Niche maps from the four callers on the same control cells."""
    pos = a.cells.set_index("cell_id")
    fig, axes = plt.subplots(1, 4, figsize=(S.DOUBLE, 45 * S.MM))
    fig.subplots_adjust(left=0.02, right=0.99, bottom=0.02, top=0.9, wspace=0.08)
    from spatialpert.pipeline import NicheModel

    nm = NicheModel.fit(a)
    ids = a.cells["cell_id"].to_numpy()
    maps = [("k-means", pd.Series(nm.labels, index=ids))]
    for f, n in (("banksy", "BANKSY"), ("spagcn", "SpaGCN"), ("cellcharter", "CellCharter")):
        t = pd.read_csv(ROOT / f"data/interim/a6/{f}_labels.csv").set_index("cell_id").iloc[:, 0]
        maps.append((n, t))
    for k, (ax, (n, lab)) in enumerate(zip(axes, maps, strict=True)):
        lab = lab[lab.index.isin(pos.index) & (lab >= 0)]
        xy = pos.loc[lab.index, ["x", "y"]].to_numpy()
        cols = np.array(S.NICHE_COLOR)[lab.to_numpy() % len(S.NICHE_COLOR)]
        ax.scatter(xy[:, 0], xy[:, 1], s=0.03, c=cols, rasterized=True, lw=0)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{n} ({lab.nunique()} niches)", loc="left", x=0.08)
        S.scale_bar(ax, 500, "500 µm")
        S.letter(ax, "abcd"[k], x=-0.05, y=1.02)
    S.save(fig, F / "ed3_callers")


def ed4(a, nm) -> None:  # noqa: ANN001
    """Seed placement at 20 lesions vs control-cell density, and host fraction by niche."""
    seeds = pd.read_csv(R / "step19_seed_centroids_L20.csv")
    ok = nm.labels >= 0
    ctrl = a.is_control & ok
    xy = nm.coords
    fig, axes = plt.subplots(1, 3, figsize=(S.DOUBLE, 55 * S.MM), gridspec_kw={"width_ratios": [1.2, 1, 0.9]})
    ax = axes[0]
    ax.scatter(xy[ctrl, 0], xy[ctrl, 1], s=0.03, c="#CCCCCC", rasterized=True, lw=0)
    ax.scatter(seeds.x, seeds.y, s=0.6, c=S.OI["vermillion"], lw=0, rasterized=True)
    ax.set_aspect("equal")
    ax.axis("off")
    S.scale_bar(ax, 1000, "1 mm")
    ax.set_title("Lesion centroids, 20 lesions (60 null draws)", loc="left", x=0.08)
    # b: radial distribution - distance to section centroid, seeds vs control cells
    ax = axes[1]
    c0 = xy[ctrl].mean(0)
    dc = np.hypot(*(xy[ctrl] - c0).T)
    ds = np.hypot(*(seeds[["x", "y"]].to_numpy() - c0).T)
    bins = np.linspace(0, max(dc.max(), ds.max()), 30)
    ax.hist(dc, bins=bins, density=True, histtype="step", color=S.OI["grey"], label="control cells")
    ax.hist(ds, bins=bins, density=True, histtype="step", color=S.OI["vermillion"], label="lesion centroids")
    ks = stats.ks_2samp(ds, dc)
    ax.set_xlabel("distance to section centre (µm)")
    ax.set_ylabel("density")
    ax.set_title(f"Placement vs tissue (KS D = {ks.statistic:.2f})", loc="left")
    ax.legend(loc="upper right", fontsize=5)
    src("ed4b_seed_radial", pd.DataFrame({"set": ["seed"] * ds.size + ["control"] * dc.size, "distance_um": np.r_[ds, dc]}))
    # c: host fraction by niche (violins)
    ax = axes[2]
    hf = a.cells["host_frac50"].to_numpy() if "host_frac50" in a.cells else a.cov[:, 1]
    data = [hf[ctrl & (nm.labels == d)] for d in range(4)]
    vp = ax.violinplot(data, showmedians=True, widths=0.8)
    for d, b in enumerate(vp["bodies"]):
        b.set_facecolor(S.NICHE_COLOR[d])
        b.set_alpha(0.8)
        b.set_linewidth(0)
    for k in ("cbars", "cmins", "cmaxes", "cmedians"):
        vp[k].set_color("black")
        vp[k].set_linewidth(0.5)
    ax.set_xticks(range(1, 5), [f"N{d}" for d in range(4)])
    ax.set_ylabel("host-cell fraction within 50 µm")
    ax.set_title("Host fraction by niche", loc="left")
    src("ed4c_host_fraction", pd.DataFrame({"niche": np.concatenate([[f"N{d}"] * len(x) for d, x in enumerate(data)]), "host_frac50": np.concatenate(data)}))
    for k, ax in enumerate(axes):
        S.letter(ax, "abc"[k], x=-0.2 if k else 0.0, y=1.02)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.2, top=0.88, wspace=0.4)
    S.save(fig, F / "ed4_placement_niches")


def main() -> None:
    sys.path.insert(0, str(ROOT / "src"))
    from spatialpert.pipeline import NicheModel, load_xenium_analysis

    fig3()
    fig4()
    fig5()
    fig6()
    ed1()
    ed2()
    a = load_xenium_analysis(ROOT)
    fig2(a, NicheModel.fit(a))
    ed3(a)
    ed4(a, NicheModel.fit(a))
    print("written", F)


if __name__ == "__main__":
    main()
