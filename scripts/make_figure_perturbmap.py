"""Perturb-map figure (AMENDMENT-010): lesions per niche, niche dependence Psi with CR2 intervals,
main effects (sanity check; mouse-level t-test, exploratory) and permutation vs calibrated p.
Writes reports/figures/v4/fig_perturbmap.{pdf,png} and source data; reports/step36b_main_effects.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nm_style as S  # noqa: E402

S.setup()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
F = R / "figures" / "v4"
SD = F / "source_data"
SD.mkdir(parents=True, exist_ok=True)


def main_effects(y: str) -> dict[str, dict[str, float]]:
    """Exploratory: per-mouse mean residual (knockout minus same-mouse same-niche reference), t-test across mice."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from step36_perturbmap_niche import build, residual  # noqa: PLC0415

    les = build()
    out = {}
    for g in json.loads((R / "step36_perturbmap.json").read_text())["design"]["eligible"]:
        e = residual(les, g, y)
        k = les[les.gene == g].assign(e=e).dropna(subset=["e"])
        m = k.groupby("mouse").e.mean()
        t = stats.ttest_1samp(m, 0.0)
        half = stats.t.ppf(0.975, len(m) - 1) * m.std(ddof=1) / np.sqrt(len(m))
        out[g] = {"mean_of_mouse_means": float(m.mean()), "lo": float(m.mean() - half), "hi": float(m.mean() + half),
                  "p": float(t.pvalue), "mice": int(len(m))}
    return out


def main() -> None:
    d = json.loads((R / "step36_perturbmap.json").read_text())
    res = d["results"]["y_CD8"]
    genes = sorted(res, key=lambda g: res[g]["psi"])
    me = main_effects("y_CD8")
    (R / "step36b_main_effects.json").write_text(json.dumps({"note": "exploratory, not in AMENDMENT-010", "y_CD8": me}, indent=2))
    fig, axes = plt.subplots(1, 4, figsize=(S.DOUBLE, 62 * S.MM), gridspec_kw={"width_ratios": [1, 1.1, 1.1, 0.9]})
    # a: lesions per niche
    ax = axes[0]
    y = np.arange(len(genes))
    rich = np.array([res[g]["immune_rich"] for g in genes])
    tot = np.array([res[g]["lesions"] for g in genes])
    ax.barh(y, tot - rich, color=S.NICHE_COLOR[0], label="immune-poor")
    ax.barh(y, rich, left=tot - rich, color=S.NICHE_COLOR[1], label="immune-rich")
    ax.set_yticks(y, genes, style="italic")
    ax.set_xlabel("lesions (8 mice)")
    ax.set_title("Lesions per niche", loc="left")
    ax.legend(fontsize=5, loc="upper right")
    pd.DataFrame({"gene": genes, "lesions": tot, "immune_rich": rich}).to_csv(SD / "pm_a_lesions.csv", index=False)
    # b: Psi with CR2 (mouse-clustered) 95% CI
    ax = axes[1]
    psi = np.array([res[g]["psi"] for g in genes])
    lo = np.array([res[g]["ci_lo"] for g in genes])
    hi = np.array([res[g]["ci_hi"] for g in genes])
    ax.errorbar(psi, y, xerr=[psi - lo, hi - psi], fmt="o", c=S.OI["blue"], capsize=1.5, ms=3)
    for i, g in enumerate(genes):
        ax.text(max(hi.max(), 1) * 1.02, i, f"q={res[g]['q_calibrated']:.2f}", va="center", fontsize=4.5)
    ax.axvline(0, c="black", lw=0.4, ls=":")
    ax.set_yticks(y, [""] * len(genes))
    ax.set_xlabel("Ψ, CD8 (immune-rich − poor)")
    ax.set_title("Niche dependence (none at q < 0.1)", loc="left")
    pd.DataFrame({"gene": genes, "psi": psi, "ci_lo": lo, "ci_hi": hi, "p_perm": [res[g]["p_perm"] for g in genes],
                  "p_calibrated": [res[g]["p_calibrated"] for g in genes], "q": [res[g]["q_calibrated"] for g in genes]}).to_csv(SD / "pm_b_psi.csv", index=False)
    # c: main effect (sanity check)
    ax = axes[2]
    mm = np.array([me[g]["mean_of_mouse_means"] for g in genes])
    mlo = np.array([me[g]["lo"] for g in genes])
    mhi = np.array([me[g]["hi"] for g in genes])
    ax.errorbar(mm, y, xerr=[mm - mlo, mhi - mm], fmt="s", c=S.OI["vermillion"], capsize=1.5, ms=3)
    ax.axvline(0, c="black", lw=0.4, ls=":")
    ax.set_yticks(y, [""] * len(genes))
    ax.set_xlabel("Δ, CD8 vs other knockouts")
    ax.set_title("Main effect (exploratory)", loc="left")
    pd.DataFrame({"gene": genes, **{k: [me[g][k] for g in genes] for k in ("mean_of_mouse_means", "lo", "hi", "p", "mice")}}).to_csv(SD / "pm_c_main.csv", index=False)
    # d: permutation vs calibrated p, all outcomes
    ax = axes[3]
    rows = []
    for yk, mk in (("y_CD8", "o"), ("y_CD4", "s"), ("y_B220", "^")):
        for g, v in d["results"][yk].items():
            rows.append({"outcome": yk[2:], "gene": g, "p_perm": v["p_perm"], "p_calibrated": v["p_calibrated"]})
            ax.plot(v["p_perm"], v["p_calibrated"], mk, ms=2.5, c=S.OI["green"], mfc="white" if yk != "y_CD8" else S.OI["green"],
                    label=yk[2:] if g == next(iter(d["results"][yk])) else None)
    ax.legend(fontsize=5, loc="lower right")
    ax.plot([0.001, 1], [0.001, 1], c="black", lw=0.4, ls=":")
    ax.axhline(0.05, c=S.OI["grey"], lw=0.4)
    ax.axvline(0.05, c=S.OI["grey"], lw=0.4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("within-mouse permutation p")
    ax.set_ylabel("calibrated p")
    ax.set_title("Calibration", loc="left")
    pd.DataFrame(rows).to_csv(SD / "pm_d_calibration.csv", index=False)
    for i, ax in enumerate(axes):
        S.letter(ax, "abcd"[i], x=-0.45 if i == 0 else -0.15)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.18, top=0.88, wspace=0.35)
    S.save(fig, F / "fig_perturbmap")
    print(json.dumps(me, indent=1))


if __name__ == "__main__":
    main()
