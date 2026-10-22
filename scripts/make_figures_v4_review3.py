"""Figure for review 3 (Priorities 1-2): why tests fail in brain, calibration, intervals, and
learned context. Writes reports/figures/v4/fig_calibration_context.{pdf,png} + source data.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nm_style as S  # noqa: E402

S.setup()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
F = R / "figures" / "v4"
SD = F / "source_data"
SD.mkdir(parents=True, exist_ok=True)
TISSUE = {"tumour": S.OI["blue"], "brain": S.OI["vermillion"]}


def wilson(k: float, n: int) -> tuple[float, float]:
    z = 1.959964
    p = k / n
    den = 1 + z**2 / n
    mid = (p + z**2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return mid - half, mid + half


def main() -> None:  # noqa: PLR0915 - linear figure assembly
    fig, axes = plt.subplots(1, 4, figsize=(S.DOUBLE, 52 * S.MM))
    # a: spread of lesion means by niche (heteroscedasticity) in both tissues
    ax = axes[0]
    rows = []
    for t in ("tumour", "brain"):
        d = json.loads((R / f"step29b_{t}_lesion_means.json").read_text())
        nic = sorted(d["per_niche"], key=int)
        sds = [d["per_niche"][k]["sd_lesion_means"] for k in nic]
        x = np.arange(len(nic)) + (0.18 if t == "brain" else -0.18)
        ax.bar(x, sds, width=0.34, color=TISSUE[t], label=f"{t} (var. ratio {d['between_lesion_var'] / (d['sampling_var_of_lesion_mean'] or np.nan):.1f})" if t == "brain" else t)
        rows += [{"tissue": t, "niche": f"N{k}", "sd_lesion_means": s, "n_lesions": d["per_niche"][k]["n_lesions"]} for k, s in zip(nic, sds, strict=True)]
    ax.set_xticks(range(4), [f"N{k}" for k in range(4)])
    ax.set_ylabel("SD of lesion-mean score")
    ax.set_title("Lesion baselines by niche", loc="left")
    ax.legend(fontsize=5, loc="upper left")
    pd.DataFrame(rows).to_csv(SD / "panel_a_lesion_spread.csv", index=False)
    # b: type-I error, raw vs calibrated (lesion permutation), both tissues
    ax = axes[1]
    cal = json.loads((R / "step31_calibrated_p.json").read_text())
    rows = []
    for t, key, test in (("tumour", "tumour_step19", "loo_lesion"), ("brain", "brain_step22", "lesion")):
        recs = cal[key][test]
        Ls = sorted(int(L) for L in recs)
        raw = [recs[str(L)]["raw_type1"] for L in Ls]
        cv = [recs[str(L)]["calibrated_type1"] for L in Ls]
        lo = [recs[str(L)]["calibrated_type1_ci"][0] for L in Ls]
        hi = [recs[str(L)]["calibrated_type1_ci"][1] for L in Ls]
        ax.plot(Ls, raw, ls="--", marker="o", mfc="white", c=TISSUE[t], label=f"{t}, permutation")
        ax.errorbar(Ls, cv, yerr=[np.subtract(cv, lo), np.subtract(hi, cv)], marker="o", c=TISSUE[t], capsize=1.5, label=f"{t}, calibrated")
        rows += [{"tissue": t, "lesions": L, "type1_permutation": r_, "type1_calibrated": c_, "ci_lo": a_, "ci_hi": b_}
                 for L, r_, c_, a_, b_ in zip(Ls, raw, cv, lo, hi, strict=True)]
    ax.axhline(0.05, c="black", lw=0.4, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([3, 5, 10, 20], [3, 5, 10, 20])
    ax.minorticks_off()
    ax.set_xlabel("lesions per knockout")
    ax.set_ylabel("type-I error")
    ax.set_title("Matched-null calibration", loc="left")
    ax.legend(fontsize=4.5, loc="upper left")
    pd.DataFrame(rows).to_csv(SD / "panel_b_calibration.csv", index=False)
    # c: CI coverage of Psi, CR2 vs naive
    ax = axes[2]
    rows = []
    for t in ("tumour", "brain"):
        d = json.loads((R / f"step33_{t}_coverage.json").read_text())
        Ls = [5, 10, 20]
        for kind, ls in (("cov_cr2", "-"), ("cov_naive", "--")):
            ys = [d[f"L{L}_lam0.5"][kind] for L in Ls]
            ax.plot(Ls, ys, ls=ls, marker="o", c=TISSUE[t], mfc="white" if kind == "cov_naive" else TISSUE[t],
                    label=f"{t}, {'CR2 lesion-clustered' if kind == 'cov_cr2' else 'naive cell-level'}")
            rows += [{"tissue": t, "interval": kind, "lesions": L, "coverage": y, "draws": d[f"L{L}_lam0.5"]["draws"]} for L, y in zip(Ls, ys, strict=True)]
    ax.axhline(0.95, c="black", lw=0.4, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([5, 10, 20], [5, 10, 20])
    ax.minorticks_off()
    ax.set_ylim(0.4, 1.0)
    ax.set_xlabel("lesions per knockout")
    ax.set_ylabel("95% CI coverage of Ψ")
    ax.set_title("Effect intervals", loc="left")
    ax.legend(fontsize=4.5, loc="center", bbox_to_anchor=(0.5, 0.6))
    pd.DataFrame(rows).to_csv(SD / "panel_c_coverage.csv", index=False)
    # d: representation benchmark pooled over three encoder seeds
    ax = axes[3]
    runs = [json.loads((R / f).read_text()) for f in ("step35_representation.json", "step35_representation_seed1.json", "step35_representation_seed2.json")]
    reps = [("kernel", "fixed kernel", S.OI["grey"]), ("learned", "learned context", S.OI["green"]), ("coord", "coordinates", S.OI["sky"]), ("expr", "own expression", S.OI["vermillion"])]
    rows = []
    for j, (key, name, col) in enumerate(reps):
        for i, L in enumerate((10, 20)):
            k = sum(r["summary"][f"L{L}_{key}"]["size_adjusted_power"] * r["summary"][f"L{L}_{key}"]["m_alt"] for r in runs)
            n = sum(r["summary"][f"L{L}_{key}"]["m_alt"] for r in runs)
            t1 = [r["summary"][f"L{L}_{key}"]["type1"] for r in runs]
            lo, hi = wilson(k, n)
            x = i * 5 + j
            ax.bar(x, k / n, color=col, width=0.8, label=name if i == 0 else None)
            ax.errorbar(x, k / n, yerr=[[k / n - lo], [hi - k / n]], c="black", capsize=1.2, lw=0.6)
            if max(t1) > 0.06:
                ax.text(x, hi + 0.008, "×", ha="center", fontsize=7, fontweight="bold", color=S.OI["vermillion"])
            rows.append({"representation": name, "lesions": L, "size_adjusted_power": k / n, "lo": lo, "hi": hi, "draws": n, "max_type1": max(t1)})
    ax.set_xticks([1.5, 6.5], ["10 lesions", "20 lesions"])
    ax.set_ylabel("size-adjusted power")
    ax.set_title("Niche representation", loc="left")
    ax.legend(fontsize=4.5, loc="upper left")
    pd.DataFrame(rows).to_csv(SD / "panel_d_representation.csv", index=False)
    for i, ax in enumerate(axes):
        S.letter(ax, "abcd"[i], x=-0.32)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.2, top=0.86, wspace=0.5)
    S.save(fig, F / "fig_calibration_context")
    print("written", F)


if __name__ == "__main__":
    main()
