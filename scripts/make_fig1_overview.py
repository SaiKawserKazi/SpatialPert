"""Fig. 1 as an illustrated overview (graphical-abstract style).

Top: a knockout grows as clonal lesions in two niches; counting cells gives a false niche
effect, counting lesions does not (real Xenium CIITA numbers). Bottom: calibrate against
matched control pseudo-knockouts, then report or flag. Stat cards: numbers read from reports.
Writes reports/figures/v4/fig1_overview.{pdf,png}.
"""

import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
OUT = R / "figures" / "v4"
NAVY, TEAL, CORAL = "#1F3A5F", "#1F9E89", "#E0644A"
BLUE_BG, PINK_BG, PAGE = "#E3EDF8", "#F9E4E0", "#F7F9FC"
CELL, CELL_EDGE, GREY = "#F2A98F", "#D9775C", "#B8C2CC"

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none"})


def j(name: str) -> dict:  # type: ignore[type-arg]
    return json.loads((R / f"{name}.json").read_text())


def card(ax, x, y, w, h, *, fc="white", ec="#D5DCE4", lw=1.2, r=1.6, z=2):  # noqa: ANN001, PLR0913
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z))


def tick(ax, x, y, r, colour, *, filled=True, mark=None):  # noqa: ANN001, PLR0913
    """Round badge with a drawn check mark (no font glyph needed)."""
    ax.add_patch(Circle((x, y), r, fc=colour if filled else "white", ec=colour, lw=1.5, zorder=6))
    c = mark or ("white" if filled else colour)
    ax.plot([x - 0.45 * r, x - 0.1 * r, x + 0.5 * r], [y, y - 0.38 * r, y + 0.42 * r], c=c, lw=2.4, solid_capstyle="round", zorder=7)


def cross(ax, x, y, r, colour, *, filled=True):  # noqa: ANN001
    ax.add_patch(Circle((x, y), r, fc=colour if filled else "white", ec=colour, lw=1.5, zorder=6))
    c = "white" if filled else colour
    d = 0.38 * r
    ax.plot([x - d, x + d], [y - d, y + d], c=c, lw=2.4, solid_capstyle="round", zorder=7)
    ax.plot([x - d, x + d], [y + d, y - d], c=c, lw=2.4, solid_capstyle="round", zorder=7)


def tissue(ax, x, y, w, h, rng, *, lesions, ring=True, scale=1.0, cell=CELL, edge=CELL_EDGE):  # noqa: ANN001, PLR0913
    """Two-niche tissue patch with grey control cells and clonal lesions."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=2", fc="white", ec="#C9D3DD", lw=1.0, zorder=2)
    ax.add_patch(box)
    for i, col in enumerate((BLUE_BG, PINK_BG)):
        rect = Rectangle((x + i * w / 2, y), w / 2, h, fc=col, ec="none", zorder=2.1)
        rect.set_clip_path(box)
        ax.add_patch(rect)
    pts = rng.uniform([x + 1, y + 1], [x + w - 1, y + h - 1], size=(int(60 * scale), 2))
    ax.scatter(pts[:, 0], pts[:, 1], s=4 * scale, c=GREY, lw=0, zorder=3)
    for cx, cy in lesions:
        if ring:
            ax.add_patch(Circle((cx, cy), 3.2 * scale, fc="none", ec=NAVY, lw=0.9, ls=(0, (2.5, 1.5)), zorder=4))
        for _ in range(8):
            px, py = np.array([cx, cy]) + rng.normal(0, 1.1 * scale, 2)
            ax.add_patch(Ellipse((px, py), 1.4 * scale, 1.0 * scale, angle=rng.uniform(0, 180), fc=cell, ec=edge, lw=0.5, zorder=5))


def main() -> None:  # noqa: PLR0915 - linear drawing
    rng = np.random.default_rng(7)
    ad, rl, pm = j("step19_adoption"), j("step28_real_smallcluster"), j("step36_perturbmap")
    cell_fp = ad["max_type1"]["loo_cell"]
    pow20 = next(r["power_size_adjusted"] for r in ad["per_test"]["loo_lesion"] if r["lesions"] == 20)
    ciita = rl["CIITA"]
    tg = pm["results"]["y_CD8"]["Tgfbr2"]

    fig = plt.figure(figsize=(7.2, 7.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 100, 100, fc=PAGE, ec="none", zorder=0))
    for (cx, cy, w, h, col) in ((22, 72, 46, 34, "#EAF0F8"), (80, 83, 40, 22, "#FBECE8"), (80, 58, 40, 22, "#E7F4F1"), (50, 32, 96, 24, "#EEF2F8")):
        ax.add_patch(Ellipse((cx, cy), w, h, fc=col, ec="none", zorder=0.5))

    ax.text(50, 96.5, "Does a knockout act differently in different tissue niches?", ha="center", va="center", fontsize=12.5, fontweight="bold", color=NAVY)

    # --- top: tissue with clonal lesions
    tissue(ax, 6, 56, 34, 30, rng, lesions=[(13, 77), (21, 63), (30, 79), (35, 64)])
    ax.text(14.5, 87.5, "niche A", ha="center", fontsize=7.5, color=NAVY)
    ax.text(31.5, 87.5, "niche B", ha="center", fontsize=7.5, color=NAVY)
    card(ax, 3, 89.6, 21, 4.4, fc="#4A5C82", ec="none", r=1.5)
    ax.text(13.5, 91.8, "spatial CRISPR", ha="center", va="center", fontsize=8.5, fontweight="bold", color="white", zorder=3)
    ax.text(23, 52.5, "knockout cells grow as clonal lesions", ha="center", fontsize=8, color="#5B6B7C")

    # arrows to the two ways of counting
    ax.add_patch(FancyArrowPatch((42, 78), (62, 86), connectionstyle="arc3,rad=-0.25", arrowstyle="-|>", mutation_scale=16, lw=3.2, color=CORAL, zorder=3))
    ax.add_patch(FancyArrowPatch((42, 64), (62, 58), connectionstyle="arc3,rad=0.25", arrowstyle="-|>", mutation_scale=16, lw=3.2, color=TEAL, zorder=3))
    ax.text(49, 89.3, "count cells", ha="center", fontsize=8.5, fontweight="bold", color=CORAL)
    ax.text(49, 54.3, "count lesions", ha="center", fontsize=8.5, fontweight="bold", color=TEAL)

    # result: cells as replicates (false positive)
    card(ax, 63, 78, 25, 15)
    for i in range(24):
        ax.add_patch(Circle((65.5 + (i % 8) * 2.6, 89.5 - (i // 8) * 2.6), 0.75, fc=CELL, ec=CELL_EDGE, lw=0.4, zorder=3))
    ax.text(75.5, 80.8, f"{ciita['cells']} cells → p = {ciita['cell_p']:.3f}", ha="center", fontsize=8, color=NAVY, zorder=3)
    cross(ax, 93, 85.5, 3.6, CORAL)
    ax.text(80, 73.6, "apparent niche effect", ha="center", fontsize=8.5, fontweight="bold", color=CORAL)

    # result: lesions as replicates (honest answer)
    card(ax, 63, 51, 25, 15)
    for i in range(min(ciita["lesions"], 10)):
        cx, cy = 66.5 + (i % 5) * 4.4, 62.3 - (i // 5) * 4.4
        ax.add_patch(Circle((cx, cy), 1.6, fc="none", ec=NAVY, lw=0.8, ls=(0, (2, 1.2)), zorder=3))
        ax.add_patch(Circle((cx, cy), 0.8, fc=CELL, ec=CELL_EDGE, lw=0.4, zorder=3))
    ax.text(75.5, 53, f"{ciita['lesions']} lesions → p = {ciita['lesion_p']:.2f}", ha="center", fontsize=8, color=NAVY, zorder=3)
    tick(ax, 93, 58.5, 3.6, TEAL)
    ax.text(80, 46.8, "the lesion is the replicate", ha="center", fontsize=8.5, fontweight="bold", color=TEAL)
    ax.text(75.5, 69.6, "same knockout (CIITA, Xenium)", ha="center", fontsize=7, style="italic", color="#5B6B7C")

    ax.plot([6, 94], [43.5, 43.5], c="#C9D3DD", lw=1.5)

    # --- bottom: calibrate before trusting
    ax.text(50, 40, "Before trusting a niche effect, calibrate it in the same tissue", ha="center", fontsize=11, fontweight="bold", color=NAVY)
    for k in range(4):
        x0 = 5 + k * 10.5
        les = [(x0 + rng.uniform(2, 7.5), 25 + rng.uniform(2, 8)) for _ in range(3)]
        tissue(ax, x0, 24.5, 9.5, 10.5, rng, lesions=les, ring=False, scale=0.38, cell="#9DB3D1", edge="#6F89AE")
    ax.text(25.5, 21.2, "control pseudo-knockouts with the same lesion sizes", ha="center", fontsize=7.5, color="#5B6B7C")
    ax.add_patch(FancyArrowPatch((47, 30), (53, 30), arrowstyle="-|>", mutation_scale=14, lw=3, color="#8DA2C8", zorder=3))
    # gauge funnel
    ax.add_patch(Polygon([[55, 36.5], [70, 36.5], [66.5, 23], [58.5, 23]], closed=True, fc="#E6ECF5", ec=NAVY, lw=1.8, zorder=3))
    for i, col in enumerate(("#E0644A", "#EFA25B", "#E9C85A", "#A9D27A", "#4BAF8A")):
        ax.add_patch(Rectangle((56.6 + i * 2.4, 32.6), 2.3, 2.2, fc=col, ec="none", zorder=4))
    ax.plot([62.5, 66], [27, 31.6], c="#1B1B1B", lw=2.2, solid_capstyle="round", zorder=5)
    ax.add_patch(Circle((62.5, 27), 0.9, fc="#1B1B1B", zorder=6))
    ax.text(62.5, 19.8, "rank the knockout's p\namong pseudo-knockouts", ha="center", fontsize=7.5, fontweight="bold", color=NAVY, linespacing=1.2)
    ax.add_patch(FancyArrowPatch((71.5, 32.5), (77, 34.5), arrowstyle="-|>", mutation_scale=14, lw=3, color=TEAL, zorder=3))
    ax.add_patch(FancyArrowPatch((71.5, 27.5), (77, 25.5), arrowstyle="-|>", mutation_scale=14, lw=3, color=CORAL, zorder=3))
    card(ax, 77.5, 31.2, 19, 6, fc=TEAL, ec="none", r=1.6)
    tick(ax, 80.6, 34.2, 1.8, "white", mark=TEAL)
    ax.text(89, 34.2, "report Ψ ± CI", ha="center", va="center", fontsize=8.5, fontweight="bold", color="white", zorder=7)
    card(ax, 77.5, 22.8, 19, 6, fc="white", ec=CORAL, lw=1.8, r=1.6)
    cross(ax, 80.6, 25.8, 1.8, CORAL, filled=False)
    ax.text(89.3, 25.8, "flag low power", ha="center", va="center", fontsize=8.5, fontweight="bold", color=CORAL, zorder=7)

    # --- stat cards
    stats = [
        (f"{cell_fp * 100:.0f}%", "of null draws called niche-\ndependent by cell-level tests", NAVY),
        ("20 lesions", f"per knockout: {pow20 * 100:.0f}% power for a\nmoderate niche effect", TEAL),
        (f"{tg['p_perm']:.3f} → {tg['p_calibrated']:.3f}", "Tgfbr2 p before vs after\ncalibration (Perturb-map)", CORAL),
    ]
    for i, (big, small, col) in enumerate(stats):
        x0 = 4 + i * 31.5
        card(ax, x0, 2.5, 29, 15, ec=col, lw=2.0, r=2)
        ax.text(x0 + 14.5, 13.2, big, ha="center", va="center", fontsize=15 if len(big) < 12 else 12.5, fontweight="bold", color=col, zorder=3)
        ax.plot([x0 + 9, x0 + 20], [10.2, 10.2], c=col, lw=2.6, zorder=3)
        ax.text(x0 + 14.5, 6.2, small, ha="center", va="center", fontsize=7.5, color="#2C3E50", linespacing=1.25, zorder=3)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig1_overview.png", dpi=300)
    fig.savefig(OUT / "fig1_overview.pdf")
    print("written", OUT / "fig1_overview.png")


if __name__ == "__main__":
    main()
