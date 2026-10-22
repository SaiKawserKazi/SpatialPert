"""Shared figure style for the Nature Methods submission (Nature figure guide).

89 mm single / 183 mm double column; Arial 5-7 pt; panel letters 8 pt bold lowercase;
line weights 0.25-1 pt; Okabe-Ito colour-blind-safe palette; one fixed colour per entity.
"""

from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM

# Okabe-Ito
OI = {
    "black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999",
}

# One colour per test, used in every figure
METHOD_COLOR = {
    "cell": OI["vermillion"], "anova": OI["orange"], "sandwich": OI["yellow"], "tile": OI["purple"],
    "pseudobulk": OI["grey"], "patch30": OI["sky"], "stud_patch": "#7FB8D8", "lesion": OI["blue"],
    "stud_lesion": "#3B6EA5", "wcb": OI["purple"], "cr2": "#8C6BB1", "mixedlm": OI["grey"],
    "signflip": OI["green"], "gmn": OI["black"],
}
METHOD_NAME = {
    "cell": "Cell-level permutation", "anova": "Per-gene ANOVA", "sandwich": "Cluster-robust (CR1)",
    "tile": "Tile permutation", "pseudobulk": "Pseudobulk ANOVA", "patch30": "Patch permutation",
    "stud_patch": "Studentised patch", "lesion": "Lesion permutation", "stud_lesion": "Studentised lesion",
    "wcb": "Wild cluster bootstrap", "cr2": "CR2 (Bell-McCaffrey)", "mixedlm": "Mixed model",
    "signflip": "Within-lesion sign-flip", "gmn": "Geometry-matched null",
}
# Niches: a separate sequential-free palette never reused for methods
NICHE_COLOR = ["#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100"]
KO_COLOR = {"CIITA": "#E69F00", "DIO3": "#56B4E9", "ESAM": "#009E73", "PRSS23": "#CC79A7"}


def setup() -> None:
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
        "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.9, "lines.markersize": 3,
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 600, "legend.frameon": False,
    })


def letter(ax, s: str, x: float = -0.18, y: float = 1.06) -> None:  # noqa: ANN001
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def scale_bar(ax, length_um: float, label: str, *, loc: tuple[float, float] = (0.05, 0.05)) -> None:  # noqa: ANN001
    """Draw a horizontal scale bar in data units (µm) at an axes-fraction location."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs = x0 + loc[0] * (x1 - x0)
    ys = y0 + loc[1] * (y1 - y0)
    ax.plot([xs, xs + length_um], [ys, ys], color="black", lw=1.0, solid_capstyle="butt")
    ax.text(xs + length_um / 2, ys + 0.02 * (y1 - y0), label, ha="center", va="bottom", fontsize=5.5)


def save(fig, path_stem) -> None:  # noqa: ANN001
    fig.savefig(f"{path_stem}.pdf")
    fig.savefig(f"{path_stem}.png", dpi=600)
    plt.close(fig)
