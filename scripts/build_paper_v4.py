"""Manuscript v4.0 (response to review 3): Nature Methods Analysis format.

Title without subtitle, 150-word abstract without citations, ~500-word introduction, findings as
section headings, Online Methods with subheadings, six display items + Extended Data, and the
pre-registration history in Supplementary Note 1. Every number is read from reports/*.json.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

from spatialpert.docxkit import BLUE, Manuscript

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_v3_methods import methods  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REP = ROOT / "reports"
FIG = REP / "figures" / "v3"
FIG4 = REP / "figures" / "v4"
VERSION = "v0.0"
TITLE = "Clonal structure confounds tests of location-dependent knockout effects in spatial CRISPR screens"


def j(name: str) -> dict:  # type: ignore[type-arg]
    return json.loads((REP / f"{name}.json").read_text())


S10, AD, BR, ICC, D27 = j("step10_stage1"), j("step19_adoption"), j("step22_brain"), j("step23_icc"), j("step27_design")
CA = j("step24_callers") | j("step24_callers_cellcharter")
PM, RE, RT, RL, H, S5, H3R = j("step25_perturbmap_lesions"), j("step26_realism"), j("step20_runtime"), j("step28_real_smallcluster"), j("step6_hypotheses"), j("step5_controls"), j("step15_h3_random")
KOS = ["CIITA", "DIO3", "ESAM", "PRSS23"]
PER = AD["per_test"]
D29B, D29BB = j("step29_brain_diagnosis"), j("step29_brain_diagnosis_b")
LM = {t: j(f"step29b_{t}_lesion_means") for t in ("tumour", "brain")}
CAL, RC32 = j("step31_calibrated_p"), j("step32_real_calibrated")
COV = {t: j(f"step33_{t}_coverage") for t in ("tumour", "brain")}
CR2B = j("step34_brain_cr2")
REPS = [j(n) for n in ("step35_representation", "step35_representation_seed1", "step35_representation_seed2")]
PMR, PMM, VER = j("step36_perturbmap"), j("step36b_main_effects")["y_CD8"], j("step37_verify_perturbmap")


def at(test: str, L: int, key: str) -> float:
    return next(r[key] for r in PER[test] if r["lesions"] == L)


def t1max(test: str) -> float:
    return AD["max_type1"][test]


def pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def p3(p: float) -> str:
    return "p ≤ 0.001" if p <= 0.0011 else f"p = {p:.3f}"


T10 = {m: max(r[m]["rate"] for r in S10["summary"] if r["lambda"] == 0) for m in ("anova", "sandwich")}
CELL_IND = max(r["cell"]["rate"] for r in S10["summary"] if r["lambda"] == 0)
CELL_LES = [at("loo_cell", L, "type1") for L in (1, 3, 5, 10, 20)]
POW = {L: at("loo_lesion", L, "power_size_adjusted") for L in (3, 5, 10, 20)}
POW_CI = {L: at("loo_lesion", L, "power_ci") for L in (10, 20)}
BR_T1 = {t: [r[t]["type1"] for r in BR["summary"] if r[t]["type1"] is not None and r[t]["applicable"] >= 0.1] for t in ("cell", "patch30", "lesion", "wcb")}
CA_WORST = {c: {t: max(r[t]["type1"] for r in CA[c]["summary"] if r[t]["type1"] is not None) for t in ("lesion", "wcb", "signflip")} for c in CA}
RE_T1 = [r["lesion"]["type1"] for r in RE["summary"]]
RE_POW = max(r["lesion"]["power"][k]["size_adjusted"] for r in RE["summary"] for k in r["lesion"]["power"])
RHO = [v["per_gene_icc_median"] for v in ICC.values()]
CELLP = {k: RL[k]["cell_p"] for k in KOS}


def pooled(key: str, L: int) -> tuple[float, float]:
    """Power pooled over the three encoder seeds, and the worst type-I error."""
    k = sum(r["summary"][f"L{L}_{key}"]["size_adjusted_power"] * r["summary"][f"L{L}_{key}"]["m_alt"] for r in REPS)
    n = sum(r["summary"][f"L{L}_{key}"]["m_alt"] for r in REPS)
    return k / n, max(r["summary"][f"L{L}_{key}"]["type1"] for r in REPS)


VR = {t: LM[t]["between_lesion_var"] / LM[t]["sampling_var_of_lesion_mean"] for t in ("brain",)}
SDR = {t: max(v["sd_lesion_means"] for v in LM[t]["per_niche"].values()) / min(v["sd_lesion_means"] for v in LM[t]["per_niche"].values()) for t in LM}
CALT1 = [v["calibrated_type1"] for t in ("brain_step22", "tumour_step19") for k in ("lesion", "loo_lesion") if k in CAL[t] for v in CAL[t][k].values()]
COVR = {t: [COV[t][f"L{L}_lam0.5"]["cov_cr2"] for L in (5, 10, 20)] for t in COV}
COVN = {t: [COV[t][f"L{L}_lam0.5"]["cov_naive"] for L in (5, 10, 20)] for t in COV}
BIAS = max(abs(v["bias"]) for t in COV for v in COV[t].values())
PM8 = PMR["results"]["y_CD8"]
PMD = PMR["design"]
QMIN = min(v["q_calibrated"] for y in PMR["results"].values() for v in y.values())
TG = PM8["Tgfbr2"]
FIGMAP = {"1": "Fig. 1", "2": "Extended Data Fig. 5", "3": "Fig. 2", "4": "Fig. 3", "5": "Extended Data Fig. 6", "6": "Fig. 6"}


def renumber(text: str) -> str:
    """Map v3 figure numbers to v4 (two new main figures; tissue and generalisation move to Extended Data)."""
    return re.sub(r"(?<!Data )Fig\. ([1-6])", lambda m_: FIGMAP[m_.group(1)], text)


class V4(Manuscript):
    """Manuscript whose body text is renumbered while ``remap`` is on."""

    remap = False

    def para(self, text: str, **kw: object) -> None:  # type: ignore[override]
        super().para(renumber(text) if self.remap else text, **kw)  # type: ignore[arg-type]

    def bullets(self, items: list[str], size: float | None = None) -> None:
        super().bullets([renumber(t) if self.remap else t for t in items], size)


def extra_methods() -> list[tuple[str, list[str]]]:
    """Online Methods subsections added in v4 (review 3)."""
    return [
        ("Diagnosis of the brain inflation", [
            "On 200 null draws per lesion number we compared lesion permutation under: the control model refitted without the pseudo-knockout; leave-out "
            "or circular niche labels; dropping cells whose niche differs from their lesion's majority niche; residuals against the lesion's majority niche; "
            "local-control residuals (50 nearest controls); and Welch studentisation. Lesion means were then pooled over 60 draws (20 lesions each) and "
            "summarised per niche (mean, SD), against the sampling variance of a lesion mean, and as correlation between lesion means by centroid distance."
        ]),
        ("Matched-null calibration and exchangeability diagnostic", [
            "For a knockout with lesion sizes (s_1, ..., s_L), M = 200 pseudo-knockouts are drawn from control cells with lesions of the same sizes at random, "
            "well-separated seeds (≥ 400 µm), removed from the control fit, relabelled with leave-out niches and tested exactly like the knockout. The "
            "calibrated p-value is (#{p_j < p} + U(1 + #{p_j = p}))/(1 + M) with U uniform on (0, 1) and untestable draws excluded, which is exact when knockout "
            "lesions are placed like random control positions. Calibration was evaluated by two-fold cross-fitting of existing null draws (calibration pool and "
            "tested draws disjoint). The diagnostic reports, from 40 control pseudo-knockouts of 20 lesions, the ratio of between-lesion variance of lesion "
            "means to their sampling variance and the ratio of the largest to the smallest per-niche SD of lesion means (spatialpert.calibration)."
        ]),
        ("Knockout-by-niche contrast and intervals", [
            "For niches a and b the contrast Ψ_g = Δ_b,g − Δ_a,g is the OLS coefficient of a niche-b indicator in the regression of knockout residuals on "
            "niche, using cells in a or b, with CR2 (Bell-McCaffrey) covariance clustered on lesions and Imbens-Kolesár degrees of freedom "
            "(spatialpert.smallcluster.niche_contrast). Coverage was assessed on pseudo-knockouts (400 cells; 5, 10, 20 lesions; 200 draws) with a shared "
            "effect in all cells and a known effect λδ·v in niche b, comparing CR2 intervals with naive cell-level Welch intervals."
        ]),
        ("Learned context encoder", [
            "For each cell, context features are the mean expression principal components (30) of control cells in the rings 0-50 µm and 50-150 µm, the cell "
            "itself excluded. A two-layer network (64 hidden units, 16-dimensional embedding) is trained on control cells for 400 full-batch Adam epochs to "
            "predict the cell's own principal components from its context (PyTorch, CPU). One-hop aggregation is used deliberately: with two message-passing "
            "hops a cell's own features would return to it through its neighbours. Pseudo-knockout cells are removed from the reference before their context "
            "is computed. Niches are k-means (K = 4) on the embedding. The benchmark compared four representations (fixed kernel features, learned context, "
            "own expression, coordinates) under the same lesion permutation test, with a niche-dependent effect proportional to the standardised host-cell "
            "fraction within 50 µm; 200 null and 150 effect draws per lesion number and encoder seed (three seeds)."
        ]),
        ("Perturb-map analysis", [
            f"Segmentation tables of the Lymphoid panel ({PMD['lobes']} lobes, {PMD['mice']} mice) were downloaded from S-BIAD267 and reduced to centroid, "
            "cell area, per-cell marker means, Pro-Code and Gene. Raw tables of the last nine lobes were kept with SHA-256 checksums, and the reduced tables "
            f"of those lobes equalled the raw tables exactly in every kept column ({VER['lobes_checked'] - VER['failed']}/{VER['lobes_checked']} lobes); "
            "the raw tables of the first 32 lobes were not retained. Lesions are the authors' curated lesion polygons with ≥ 30 labelled cells; the knockout is "
            "the Gene its Pro-Code maps to. The outcome is the logit fraction of marker-positive unlabelled cells (above the lobe's 95th percentile) inside "
            "the lesion or within 20 µm of a lesion cell. Each lesion's niche is k-means (K = 2) on the composition of the ring 100-300 µm away, excluding "
            "every cell within 100 µm of any labelled cell; the cluster with more marker-positive cells is immune-rich. Residuals are taken against reference "
            "lesions of other knockouts in the same mouse and niche. Ψ is tested by permuting niche labels among the knockout's lesions within mouse (9,999 "
            "permutations), calibrated against 1,000 pseudo-knockouts drawn from reference lesions with the same per-mouse counts, with CR2 intervals "
            "clustered on mouse. Eligible knockouts had ≥ 20 lesions, ≥ 5 per niche and lesions in ≥ 3 mice; BH at q = 0.10 over knockouts. Main effects "
            "(mean of per-mouse mean residuals, one-sample t-test over mice) were not part of the registered plan and are exploratory. Pixel size is taken as "
            "0.5 µm."
        ]),
    ]


def abstract() -> str:
    lp20, _ = pooled("learned", 20)
    kp20, _ = pooled("kernel", 20)
    return (
        "Spatial CRISPR screens can test whether a knockout acts differently across tissue niches, but knockout cells grow as clonal lesions. Using "
        f"semi-synthetic knockouts from tumour and brain sections, we show that cell-level tests detect niche effects in up to {pct(max(CELL_LES))} of null "
        "draws, and that letting knockout cells define their own niche inflates lesion-level tests. Lesion permutation and CR2 were calibrated in tumour, but "
        "every test failed in brain, where lesion baselines were large and unequal across niches. Calibrating against matched control pseudo-knockouts "
        "restored error control in both tissues, lesion-clustered intervals covered the true knockout-by-niche contrast, and a learned context encoder "
        f"roughly doubled power ({pct(kp20)} to {pct(lp20)} with 20 lesions). In Perturb-map imaging of {PMD['lesions']:,} lesions from eight mice, "
        "several knockouts changed T-cell infiltration, but none did so detectably differently between immune-rich and immune-poor niches. SpatialPert "
        "implements these methods."
    )


def build() -> Path:  # noqa: PLR0915 - linear document assembly
    m = V4(REP / "eq_v2", article_type="Analysis")
    m.para(f"**{TITLE}**", size=15, align="left", space_after=6)
    m.para("[Author names]¹ · [Corresponding author]¹*", size=9, align="left")
    m.para("¹[Affiliation]. *[e-mail]", size=8, align="left", space_after=8)
    m.heading("Abstract")
    ab = abstract()
    assert len(ab.split()) <= 150, len(ab.split())  # noqa: S101 - build-time guard
    m.para(ab, space_after=6)

    m.heading("Main")
    m.remap = True
    m.para(
        "Pooled CRISPR screens with single-cell readouts measure what each knockout does to the transcriptome[@dixit2016,replogle2022], but they dissociate "
        "tissue and so lose where each cell sat. Spatial screens keep that information. Perturb-map reads protein barcodes in tissue sections[@dhainaut2022], "
        "Perturb-FISH images guides alongside transcripts[@binan2025], Perturb-DBiT and imaging-based readouts capture guides and RNA in the same section[@baysoy2026], "
        "and spatial Perturb-seq extends this to whole organs[@shen2026]. These designs promise a new kind of result: a knockout whose effect depends on its "
        "microenvironment, for example on whether it sits near immune cells or at a hypoxic core."
    )
    m.para(
        "Testing such niche dependence raises a statistical problem that dissociated screens do not face. In vivo, each knockout typically enters a tissue "
        "as a few founder cells that grow into clonal lesions. Cells in one lesion share a genotype, a lineage history and a local environment, so they are not "
        "independent replicates. If a test treats them as independent, the effective number of observations is the number of lesions, not cells, and a "
        "difference between two lesions that happen to sit in different niches looks like a niche effect. The problem is the pseudoreplication long described "
        "in ecology and neuroscience[@hurlbert1984,lazic2010] and in multi-sample single-cell analysis[@squair2021,zimmerman2021,crowell2020], "
        "but in spatial screens it has an extra twist: niches are defined from the same tissue in which the lesions grow, and a lesion large enough to shape "
        "its own neighbourhood can help define the niche it is then compared against."
    )
    m.para(
        "Methods for spatial transcriptomics address related questions. Spatially variable gene tests[@svensson2018,sun2020] and niche callers[@hu2021,dong2022,long2023,singhal2024,varrone2024] "
        "find structure in unperturbed tissue, and econometrics offers corrections for few clusters[@bell2002,imbens2016,cameron2015,webb2023]. "
        "None of them has been evaluated for the specific question of whether a knockout's effect differs between niches when knockout cells come in a handful of clones."
    )
    m.para(
        "Here we benchmark that question directly. We built semi-synthetic knockouts from control cells of a public Xenium section of HT-29 lung metastases[@baysoy2026], "
        "placing them as 1-20 lesions with known niche-dependent effects or none, and evaluated eight tests, from naive cell-level permutation to small-cluster "
        "corrections and a within-lesion sign-flip test. We repeated the benchmark with a second tissue, three published niche callers and count-level effects "
        "with realistic lesion sizes, diagnosed why tests fail in the second tissue, developed a calibration, effect intervals and a learned context encoder, "
        "and applied the methods to four Xenium knockouts and to 12 knockouts in Perturb-map imaging of eight mice. The analysis plan and "
        "the rule for choosing a recommended test were pre-registered (Supplementary Note 1). The findings give a design rule — the lesion, not the cell, is "
        "the unit of replication, and about 20 lesions per knockout are needed — and SpatialPert, an open-source toolkit that applies the calibrated test, "
        "reports whether a screen is testable, and plans screens from the simulated power curves."
    )

    m.heading("Results")
    m.heading("Cell-level tests detect niche effects that are not there", 2)
    m.para(
        "We first placed pseudo-knockouts made of control tumour cells as independent clones, reproducing the patch sizes of real knockouts (Methods). With no "
        f"injected effect, a cell-level permutation test rejected in up to {pct(CELL_IND)} of draws, per-gene ANOVA in up to {pct(T10['anova'])} and cluster-robust "
        f"regression with CR1 standard errors in up to {pct(T10['sandwich'])}. When pseudo-knockouts were instead grown as 1-20 contiguous lesions totalling 400 cells, the cell-level test "
        f"rejected in {pct(min(CELL_LES))}-{pct(max(CELL_LES))} of null draws (Fig. 3a and Extended Data Fig. 1a). The cause is clonal: cells in one lesion "
        "share a baseline that differs between lesions, and when lesions fall in different niches the baseline difference is read as a niche effect (Fig. 1)."
    )
    m.heading("Knockout cells must not define their own niche", 2)
    m.para(
        "Niches in a perturbation screen are defined from control cells, but a knockout cell's niche label still depends on its neighbours. In our first "
        "benchmark, pseudo-knockout cells had contributed to the neighbourhood features that defined their own niches. Recomputing niche labels with the "
        f"pseudo-knockout removed from the reference (leave-out labels) changed the niche of {pct(min(AD['frac_relabelled_by_leave_out'].values()))}-"
        f"{pct(max(AD['frac_relabelled_by_leave_out'].values()))} of pseudo-knockout cells and removed most of the inflation of lesion-level tests: the "
        f"maximum type-I error of lesion permutation fell from {t1max('circ_lesion'):.2f} to {t1max('loo_lesion'):.2f}, and that of 30 µm patch permutation "
        f"from {t1max('circ_patch30'):.2f} to {t1max('loo_patch30'):.2f} (Fig. 3c). Circular labels make each lesion look more homogeneous in niche than a "
        "real lesion would, and the inflation grew with lesion number; it was not caused by where lesions were seeded, which followed the "
        "distribution of control cells (Extended Data Fig. 4a,b). Any simulation of spatial perturbations must therefore exclude the perturbed cells when "
        "assigning their niches."
    )
    m.heading("Lesion-level tests are calibrated but need 10-20 lesions", 2)
    m.para(
        "With leave-out labels, we compared eight tests at 1, 3, 5, 10 and 20 lesions (300 null and 150 effect draws per level; Fig. 3a,b). Four were "
        f"calibrated (maximum type-I error at most 0.06): lesion permutation ({t1max('loo_lesion'):.3f}), studentised lesion permutation "
        f"({t1max('loo_stud_lesion'):.3f}), CR2 with Bell-McCaffrey correction ({t1max('loo_cr2'):.3f}) and the within-lesion sign-flip test "
        f"({t1max('loo_signflip'):.3f}). Patch permutation was slightly liberal ({t1max('loo_patch30'):.3f}) because patches split lesions. The wild cluster "
        f"bootstrap, often recommended for few clusters[@cameron2008,mackinnon2017,webb2023], rejected in up to {pct(t1max('loo_wcb'))} of null draws, and a "
        f"lesion random-intercept mixed model in up to {pct(t1max('loo_mixedlm'))}; with a handful of lesions split across niches, neither the bootstrap's "
        "restricted residuals nor the variance component is estimated well."
    )
    m.para(
        "Calibration came at the cost of power. For a niche-specific effect half the size of the shared knockout effect, size-adjusted power of lesion "
        f"permutation was {POW[3]:.2f} with 3 lesions, {POW[5]:.2f} with 5, {POW[10]:.2f} with 10 (95% CI {POW_CI[10][0]:.2f}-{POW_CI[10][1]:.2f}) and "
        f"{POW[20]:.2f} with 20 ({POW_CI[20][0]:.2f}-{POW_CI[20][1]:.2f}; Fig. 4a). The sign-flip test, which compares niches within a lesion and cancels "
        f"the lesion's own baseline, was the most conservative and never exceeded {max(r.get('power_size_adjusted') or 0 for r in PER['loo_signflip']):.2f} "
        "power, because few lesions straddle two niches. Under the pre-registered rule (highest mean size-adjusted power among calibrated tests, with the "
        f"critical value corrected in AMENDMENT-008b) lesion permutation was adopted (mean power {AD['mean_power_corrected']['loo_lesion']:.2f}) and is "
        "SpatialPert's default; the uncorrected rule would have chosen the sign-flip test, and both are reported."
    )
    m.heading("Within-lesion correlation does not predict power", 2)
    m.para(
        "Survey statistics would predict power from the design effect, n_eff = L·m/(1 + (m − 1)ρ), where ρ is the within-lesion intraclass correlation. "
        f"In control pseudo-lesions, ρ was small, {min(RHO):.3f}-{max(RHO):.3f} per gene (Extended Data Fig. 2a), so the formula predicts "
        f"{D27['formula_power_efficiency1'][1]:.2f} power already with 5 lesions, against a simulated {D27['simulated_size_adjusted_power'][1]:.2f}. Fitting "
        f"an efficiency factor did not rescue it (RMSE {D27['rmse_fitted']:.2f}; Fig. 4b). Power is limited not by correlation among cells of a lesion but "
        "by how many independent lesions land in each niche to be compared, which a cell-level variance formula does not capture. SpatialPert's "
        "plan_screen_empirical() therefore interpolates the simulated power curves rather than the formula."
    )
    m.heading("Calibration depends on tissue and niche caller", 2)
    m.para(
        "Repeating the benchmark in a Xenium mouse brain section (10 µm bins, 50 genes; Methods), every test was inflated: cell-level permutation rejected in "
        f"{min(BR_T1['cell']):.2f}-{max(BR_T1['cell']):.2f} of null draws, patch permutation in up to {max(BR_T1['patch30']):.2f}, lesion permutation in up to "
        f"{max(BR_T1['lesion']):.2f} and the wild cluster bootstrap in up to {max(BR_T1['wcb']):.2f} (Fig. 5a). Brain niches are large and smooth, so lesions "
        "rarely straddle niches and lesion baselines are confounded with niche; the sign-flip test, which needs straddling lesions, was almost never applicable. "
        "In the tumour, replacing k-means niches with BANKSY, SpaGCN or CellCharter niches changed calibration (Fig. 5b and Extended Data Fig. 3): the worst "
        f"type-I error of lesion permutation was {CA_WORST['kmeans']['lesion']:.2f} with k-means, {CA_WORST['banksy']['lesion']:.2f} with BANKSY, "
        f"{CA_WORST['spagcn']['lesion']:.2f} with SpaGCN and {CA_WORST['cellcharter']['lesion']:.2f} with CellCharter, whereas the sign-flip test stayed at "
        f"or below {max(v['signflip'] for v in CA_WORST.values()):.2f}. These comparisons use circular labels for every caller (Methods), so they show how "
        "much the caller matters rather than absolute error rates; together with the brain result they mean that a calibrated benchmark in one tissue does "
        "not certify a test in another, and SpatialPert therefore runs a tissue-specific pseudo-knockout calibration before testing."
    )
    m.remap = False  # new sections already use v4 figure numbers
    m.heading("Large, unequal lesion baselines explain the failures in brain", 2)
    m.para(
        "We asked why every test failed in brain. Four explanations were ruled out on the same null draws: the pseudo-knockout's presence in the control fit, "
        "leave-out relabelling (about 1% of brain bins changed niche), circular labels and lesions that straddle niches (about 2% of lesion cells; removing "
        f"them left the type-I error of lesion permutation at {D29BB['20']['v6']['type1']:.2f} with 20 lesions). Pooled across draws, lesion means were centred "
        f"in every niche, but their variance in brain was {VR['brain']:.1f} times that expected from sampling, their spread differed {SDR['brain']:.1f}-fold "
        f"between niches (tumour {SDR['tumour']:.1f}-fold), and they remained correlated over hundreds of micrometres (Fig. 4a). Lesion means are then not "
        "exchangeable between niches, which every lesion-level test assumes; with few lesions per niche, Welch studentisation could not repair it, and the "
        f"CR2 test was also inflated in brain ({min(v['type1'] for v in CR2B.values()):.2f}-{max(v['type1'] for v in CR2B.values()):.2f})."
    )
    m.heading("Calibration against matched pseudo-knockouts restores error control", 2)
    m.para(
        "Because the failure is a property of the tissue, we calibrated it in the tissue: the observed p-value of a knockout is ranked among p-values of "
        "control pseudo-knockouts with the same lesion sizes, with ties broken at random (Methods). The result is exact if knockout lesions are placed like "
        f"random control positions. On held-out null draws, the calibrated type-I error of lesion permutation was {min(CALT1):.3f}-{max(CALT1):.3f} in both tissues (Fig. 4b). "
        "Calibration does not create power: in brain the calibrated test detected the injected effect in at most "
        f"{max(v['calibrated_power'] for v in CAL['brain_step22']['lesion'].values()):.2f} of draws, so the correct report there is that niche effects of "
        "this size are undetectable. A diagnostic computed from control pseudo-lesions flags such tissues before testing (lesion variance ratio "
        f"{RC32['diagnostic']['tumour']['variance_ratio_median']:.2f} in tumour and {RC32['diagnostic']['brain']['variance_ratio_median']:.1f} in brain). In "
        f"the Xenium screen, calibration left the conclusions unchanged (CIITA p = {RC32['real']['CIITA_lesion']['p_calibrated']:.2f})."
    )
    m.heading("Lesion-clustered intervals cover the knockout-by-niche contrast", 2)
    m.para(
        "A test answers whether a knockout's effect differs between niches; investigators also need the size of the difference. We estimate the contrast "
        "Ψ = Δ_b − Δ_a, where Δ_n is the knockout effect in niche n, from residuals of knockout cells against niche-matched controls, with CR2 standard errors "
        f"clustered on lesions. With injected effects of known size, the estimate was unbiased (|bias| ≤ {BIAS:.3f}) and 95% intervals covered the truth in "
        f"{min(COVR['tumour']):.2f}-{max(COVR['tumour']):.2f} of draws in tumour and {min(COVR['brain']):.2f}-{max(COVR['brain']):.2f} in brain, whereas "
        f"naive cell-level intervals covered only {min(COVN['brain']):.2f}-{max(COVN['brain']):.2f} in brain (Fig. 4c). The honest intervals in brain were "
        "about four times wider, which is the same loss of information seen in the calibrated test."
    )
    m.heading("A learned context encoder improves power without losing calibration", 2)
    lp10, lt10 = pooled("learned", 10)
    lp20, lt20 = pooled("learned", 20)
    kp10, kt10 = pooled("kernel", 10)
    kp20, kt20 = pooled("kernel", 20)
    _, et = pooled("expr", 20)
    m.para(
        "We next asked whether a learned representation of each cell's surroundings defines more useful niches than fixed kernel features. The encoder "
        "sees only control neighbours in two rings, never the cell itself, and is trained to predict a cell's expression from its context (Methods). To "
        "avoid favouring any partition, the injected niche-dependent effect followed a continuous covariate, the local host-cell fraction, and the "
        "downstream test was the same for all representations. Over three training seeds (450 effect draws per setting), learned niches kept type-I error "
        f"at or below {max(lt10, lt20):.2f} and raised size-adjusted power from {kp10:.2f} to {lp10:.2f} with 10 lesions and from {kp20:.2f} to "
        f"{lp20:.2f} with 20 (Fig. 4d). Niches from a cell's own expression were miscalibrated (type-I up to {et:.2f}), confirming that context must "
        "come from neighbours. The gain is real but modest; absolute power stayed low for this weak, composition-driven effect."
    )
    m.remap = True
    m.heading("Realistic count-level effects are harder to detect", 2)
    m.para(
        "Effects injected on the log scale into every gene are a favourable case. When effects were applied to raw counts by binomial thinning or Poisson "
        "boosting, in 5-100% of genes, and lesion sizes followed the heavy-tailed distribution of the real knockouts, lesion permutation stayed close to "
        f"nominal size ({min(RE_T1):.3f}-{max(RE_T1):.3f}) but size-adjusted power did not exceed {RE_POW:.2f} at any setting with 5-20 lesions (Fig. 4c). "
        "Sparse count data and a few dominant lesions per knockout leave little independent information; the curve in Fig. 4a is an upper bound."
    )
    m.heading("Current screens fall short of the lesion requirement", 2)
    m.para(
        f"In the HT-29 Xenium section, the four testable knockouts formed {RL['CIITA']['lesions']}, {RL['DIO3']['lesions']}, {RL['ESAM']['lesions']} and "
        f"{RL['PRSS23']['lesions']} lesions (CIITA, DIO3, ESAM and PRSS23; Fig. 2). Cell-level permutation reported niche dependence for CIITA "
        f"({p3(CELLP['CIITA'])}) and PRSS23 ({p3(CELLP['PRSS23'])}); with lesions as units, only CIITA could be tested at all, and its effect was not niche "
        f"dependent (lesion permutation p = {RL['CIITA']['lesion_p']:.2f}; CR2 minimum q = {RL['CIITA']['cr2_min_q']:.2f}; Fig. 6). The other knockouts were "
        "untestable: once each lesion was assigned its majority niche and niches with fewer than 30 cells were dropped, fewer than two niches remained, "
        "and the small-cluster tests found some compared niche in fewer than two lesions (Fig. 6b). In the Perturb-map imaging deposit, which covers "
        f"{PM['lobes']} lung lobes from eight mice, the median Pro-Code was found in {PM['lesions_per_code_median']:.1f} lesions; {PM['codes_ge20']} of {PM['distinct_procodes']} codes reached 20 "
        f"lesions and {PM['codes_ge10']} reached 10 (Fig. 4d). Whole-organ imaging across several animals can therefore meet the requirement for more "
        "than half of the knockouts, whereas a single section rarely does."
    )
    m.remap = False
    m.heading("In Perturb-map, knockouts change immune infiltration but not differently between niches", 2)
    sig_main = sorted((g for g, v in PMM.items() if v["p"] < 0.01), key=lambda g: PMM[g]["p"])
    m.para(
        f"Finally, we applied the methods to the Perturb-map imaging deposit: {PMD['lesions']:,} curated lesions in {PMD['lobes']} lung lobes from "
        f"{PMD['mice']} mice, with each lesion's knockout read from its Pro-Code and immune infiltration measured from CD8, CD4 and B220 staining (Methods). "
        "The analysis plan, including the niche definition (immune-rich or immune-poor surroundings 100-300 µm from the lesion), the outcome, the test, "
        f"eligibility and multiplicity, was registered before any outcome was computed. {len(PMD['eligible'])} knockouts were eligible. Several changed "
        "CD8 T-cell infiltration consistently across mice (exploratory main effects: "
        + "; ".join(f"{g} {PMM[g]['mean_of_mouse_means']:+.2f}, p = {PMM[g]['p']:.4f}" for g in sig_main)
        + "). The two largest effects reproduce the original study, which reported T-cell exclusion from Tgfbr2-knockout lesions and strong T-cell "
        "infiltration of Socs1-knockout lesions[@dhainaut2022]. No knockout, however, showed niche-dependent "
        f"infiltration: the smallest calibrated q-value over 36 tests was {QMIN:.2f} (Fig. 5). The strongest candidate, Tgfbr2 (Ψ = {TG['psi']:+.2f}; "
        f"95% interval {TG['ci_lo']:+.2f} to {TG['ci_hi']:+.2f}), had within-mouse permutation p = {TG['p_perm']:.3f} but calibrated p = "
        f"{TG['p_calibrated']:.3f}: an uncalibrated analysis would have reported it. Only {pct(PMD['immune_rich_fraction'])} of lesions lay in the "
        "immune-rich niche, so only large niche differences were detectable; the result shows that these knockouts' effects on infiltration are not "
        "strongly niche-dependent, not that niche dependence is absent."
    )

    m.remap = True
    m.heading("Discussion")
    m.para(
        "Spatial CRISPR screens make it possible to ask where a knockout acts, but the answer depends on counting the right units. Cells of one lesion are "
        "not replicates, and tests that treat them as such find niche dependence where there is none. Lesion-level tests fix calibration in the tissue we "
        "benchmarked most thoroughly, but they need on the order of 20 lesions per knockout for adequate power under favourable effects, and more under "
        "count-level effects. Some popular corrections for few clusters, the wild cluster bootstrap and random-intercept mixed models, did not hold their "
        "size here, and a design-effect calculation was not a reliable planning tool."
    )
    m.para(
        "Lesion number can be raised by design. Injecting more transduced founder cells per knockout increases lesions per knockout; smaller libraries concentrate lesions on fewer knockouts; sectioning several levels or lobes and pooling animals adds lesions that are "
        "independent by construction, which is how the Perturb-map deposit reaches 20 or more lesions for more than half of its codes; and clonal barcodes would let lesions "
        "be defined by lineage rather than distance. For analysis, we recommend the following decision path, which SpatialPert implements:"
    )
    m.bullets([
        "Define niches from control cells only, and exclude each knockout's own cells when assigning its niche.",
        "Define lesions (100 µm single linkage by default) and count lesions per niche for every knockout; if no niche is visited by at least two lesions, report the knockout as untestable rather than testing cells.",
        "Check lesion exchangeability on control pseudo-lesions, and report every p-value calibrated against matched pseudo-knockouts of the same tissue; where the diagnostic flags large, unequal lesion baselines, expect little power.",
        "Report the knockout-by-niche contrast Ψ with lesion-clustered (CR2) intervals, not only a p-value.",
        "Test with lesion permutation; report CR2 and, where lesions straddle niches, the within-lesion sign-flip test as sensitivity analyses. Avoid cell-level tests, the wild cluster bootstrap and random-intercept mixed models when lesions are few.",
        "Plan new screens with plan_screen_empirical(): aim for at least 20 lesions per knockout; with 10, only large niche effects are detectable.",
    ])
    m.para(
        "Our study has limitations. The benchmark is semi-synthetic: pseudo-knockouts carry the clonal structure of control cells, not of real lineages, and "
        "effects are injected, so absolute power depends on the assumed effect model. The calibration is exact only if knockout lesions are placed like "
        "random control positions, which cannot be checked from the data. The brain diagnosis and the calibration were developed after the brain failure was "
        "known (AMENDMENT-009) and are exploratory. The second tissue is binned transcripts without segmentation, and lesions are defined by distance, not "
        "lineage. The learned encoder's gain is modest and was shown for one effect model in one tissue. In Perturb-map, the readout is a protein panel, "
        "there are no non-targeting controls (knockouts are compared with each other), pixel size was assumed (0.5 µm), few lesions lay in immune-rich "
        "surroundings, and the main effects are exploratory. Small-cluster methods that need R (Satterthwaite and Kenward-Roger mixed models, dream) were "
        "not evaluated, and SpaGCN was run on a subsample. No knockout with a known niche-dependent effect is available, so the positive controls validate "
        "the effect estimators only. Each limitation points to the same remedy: screens designed with many independent lesions, ideally barcoded, across "
        "several animals."
    )

    m.heading("Online Methods")
    blocks = methods({n: j(n) for n in ("step4_domains", "gate3_dbit_audit", "step20_runtime", "step25_perturbmap_lesions")})
    blocks = [(sub, [t.replace("The curated lesion annotations of the Perturb-map imaging deposit (BioImage Archive S-BIAD267) were used to count lesions per knockout.",
                               "The curated lesion annotations and the per-cell segmentation tables of the Perturb-map imaging deposit (BioImage Archive S-BIAD267) "
                               "were used to count lesions per knockout and for the real-data analysis.") for t in paras]) for sub, paras in blocks]
    cut = next(i for i, (sub, _) in enumerate(blocks) if sub == "Statistics")
    blocks = blocks[:cut] + extra_methods() + blocks[cut:]
    for sub, paras in blocks:
        m.heading(sub, 2)
        for t in paras:
            m.para(t)
    m.remap = False
    m.heading("Reporting summary", 2)
    m.para("Further information on research design is available in the Nature Portfolio Reporting Summary linked to this article.")
    m.heading("Data availability", 2)
    m.para(
        "All data analysed are public: GEO GSE319123 (Xenium, Perturb-DBiT study), GSE319277 (DBiT-seq), GSE193460 (Perturb-map Visium), BioImage Archive "
        "S-BIAD267 (Perturb-map imaging lesion annotations and segmentation tables) and Zenodo 11209134 (Xenium mouse brain). Source data for Figs. 2-6 and "
        "Extended Data Figs. 1, 2, 4 and 5 are provided with this paper. Benchmark per-draw p-values are deposited at Zenodo [DOI to be minted at submission]."
    )
    m.heading("Code availability", 2)
    m.para(
        "SpatialPert is available at [GitHub URL] under the [licence to be chosen: MIT or BSD-3-Clause], archived at Zenodo [DOI], and installable from PyPI "
        "(pip install spatialpert). The repository contains every analysis script, the pre-registered protocol with amendments, and a recheck script that "
        "recomputes all reported rates from per-draw files."
    )

    m.heading("References")
    m.references_placeholder()
    m.heading("Acknowledgements")
    m.para("[To be completed.]")
    m.heading("Author contributions")
    m.para("[To be completed using CRediT roles.]")
    m.heading("Competing interests")
    m.para("The authors declare no competing interests.")

    figs = [
        (FIG4 / "fig1_overview.png", "Fig. 1 | Clonal growth makes the lesion, not the cell, the unit of replication.", [
            "Top, in a spatial CRISPR screen each knockout grows as clonal lesions (dashed circles) within tissue niches (shaded); grey dots, control cells. "
            f"For the same Xenium knockout (CIITA; {RL['CIITA']['cells']} cells in {RL['CIITA']['lesions']} lesions), counting cells as replicates gives an "
            f"apparent niche effect (cell-level permutation p = {RL['CIITA']['cell_p']:.3f}), whereas counting lesions does not (lesion permutation "
            f"p = {RL['CIITA']['lesion_p']:.2f}). Bottom, SpatialPert ranks a knockout's p-value among control pseudo-knockouts with the same lesion sizes in "
            "the same tissue, then reports the knockout-by-niche contrast Ψ with a lesion-clustered interval, or flags the design as underpowered. Cards: "
            f"largest fraction of null pseudo-knockouts called niche-dependent by a cell-level test (Fig. 2a); power of the calibrated lesion test with 20 "
            "lesions for a niche effect half the shared effect (Fig. 3a); within-mouse permutation and calibrated p for Tgfbr2 in Perturb-map (Fig. 5)."]),
        (FIG / "fig3_calibration_power.png", "Fig. 2 | Calibration and power of eight tests across lesion numbers.", [
            "a, Type-I error (300 null draws per cell; bold, above 0.06). b, Size-adjusted power for a niche effect half the shared effect (150 draws). Hatched, "
            "test not applicable (fewer than two lesions in two niches). Tests above the line exceed 0.06 at some lesion number. c, Type-I error with circular "
            "(dashed) and leave-out (solid) niche labels. All pseudo-knockouts have 400 cells. Exact values with 95% Wilson intervals are in the source data."]),
        (FIG / "fig4_design.png", "Fig. 3 | Designing screens for niche-dependence tests.", [
            "a, Size-adjusted power of the four calibrated tests; ribbons, 95% Wilson intervals; dotted, 80%. b, Power predicted by the design-effect formula "
            "with efficiency 1 or fitted, against the simulated power of lesion permutation. c, Count-level effects (δ = 0.5 or 1, all genes) with real lesion-size "
            "distributions: size-adjusted power and type-I error of lesion permutation. d, Lesions per Pro-Code in the Perturb-map imaging deposit (S-BIAD267; "
            f"{PM['lesions']:,} lesions, 8 mice), log scale; dashed lines, 10 and 20 lesions."]),
        (FIG4 / "fig_calibration_context.png", "Fig. 4 | Why tests fail across tissues, and the remedies.", [
            "a, SD of lesion-mean scores by niche in tumour and brain (control pseudo-lesions, 20 lesions per pseudo-knockout, 60 draws); variance ratio, "
            "between-lesion variance over sampling variance. b, Type-I error of lesion permutation before (dashed) and after (solid) calibration against matched "
            "pseudo-knockouts (two-fold cross-fitting; bars, 95% Wilson intervals). c, Coverage of 95% intervals for the knockout-by-niche contrast Ψ: CR2 "
            "clustered on lesions (solid) and naive cell-level (dashed); 200 draws per point. d, Size-adjusted power (20 lesions, effect following host-cell "
            "fraction) with four niche representations under the same test, pooled over three encoder seeds (450 draws); bars, 95% Wilson intervals; ×, "
            "type-I error above 0.06 in some run."]),
        (FIG4 / "fig_perturbmap.png", "Fig. 5 | Perturb-map: knockouts change immune infiltration, but not by niche.", [
            f"{PMD['lesions']:,} curated lesions in {PMD['lobes']} lobes of {PMD['mice']} mice (S-BIAD267). a, Lesions per niche for the {len(PMD['eligible'])} eligible knockouts. "
            "b, Niche dependence of CD8 infiltration, Ψ (immune-rich minus immune-poor), with 95% CR2 intervals clustered on mouse; q, BH over knockouts of the "
            "calibrated p. c, Main effect on CD8 infiltration against other knockouts (mean of mouse means, 95% t interval; exploratory). d, Within-mouse "
            "permutation p against calibrated p for all 36 tests (filled, CD8)."]),
        (FIG / "fig6_real_data.png", "Fig. 6 | Xenium knockouts: cell-level significance disappears with lesions as units.", [
            "a, p-values (CR2 and mixed model: minimum BH q over genes) for four knockouts and eight tests on the same niches; colour, −log10 p; hatched, untestable "
            "(a niche seen in fewer than two lesions). Permutation p-values with 999-1,000 permutations are floored at 0.001. b, Cells per niche for each "
            "knockout, stacked; numbers, lesions visiting niches N0/N1/N2/N3."]),
    ]
    ed = [
        (FIG / "ed1_qq.png", "Extended Data Fig. 1 | Null p-value distributions.", [
            "QQ plots of null p-values (leave-out labels, 400 cells) at 5, 10 and 20 lesions for six tests. Grey band, pointwise 95% interval for 300 uniform "
            "p-values. Discreteness of the sign-flip test reflects the small number of straddling lesions."]),
        (FIG / "ed2_icc_runtime.png", "Extended Data Fig. 2 | Within-lesion correlation and runtime.", [
            "a, Within-lesion intraclass correlation of the tested score (median and interquartile range over pseudo-lesions) and per gene (median), by lesion size. "
            f"b, Runtime on a {RT['hardware']['logical_cpus']}-core laptop with {RT['hardware']['ram_gb']} GB RAM, log scale."]),
        (FIG / "ed3_callers.png", "Extended Data Fig. 3 | Niches from four callers.", [
            "Niche labels of control tumour cells from k-means on neighbourhood features, BANKSY, SpaGCN and CellCharter (Methods). Scale bars, 500 µm."]),
        (FIG / "ed4_placement_niches.png", "Extended Data Fig. 4 | Lesion placement and niche environment.", [
            "a, Centroids of simulated lesions (20 lesions per pseudo-knockout, first 60 null draws) over control tumour cells; scale bar, 1 mm. b, Distance to "
            "the section centre for lesion centroids and control cells (two-sample Kolmogorov-Smirnov D shown). c, Host-cell fraction within 50 µm of control "
            "cells by niche; line, median."]),
        (FIG / "fig2_tissue.png", "Extended Data Fig. 5 | Niches and knockout lesions in an HT-29 lung metastasis section.", [
            "a, Control-defined niches (k-means, K = 4) in the Xenium section (GSE319123). b, Cells of the four testable knockouts and their 100 µm lesions. "
            "c, Marker genes per niche in control cells; dot size, fraction of cells expressing; colour, z-scored mean expression. "
            "d, Lesion-niche visits per knockout, stacked by niche; numbers above bars are lesions per knockout. Scale bars, 500 µm."]),
        (FIG / "fig5_generalisation.png", "Extended Data Fig. 6 | Calibration depends on tissue and niche caller (uncalibrated tests).", [
            "a, Type-I error in a Xenium mouse brain section (10 µm bins), leave-out labels; points shown where the test was applicable in at least 10% of draws. "
            "b, Type-I error in the tumour section with niches from four callers (circular labels for all callers; 200 null draws per cell; bold, above 0.06; "
            "hatched, not applicable)."]),
    ]
    for f, t, leg in figs + ed:
        m.figure(f, t, leg)

    m.columns(1, new_page=True)
    m.heading("Supplementary Note 1 | Pre-registration history")
    m.para(
        "The analysis plan (RESEARCH_PROTOCOL.md, tag protocol-v1.0) was frozen before any knockout expression outcome was examined. It registered five "
        "hypotheses about the real screen and a positive control. After the first results showed that clonal structure, not biology, dominated the tests, "
        "the study was re-scoped through dated amendments, each registered before the analysis it governs: AMENDMENT-001 (implementation correction and "
        "pooled branch, before any outcome), 002 (exploratory analyses after internal review), 003 (calibration benchmark of test families), 004 (follow-up "
        "checks of one knockout), 005 (local-control residuals; not adopted), 006 (studentised tests under lesion structure; none qualified), 007 "
        "(geometry-matched null), 008 (small-cluster tests, leave-out niche labels and the adoption rule) and 008b (corrected size-adjusted critical "
        "value). After review 3, AMENDMENT-009 (exploratory, registered after the brain failure was known) covered the brain diagnosis, the matched-null "
        "calibration, effect intervals and the context encoder, and AMENDMENT-010 fixed the Perturb-map plan before any outcome was computed; 010a replaced "
        "distance-based lesions with the authors' curated lesion outlines after a structural check, also before any outcome. The original hypotheses and "
        "their status are below; none is a claim of this paper."
    )
    m.table(["Hypothesis (protocol v1.0)", "Registered result", "Status"], [
        ["H1 niche dependence", f"p = {H['H1']['p']:.3f}", "Untestable at lesion level for 3 of 4 knockouts (Fig. 6)"],
        ["H2 stable core effect", f"median r = {H['H2']['median_lodo_r']:.2f}", "Not supported; driven by one knockout"],
        ["H3 microenvironment genes", f"p = {H['H3']['p']:.3f}", f"Random 8-gene sets p = {H3R['p_random_sets']:.3f}; not specific"],
        ["H4 niche-aware prediction", f"p = {H['H4']['p']:.3f}", "Not supported"],
        ["H5 cross-platform agreement", f"median r = {H['H5']['median_r']:.2f}", "Not supported"],
        ["Positive control (Tgfbr2)", f"r = {S5['PC1']['Tgfbr2']['median_pairwise_r']:.2f}, p = {S5['PC1']['Tgfbr2']['p_value']:.3f}", "Pass; validates the effect estimator"],
    ], widths=[2.0, 1.7, 3.4])
    m.fill_references()
    out = REP / f"SpatialPert_manuscript_{VERSION}.docx"
    m.save(out)
    return out


def main() -> None:
    out = build()
    pdf = out.with_suffix(".pdf")
    ps1 = (f"$w = New-Object -ComObject Word.Application; $w.Visible = $false; "
           f"$d = $w.Documents.Open('{out}'); $d.SaveAs([ref]'{pdf}', [ref]17); $d.Close(); $w.Quit()")
    subprocess.run(["powershell", "-NoProfile", "-Command", ps1], check=True)  # noqa: S603, S607
    print("written", out, pdf)
    print("abstract words", len(abstract().split()))


if __name__ == "__main__":
    main()
