# SpatialPert research protocol (v1.0, frozen 2026-10-20)

**Question.** Does the transcriptional effect of a gene knockout depend on where in a tumour
the knocked-out cell sits, beyond what its local cellular neighbourhood explains?

Frozen after Gates 1–3, which inspected guide counts only. No expression outcome of any
knockout was examined before this freeze. Every change after the freeze is a numbered amendment.

## 1. Data and roles

| Role | Dataset | Unit |
|---|---|---|
| Discovery | Xenium HT-29 lung metastasis, GSE319123 (Perturb-DBiT paper) | single cell |
| Replication (cross-platform) | DBiT-seq HT-29 lung, GSE319277: Bru.LC.6, Bru.LC.8, Bru.LC.1.33.3 | 20–50 µm pixel |
| Positive control | Perturb-map Visium, GSE193460: Tgfbr2 (3 sections), Ifngr2 (2) | 55 µm spot |

Perturb-FISH (Binan 2025) is dropped: its tumour tables are not public (Gate 2).

## 2. Pre-processing (Xenium)

1. Remove guide probes, control codewords and cells with < 20 total counts from the
   expression matrix. Normalise to the median total, then apply log1p.
2. **Tumour cells.** Tumour score = mean z-score of EPCAM, KRT8, KRT18, KRT19, KRT20, CDX2
   and CEACAM5 (whichever are in the panel). Fit a 2-component Gaussian mixture on the score.
   Tumour cell = posterior of the high component > 0.9. All effects are estimated in tumour
   cells only.
3. **Knockout calls.** A tumour cell is KO-k if guide probe k has ≥ 1 UMI and no other guide
   probe has any. Cells with 2 or more different guides are excluded. Control = tumour cells
   with no guide UMI. There are no non-targeting probes, so controls may carry undetected
   guides. This biases effects towards zero (conservative).
4. **Analysis genes.** Panel genes with mean normalised expression ≥ 0.1 in control tumour
   cells.

## 3. Spatial domains

Domains are called **without perturbed cells** so that a knockout cannot define its own
context.

1. Compute 30 PCs on control tumour cells.
2. **Neighbourhood features.** For every tumour cell, take the SpaGCN kernel-weighted mean of
   the PCs of its control tumour neighbours within 100 µm, **excluding the cell itself**.
   Weight: w = exp(−d²/2l²), with l = 50 µm.
3. Run k-means on the neighbourhood features of control cells, K ∈ {3, 4, 5, 6}. Choose K
   by silhouette on a random 20,000-cell subsample, seed 0.
4. Assign every tumour cell (control or KO) to its nearest centroid using its own
   neighbourhood features.

**Covariates per cell:** log total counts; host fraction (share of non-tumour cells within
50 µm); local cell density (cells within 50 µm).

**Full gate (first check of Step 4).** A knockout is *testable* if it has ≥ 30 KO tumour
cells in each of ≥ 2 domains.
- 10 or more testable: proceed as planned.
- 5–9: proceed, with H1 labelled as limited power.
- Fewer than 5: H1–H4 are reported for the pooled set of all knockouts only.

## 4. Effect model

For cell i in domain d, the residual is r_i = y_i − μ̂_d − β̂ᵀc_i. Here μ̂_d is the
control-cell mean in domain d, and β̂ comes from regressing control cells on the covariates
c_i. The domain effect of knockout k is e_kd = mean of r_i over KO-k cells in d, giving one
vector over the analysis genes.

**Clone patches.** KO-k cells joined by single linkage at 30 µm form patches. A patch takes
the majority domain of its cells.

## 5. Hypotheses and decision rules

**H1. Location dependence.**
- Statistic: T_k = Σ_g Σ_d n_kd (e_kdg − ē_kg)², over testable domains.
- Null: permute domain labels **among patches** of KO-k (residuals fixed), 1,000
  permutations. Benjamini–Hochberg across testable knockouts.
- Supported if FDR < 0.05 for ≥ 20% of testable knockouts.
- Cell-level permutation is reported as a sensitivity analysis only.

**H2. Stable core.**
- For each knockout and held-out domain d, compute the Pearson r over genes between e_kd and
  the pooled effect from the other domains.
- Supported if the median r ≥ 0.3 and exceeds the same statistic for size-matched random
  control patches (pseudo-KOs, 200 draws) with one-sided p < 0.05.

**H3. Microenvironment genes vary more.**
- Per-gene interaction energy Σ_d n_kd (e_kdg − ē_kg)², summed over testable knockouts.
- Compare genes in MSigDB Hallmark IFN-γ response, IFN-α response, TNF-α via NF-κB,
  hypoxia and EMT against all other analysis genes.
- Supported if one-sided Wilcoxon p < 0.05.
- Gene sets come from Enrichr `MSigDB_Hallmark_2020`, downloaded once and committed.

**H4. Domain-aware prediction.**
- Divide the section into a 4 × 4 grid of spatial blocks. For each block and testable
  knockout with ≥ 15 KO cells in the block, predict the block's effect vector two ways:
  (a) the pooled effect from other blocks, and (b) the effect from other blocks in the same
  domain.
- Supported if the paired median of r(b) − r(a) > 0 with Wilcoxon p < 0.05.

**H5. Cross-platform agreement.**
- For knockouts with ≥ 30 KO cells in Xenium and ≥ 30 guide-positive pixels in a paired DBiT
  section, correlate the pooled Xenium effect with the DBiT pixel effect over shared genes.
  DBiT pixels are normalised the same way. DBiT controls are tissue pixels without any guide.
- Supported if the median r exceeds the random-gene-permutation null (1,000 draws) with
  p < 0.05.

## 6. Controls that must pass before H1–H5 are opened

- **PC1, estimator (Perturb-map).** Tgfbr2-vs-unassigned-KP-tumour effect vectors (Visium
  spots, same residual method, section as domain) must agree across sections: median pairwise
  Pearson r ≥ 0.3, above a shuffled-label null (p < 0.05).
- **PC2, calibration (Xenium).** Run the H1 test on 200 pseudo-KOs: control patches matched
  in size and domain distribution to a real knockout. The false-positive rate at
  α = 0.05 must be ≤ 0.08 (the binomial upper limit for 200 draws).

If PC1 or PC2 fails, analysis stops and the pipeline is fixed under an amendment. No
hypothesis result is examined first.

## 7. Deviations from the proposal (2026-10-20)

| Proposal | Protocol | Reason |
|---|---|---|
| Discovery: Perturb-map | Xenium HT-29 (GSE319123) | Gate 1: 2/35 KOs labelled |
| Replication: Perturb-DBiT; external: Perturb-FISH | DBiT-seq replication; Perturb-FISH dropped | Gate 2: tables not public |
| Leave-one-slide-out (H4) | 4 × 4 spatial blocks | Single Xenium section |
| H3 classes: TME vs intrinsic genes | Hallmark microenvironment gene sets | Brunello hits have no pre-set classes |
| H5 vs dissociated Perturb-seq | Cross-platform Xenium vs DBiT | No dissociated screen of these KOs in HT-29 |
| Negative control: non-targeting guides | Pseudo-KO control patches | No NTC probes in the Xenium panel |

## 8. Software

Python 3.12, `uv`. Gate: ruff + mypy (strict) + pytest via `status.ps1`. Seeds fixed at 0.
Every result is written to `reports/` by a script, and `scripts/recheck.py` recomputes the
headline numbers independently.

## AMENDMENT-001 (2026-10-20, before any knockout expression outcome was examined)

1. **Implementation correction.** The first Step 4 run filtered analysis genes on
   expm1(mean log expression), not on mean normalised expression as §2.4 states. Corrected:
   74 analysis genes, K = 4 (silhouette 0.083), 4 testable knockouts (CIITA, DIO3, ESAM,
   PRSS23). The pre-set rule (< 5 testable) selects the **pooled branch**. The buggy run
   (5 testable, "limited" branch) is superseded and not used.
2. **Definition of the pooled branch.** The confirmatory H1 statistic is T = Σ_k T_k over
   the testable knockouts. Its null distribution sums the per-knockout patch-permutation
   statistics drawn independently (1,000 draws). H2 and H4 pool knockout-domain (H2) or
   knockout-block (H4) units across the testable knockouts. Per-knockout results are
   reported as exploratory only.
3. **New Aim 6: semi-synthetic benchmark (known truth).** Using real control tumour cells
   and the real domain map:
   - Draw pseudo-knockouts with the patch structure of the 4 testable knockouts, scaled to
     n ∈ {50, 100, 200, 400, 800} cells.
   - Add a known effect: a common core (one random direction over the analysis genes, scaled
     to effect size δ) plus a domain-specific part (λ·δ in one domain, λ ∈ {0, 0.25, 0.5, 1}).
   - Run 100 replicates per setting with seed-fixed draws.
   - Compare three tests at α = 0.05: (a) patch permutation (ours), (b) cell-level
     permutation (naive), (c) per-gene ANOVA + BH (any gene q < 0.05).
   - Report type-I error (λ = 0) and power (λ > 0), and the cells per knockout needed for
     80% power.
4. **H5 scope.** Testable overlap between Xenium and paired DBiT sections is decided by the
   §5 H5 rule. Knockouts below it are reported descriptively.

## AMENDMENT-002 (2026-10-21, after results-v0.1; all items exploratory)

Prompted by review round 1 (review/REVIEW_ROUND1.md). Confirmatory H1-H5 results stay as
reported. Everything below is labelled as exploratory or sensitivity analysis.

1. Benchmark variants (same settings as Aim 6):
   - (A) control model with the first 10 neighbourhood PCs added as covariates;
   - (B) permutation among 250 µm tiles.
   A variant is adopted for the software default only if its type-I error is ≤ 0.08 at
   every n. If both qualify, the one with higher mean power wins.
2. H1 re-run with the adopted variant (if any), and for K = 3, 5, 6.
3. Lesion counts at 100 µm linkage for each testable knockout.
4. H2 on unique units; H3 with variance-scaled energy; H5 with tumour-rich control pixels.
5. Tumour-call sensitivity at posterior thresholds 0.5 and 0.99 (counts only).

## AMENDMENT-003 (2026-10-21, after the external review "SpatialPert Review vs SpaGCN"; exploratory)

Registered before running any of the analyses below. Confirmatory H1-H5 stay as reported.

**Stage 1. Calibration.**
- Benchmark nulls (λ = 0) at 50-800 cells, 500 replicates, recording the number of patches.
- Methods compared:
  1. patch permutation (current);
  2. studentised patch permutation: per-gene Welch statistics on patch means, summed;
  3. patch pseudobulk per-gene ANOVA + BH (a pseudobulk analogue of the random-intercept LMM);
  4. per-gene OLS with patch-clustered sandwich SE (CR1) + BH;
  5. 250 µm tile permutation;
  6. cell-level permutation;
  7. per-gene cell-level ANOVA.
- Adoption rule: the method with type-I error ≤ 0.06 at every n and the highest mean power
  over λ ∈ {0.25, 0.5, 1} becomes the default. If none qualifies, the current default stays
  and the paper is reframed as a benchmark plus design guidance.

**Stage 2. Making the null informative.**
- Lesion number L ∈ {1, 3, 5, 10, 20} as a benchmark factor at a fixed 400 cells.
- λ estimate with CI and minimum detectable λ for each real knockout.
- Niche stability (bootstrap ARI) and niche biology.
- A continuous host-fraction interaction test.
- Guide calls at ≥ 2 UMI; guide spillover versus distance; on-target check (CIITA is the
  only testable target on the panel).

**Stage 3. A real positive.**
- Same real data analysed with the naive cell-level test, to show where conclusions differ.
- Spike-in aligned with a real niche-varying programme (hallmark hypoxia genes on the panel).
- Author emails drafted, not sent.

## AMENDMENT-004 (2026-10-21, exploratory follow-up of a review sensitivity result)

With ≥ 2-UMI guide calls, DIO3 gave H1 p = 0.011. In ≤ 200-cell subsamples, 65% of draws had
p < 0.05 (null ≈ 3.5%). DIO3 is a single lesion, so this is at most within-lesion niche
dependence. Before interpreting it, these checks are fixed:

1. ≥ 3-UMI calls (same subsampling).
2. Distance to the nearest host cell added as a covariate (rim versus core).
3. Depth: log total counts of DIO3 ≥ 2-UMI cells per niche, against controls in the same
   niches. Repeat with cells matched on log total count (±0.1) to controls.
4. A specificity control: the same ≥ 2-UMI analysis for PRSS23, which has three lesions.
5. Genes driving the niche difference (the top 10 by standardised between-niche difference),
   reported descriptively.

The finding is reported as exploratory in all cases.

AMENDMENT-004b (before running): check 5 showed that the largest DIO3 niche differences come
from niche 3, which has 7 cells after patch-majority assignment. Check 6: repeat the ≥ 2-UMI
subsample test with niche-3 cells removed (niches 0 vs 2 only). If the signal disappears, the
DIO3 result is reported as driven by a 7-cell niche and is not interpreted.

## AMENDMENT-005 (2026-10-21, before running; method revision after the lesion benchmark)

**Result motivating this amendment** (reports/step12_lesions.json): with lesion-structured
pseudo-knockouts, the default 30 µm patch test has type-I error 0.11-0.28 at every lesion
number (400 and 800 cells). Lesion-level units hold 0.02-0.08 for L ≤ 10 but 0.18 at L = 20.
The likely cause is lesion- and region-specific baselines that tissue-wide niche means do not
remove.

**New residual definition: local-control residuals.**
- Remove covariate effects with the control model's slopes: z = y − β'(c − centre).
- For each knockout cell, subtract the mean z of its k = 50 nearest control cells. The
  knockout set's own cells are excluded from the controls (in the benchmark the pseudo-
  knockout cells are controls, and are excluded).
- Niche means are then no longer used. Niche labels are used only to group cells for the
  test.

**Evaluation** (fixed before running):
- the lesion design (L ∈ {1, 3, 5, 10, 20}, n ∈ {400, 800}, 300 null replicates,
  λ ∈ {0.5, 1}, 150 replicates);
- the independent-patch design (n ∈ {50, 200, 800}, 300 null replicates).

The units are 30 µm patches and lesions (100 µm).

**Adoption rule:** a local-residual variant becomes the default if its type-I error is
≤ 0.06 in every null cell of both designs. Otherwise the paper reports the method as
uncalibrated under lesion structure and is reframed as a benchmark plus design guidance.

## AMENDMENT-006 (2026-10-21, before running)

**Outcome of AMENDMENT-005:** local-control residuals were **not adopted**. Patch30 type-I
error was 0.22-0.48 under lesion structure; lesion units were 0.03-0.22.

**Outcome of the Stage-1 rule (AMENDMENT-003):** the studentised patch permutation qualified
in the independent-patch design (max type-I 0.023), but it was never tested with lesion
structure.

**Final candidate, evaluated in the lesion design** (L ∈ {1, 3, 5, 10, 20}, n ∈ {400, 800},
300 null and 150 λ ∈ {0.5, 1} replicates, global residuals):
- (a) studentised permutation with 30 µm patch units;
- (b) studentised permutation with lesion units (100 µm);
- (c) plain lesion units, for reference.

**Adoption rule:** type-I error ≤ 0.06 wherever the test is applicable. Cases where it is not
applicable (fewer than 2 units per niche) are reported as untestable, not counted.
- If (b) qualifies, the default becomes lesion-level studentised permutation, and knockouts
  with fewer than 2 lesions per niche are declared untestable.
- If nothing qualifies, the paper is reframed as a benchmark plus design guidance.

## AMENDMENT-007 (2026-10-21, before running)

**Outcome of AMENDMENT-006:** no candidate qualified. Max type-I error under lesion
structure: studentised patch 0.33, studentised lesion 0.19, plain lesion 0.23. Real data: at
lesion level, no knockout has ≥ 2 lesions per niche (studentised untestable for all 4).
Patch-level results depend on linkage (20 µm: DIO3 p = 0.001; 30 µm: p = 0.79; 50 µm:
untestable). Per the AMENDMENT-006 rule, the paper is reframed as a benchmark plus design
guidance.

**Last method candidate: geometry-matched null (GMN).**
- Statistic: between-niche sum of squares of knockout-cell residuals, with cells grouped by
  their own niche.
- Null: the same statistic for M = 500 pseudo-knockouts. Each pseudo-knockout has exactly
  the knockout's lesion geometry (number of lesions and cells per lesion), with lesion seeds
  placed at random among control tumour cells at least 400 µm apart.
- p = (1 + #null ≥ obs) / (1 + M).
- Assumption, stated in the paper: under H0, knockout lesions are seeded independently of
  tissue position.

**Evaluation:**
- lesion design, L ∈ {1, 3, 5, 10, 20}, n = 400: 200 null replicates and 150 per
  λ ∈ {0.5, 1};
- one shared null pool of 500 per setting.

**Adoption:** type-I error ≤ 0.06 at every L. If adopted, it is applied to the real
knockouts with their exact lesion geometry. GMN is reported alongside the benchmark
framing, not instead of it.

**Outcome of AMENDMENT-007:** GMN type-I error was 0.085, 0.035, 0.045, 0.070 and 0.050 for
L = 1, 3, 5, 10 and 20 (200 replicates each). The 0.085 exceeds 0.06, so GMN is **not
adopted**. Real-data GMN p-values (CIITA 0.27, DIO3 0.001, ESAM 0.016, PRSS23 0.023) are
exploratory only. Pseudo-knockouts are control cells that helped define the niches, while
real knockout cells did not; this asymmetry alone can produce differences.

**Final framing:** a benchmark plus design guidance. SpatialPert is presented as a
diagnostic toolkit (testability, lesion counts, linkage sensitivity, simulation on the user's
own tissue), not as a calibrated test for few-lesion designs.

## AMENDMENT-008 (2026-10-21, response to "SpatialPert v2.0 Nature Methods Review"; before running)

**A2: circularity.** In every benchmark replicate, pseudo-knockout cells are removed from the
reference used to compute their own niche features. Each of their niche labels is then
re-predicted from the remaining controls, with fixed centroids. Real knockout cells are never
in the reference, so this makes pseudo-knockouts match the real setting. The same draws are
analysed with the circular labels and the leave-out labels. Inflation is also reported
against cells per lesion. Seed positions for L = 20 are mapped and compared with
unconstrained placement (min_sep = 0).

**A3: size-adjusted power.** For each setting, the critical value is the 5th percentile of
the test's null p-values in the matching null run; power is P(p ≤ critical value) under the
alternative. Raw and size-adjusted power are both reported, with Wilson CIs. Abstract
numbers use the size-adjusted values.

**A1: small-cluster methods,** with the lesion (100 µm) as the cluster:
- (i) wild cluster bootstrap of the summed per-gene Wald statistic, restricted model, Webb
  6-point weights, 499 draws;
- (ii) CR2 (Bell-McCaffrey) with Satterthwaite df for the contrast between the two largest
  niches, per gene, BH, any gene q < 0.05;
- (iii) a within-lesion paired test: for lesions spanning two niches, the niche contrast of
  mean residuals; sum over genes of the squared mean contrast; exact or 2,000 random
  sign-flips;
- (iv) statsmodels MixedLM on lesion × niche pseudobulk, lesion random intercept, Wald test
  per gene, BH (no Satterthwaite df available without R; stated).

Lesion-design nulls use L ∈ {3, 5, 10, 20}, 400 cells, 300 replicates; power uses λ = 0.5,
150 replicates.

**Adoption rule:** type-I ≤ 0.06 at every L where the test is applicable. Among qualifiers,
the highest mean size-adjusted power becomes the SpatialPert default.

**A7.** A design-effect model, effective n = L·m / (1 + (m − 1)ρ), with ρ = the
within-lesion intraclass correlation of the residual summary, estimated from control
pseudo-lesions. It is compared with the size-adjusted simulated power. `plan_screen()`
implements it.

**A8.** Count-level effect injection (binomial thinning or Poisson boosting of raw counts,
then the same normalisation); δ ∈ {0.25, 0.5, 1}; sparsity 5%, 20% or 100% of genes;
lesion sizes resampled from the real knockouts' lesion-size distribution. Run with the
adopted test (or lesion permutation if none qualifies).

**A9.** Runtime and memory of `location_dependence` (1,000 permutations) on the full analysis
set, with hardware stated.

**A6.** Niches from pybanksy, SpaGCN and CellCharter (installed from PyPI in isolated
environments). If any cannot run, this is stated. The core null and design curve are
repeated with each caller.

**A4.** Pseudo-knockout lesion benchmark in a second, non-tumour public single-cell spatial
dataset (no guides needed), plus the DBiT pixel data. Data requests to the Perturb-FISH and
Spatial Perturb-seq authors are drafted for the authors to send.

**A5.** Perturb-map lesion-level IFN analysis if lesion identities are obtainable; otherwise
stated as not feasible with public data.

**AMENDMENT-008 deviation (before any MixedLM benchmark result):** statsmodels MixedLM takes
about 26 s per pseudo-knockout (74 fits). Its full replication would take about 13 h. MixedLM
is therefore run on the first 40 null and 20 λ = 0.5 draws per lesion number. One timing draw
(L = 3, null) returned p = 6e-10 and was not used.

**AMENDMENT-008b (2026-10-21; correction of the A3 definition before applying the adoption
rule).** The A3 rule set the critical value to the 5th percentile of null p-values. For
conservative, discrete tests that percentile exceeds 0.05. An example is the sign-flip test
with few lesions: with 3 paired lesions the smallest p is 0.25. Using it would raise α above
0.05 and credit such tests with power they cannot have. Size adjustment is meant only to
correct over-sized tests.

**Corrected critical value:** min(5th percentile of null p, 0.05). The adoption rule (type-I
≤ 0.06 at every applicable L; highest mean size-adjusted power) is applied with the
corrected power. The literal (uncorrected) result is also reported: the sign-flip test, mean
0.40.

## AMENDMENT-009 (2026-10-21, review 3 "Can SpatialPert be published in Nature Methods?"; exploratory, registered after the brain inflation was known)

Everything in this amendment is **post hoc and exploratory**: it was planned after seeing that
every test inflates in the mouse brain (step22). No result below changes the AMENDMENT-008
adoption.

**P1a. Diagnosis of the brain inflation** (`step29`, `step29b`; 200 null draws per L).
Candidate causes tested on the same draws: (i) the control model includes the pseudo-knockout
(refit without it); (ii) leave-out relabelling; (iii) circular labels; (iv) lesions that
straddle niches (drop minority-niche cells; residualise against the majority niche);
(v) spatial baselines (local-control residuals); (vi) unequal spread (Welch studentisation).
Outcome: (i)-(iv) explain nothing (brain relabels ~1% of bins; minority cells ~2%; L = 20
type-I 0.33-0.36 under all of them). Pooled lesion means are centred in every niche
(|t| ≤ 0.4) but their between-lesion variance is 9× the tumour's and 7.5× the sampling
variance, and their SD differs 2.6-fold between niches; lesion means are correlated up to
~500 µm. Diagnosis: lesion baselines are large, spatially correlated and heteroscedastic
across niches, so lesion means are not exchangeable between niches.

**P1b. Remedies evaluated.**
- Matched-pseudo-knockout calibration (`spatialpert.calibration`): the observed p-value is
  ranked among control pseudo-knockouts with the knockout's lesion geometry (randomised ties).
  Exact if knockout lesions are placed like random control positions (stated, untestable).
  Cross-fitted evaluation (`step31`): type-I 0.043-0.067 (brain), 0.043-0.077 (tumour).
- Diagnostic `lesion_exchangeability`: variance ratio / niche SD ratio from control
  pseudo-lesions, tumour 1.18 / 1.19, brain 6.9 / 1.7 (`step32`).
- Effect estimation (`smallcluster.niche_contrast`): Psi = Delta_b − Delta_a with CR2
  lesion-clustered intervals; coverage and bias in `step33`.

**P2. Graph context encoder** (planned): self-supervised graph encoder of control-cell
neighbourhoods; compared with fixed kernel features, expression-only and coordinate-only
representations under one downstream estimator, with effects driven by a continuous
neighbourhood covariate (not by any niche partition), so no representation is favoured by
construction.

**P5. Perturb-map imaging** (S-BIAD267 segmentation tables, per-cell Gene labels, 8 mice):
lesion-level niche dependence of immune infiltration, with lesions nested in lobes and mice.
Analysis plan to be fixed in AMENDMENT-010 before any knockout outcome is examined.

## AMENDMENT-010 (2026-10-21; Perturb-map imaging analysis plan, fixed before any knockout outcome was computed)

Only the table structure has been inspected (columns, Gene/Pro-Code counts, coordinate range).

**Data.** S-BIAD267 Lymphoid-panel segmentation tables, all available mouse-lobes (8 mice).
Per cell: centroid (px), cell area, mean intensities of CD4, CD8, B220, C and the five Pro-Code
epitopes, procode, Gene. Pixel size is taken as 0.5 µm (assumption; all radii are also stated
in px). Lobes that fail to download are listed, not imputed.

**Lesions.** Cells with a Gene label, single linkage at 60 px (30 µm) within each lobe. A lesion
is kept if it has ≥ 30 labelled cells and its majority Gene covers ≥ 80% of them; its knockout
is that Gene. No non-targeting control exists in the panel, so each knockout is compared with
the lesions of all other knockouts (reference lesions).

**Outcome (primary).** CD8 infiltration: fraction of CD8-positive cells among unlabelled cells
within 40 px (20 µm) of any lesion cell (halo + interior). CD8-positive = CD8 intensity above
the 95th percentile of all cells of that lobe (per-lobe threshold, fixed before outcomes).
Transformed by logit with 0.5/(n + 1) continuity. Secondary outcomes: CD4 and B220 (same
definition).

**Niche (lesion context).** Composition of the ring 200-600 px (100-300 µm) around the lesion,
excluding every cell within 200 px of any labelled cell (so no lesion contributes to its own
context): fractions of CD4-, CD8- and B220-positive cells and cell density. Niches = k-means
(K = 2) on standardised ring features of all lesions, fitted once; the niche with higher
B220+CD4+CD8 fraction is called "immune-rich".

**Estimand and test.** For knockout k, each of its lesions gets a residual e = outcome −
mean outcome of reference lesions in the same mouse and niche. Psi_k = mean(e | immune-rich) −
mean(e | immune-poor). Test 1: permutation of niche labels among knockout-k lesions within
mouse (9,999 permutations). Test 2 (calibration): matched null from 1,000 pseudo-knockouts made
of reference lesions with the same per-mouse lesion counts; calibrated p by randomised rank.
Interval: CR2 clustered on mouse (Imbens-Kolesár df). Main effect Delta_k (all niches) reported
as a sanity check.

**Eligibility.** Knockouts with ≥ 20 lesions, ≥ 5 in each niche and lesions in ≥ 3 mice.

**Multiplicity and claims.** Primary family: all eligible knockouts, primary outcome, calibrated
p, BH at q = 0.10. Secondary pre-specified subset: IFN pathway (Jak2, Stat1, Stat2, Ifnar,
Ifngr2, Irf1, Irf9, Ptpn2, Socs1). A knockout is called niche-dependent only if q < 0.10 and its
CR2 interval excludes 0. All other results are descriptive.

**AMENDMENT-010a (2026-10-21, before any outcome; structural check only).** On the first four
lobes the 60 px single-linkage rule gave 49 lesions where the authors' curated annotations
list ~56 per lobe; linkage fragments or drops lesions. Lesions are therefore the authors'
curated lesion polygons (`*_Lymphoid_<lobe>_curated.fix.json`, same pixel space as the
segmentation tables): lesion cells = labelled cells inside the polygon; knockout = the Gene
that the polygon's Pro-Code maps to in the segmentation tables (majority over all cells with
that Pro-Code). A polygon is kept if it contains ≥ 30 labelled cells. Outcome, niche, test,
eligibility and multiplicity are unchanged.
