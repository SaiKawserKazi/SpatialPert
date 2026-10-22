# SpatialPert V5 audit: discrepancy log, risk assessment, evidence inventory

Audit of the current manuscript (`reports/SpatialPert_manuscript_v0.0.docx` / `.pdf`, the
former v4.0; renamed at the author's request) against the brief
`SPATIALPERT_V5_CLAUDE_FEEDBACK.md`. Every statement below is traced to a file. Items marked
**[verified]** were re-derived during this audit; **[open]** items still need work.

Labels:
- **S** = source-supported fact;
- **N** = required new computational analysis;
- **I** = interpretation or editorial recommendation.

## 1. How manuscript numbers are produced (provenance)

Every number in the manuscript text and figure legends is computed at build time by
`scripts/build_paper_v4.py` (with Methods text from `scripts/paper_v3_methods.py`) from
`reports/*.json`. None is typed by hand. The figures are drawn by these scripts:
- `scripts/make_figures_v3.py`
- `scripts/make_figures_v4_review3.py`
- `scripts/make_figure_perturbmap.py`
- `scripts/make_fig1_overview.py`

Each figure has source-data CSVs under `reports/figures/v3|v4/source_data/`. The
independent recheck `scripts/recheck.py` recomputes 489 reported quantities from per-draw
CSVs or raw files (`reports/recheck.json`, 489/489 OK at the start of this audit). Its
limitation is that some checks reuse pipeline code, for example the lesion table
`build()`. It verifies summaries, not the simulations themselves; those are checked by
re-running (section 4).

## 2. Discrepancy log

| ID | Location | Current statement | Source evidence | Severity | Required action | Status |
|---|---|---|---|---|---|---|
| D1 | Results (Perturb-map), Fig. 5 legend vs Methods / Fig. 3d legend | 1,377 lesions vs 2,192 polygons | `step38_lesion_attrition.json`: 2,192 named polygons → 2,180 with a mappable Pro-Code → 1,414 with ≥ 30 labelled cells → 1,377 with ≥ 20 ring cells | Medium (S) | Report the attrition table in Methods and legends; use "polygons" for annotations and "analysed lesions" for 1,377 | Resolved; text pending |
| D2 | `step36_perturbmap_niche.py` `build()` | Lesions with < 20 cells in the context ring silently dropped (37 lesions) | Code line `ok = feats.notna().all(1) & (les.ring_n >= 20)`; rule absent from AMENDMENT-010/010a and Methods | **High (S)** — undocumented protocol deviation | Declare as deviation; sensitivity analysis without the filter (N) | Open |
| D3 | Results ("Current screens…"), Fig. 3d, step25 | "median Pro-Code … 22.5 lesions; 25 of 46 codes" | 11 of the 46 "codes" are suffix variants such as `C;NWS;V5_M1` (12 polygons) that the parser did not strip; after removal: **35 codes, median 59** lesions; codes ≥ 20 = 25, ≥ 10 = 31 (unchanged) | **High (S)** — wrong number and a figure with 11 spurious bars | Fix the parser (exclude `_M` variants as ambiguous, report them), regenerate step25 and Fig. 3d | Open |
| D4 | Perturb-map analysis | 12 `_M`-suffixed polygons excluded without mention | Unmapped codes in `step38_lesion_attrition.json` | Low (S) | Document; the meaning of the suffix is unknown, so exclusion is the conservative choice | Open |
| D5 | Results, Fig. 5 legend | "smallest calibrated q-value over 36 tests was 0.57" | `step36_perturbmap.json`: BH applied **within each outcome** (12 knockouts × 3 families), not across 36 | Medium (S) | State the family: per outcome, CD8 primary; report q across all 36 as sensitivity (N, trivial) | Open |
| D6 | Abstract, Results, Fig. 1 | Calibration against "control pseudo-knockouts" | True for Xenium and the simulations. In Perturb-map there are no non-targeting controls, so the pool is **other-knockout lesions** (Methods states this; Results and Fig. 1 imply controls) | **High (I)** — the validity contract differs: other knockouts may themselves have niche effects, so the pool can be contaminated | State the Perturb-map pool explicitly and its extra assumption; sensitivity with pools excluding knockouts with main effects (N) | Open |
| D7 | Title, abstract, intro, results | "clonal lesions", "clones", "Clonal structure confounds…" | Lesions are distance-based connected components (100 µm single linkage), or curated polygons in Perturb-map; no lineage barcodes | Medium (I) | Use "spatial lesions" (operational definition); keep "clonal growth" only as the biological motivation, with a caveat | Open |
| D8 | Abstract, intro, Discussion, Fig. 1 card | "about 20 lesions per knockout are needed"; "20 lesions → 83% power" | 83% holds for one effect model (λ = 0.5, log scale, tumour, 400 cells). Count-level effects give ≤ 0.19 at 20 lesions (`step26`); brain gives ≤ 0.11 after calibration (`step31`) | **High (I)** — over-general design rule | Replace with a conditional statement plus an MDE/power grid (P3) | Open |
| D9 | Abstract, intro | Calibration, intervals and encoder presented alongside pre-registered results | AMENDMENT-009 registered **after** the brain failure (post hoc, exploratory) | Medium (I) | Label as exploratory method development in abstract and results | Open |
| D10 | Results ("Calibration…"), Fig. 4b | "calibrated type-I 0.043–0.065 in both tissues" | Cross-fitting splits draws from the **same** generator as the pool, so the evaluation is exchangeable by construction (Scenario A only); no test under non-exchangeable placement | **High (I/N)** — validity shown only where it holds by design | Scenarios B–H (P1.2); narrow the claim | Open |
| D11 | All type-I estimates | 200–300 null draws (Monte Carlo SE ≈ 0.013 at 0.05) | `step19`, `step22`, `step31` | Medium (N) | ≥ 2,000 draws for the key validity claims where feasible, else report the uncertainty | Open |
| D12 | Benchmark (mixed model) | "type-I up to 100%" | statsmodels MixedLM Wald test on lesion × niche pseudobulk, 40 null draws per level (registered deviation); no Satterthwaite/Kenward-Roger (R unavailable) | Medium (I) | Describe as "the tested implementation"; do not generalise to mixed models | Open |
| D13 | Benchmark (CR1) | "CR1 standard errors up to 94%" | `baselines.cluster_robust_any`: per-gene CR1 Wald F(k−1, G−1), BH, knockout rejected if any gene q < 0.05; independent-clone design, 500 null draws | Low (S) | State the test family; CR1 is a known liberal small-cluster estimator, and CR2 was also evaluated | Open |
| D14 | Encoder (Fig. 4d) | "roughly doubled power (8% to 17%)" | Encoder trained on **all** control cells, including those later drawn as pseudo-knockouts (their context is left out at embedding time, but the weights saw them); no train/validation/test split; one effect model; three seeds | Medium (N) | Spatial train/test split of the encoder (P2.4); keep it out of the headline unless it holds | Open |
| D15 | Callers (Ext. Data Fig. 6b) | caller comparison | SpaGCN on a 6,000-cell subsample with k-means initialisation and kNN propagation; CellCharter without scVI; circular labels for all callers | Low (S) | Already stated; keep "tested implementations" wording | Documented |
| D16 | Methods (Perturb-map) | "Pixel size is taken as 0.5 µm" | Not verified from the deposit metadata | Medium (S) | Verify from the S-BIAD267 metadata or the original paper; else keep as an assumption and give radii in pixels (already in protocol) | Open |
| D17 | Results (Xenium real data) | CIITA p = 0.003 (cell) vs 0.43 (lesion) | `step28_real_smallcluster.json`; niche labels from NicheModel (k-means K = 4), which differ from the Step-4 domains used for the original H1 (Supplementary Note) | Low (S) | State which niche definition each result uses | Open |
| D18 | Front matter, availability | Author list, affiliations, GitHub URL, licence, DOI, Zenodo | Placeholders in the manuscript | **Blocking (S)** | Author-only items; must not be invented | Open (author) |
| D19 | Recheck | "independent recheck" | Some checks reuse pipeline code (`step36.build`) | Low (I) | Say "recomputes summaries from per-draw outputs"; re-runs give simulation-level reproduction | Open |
| D20 | Real-data inference (Xenium) | calibrated p for Xenium | Pool of 200 pseudo-knockouts; only 72.5% (CIITA, lesion) testable; others untestable | Low (S) | Report the testable fraction and abstentions alongside p | Open |

## 3. Scientific risk assessment (from the brief, checked against the files)

1. **Pseudo-knockout exchangeability (high).** The calibration validity rests on knockout
   lesions being placed like random control positions (simulations, Xenium), or like other
   knockouts' lesions (Perturb-map). Nothing in the current files tests a violation (D10).
   Perturbation-dependent survival or placement is plausible in tumours.
2. **Semi-synthetic realism (medium).** Pseudo-knockouts are control cells. No lineage, no
   fitness effects, and injected effects only on the log scale or as count thinning/boosting.
3. **Lesion definition (medium).** Single linkage at 100 µm, untested against known lineage.
   A threshold sensitivity exists only for patches in real data (v2), not for the calibrated
   test.
4. **Low power and imbalance (high).** Perturb-map has 5–19 immune-rich lesions per knockout;
   brain power is ≤ 0.11 after calibration. The design rule is over-general (D8).
5. **Post hoc development (medium).** AMENDMENT-009 is exploratory (D9). The Perturb-map
   plan (010/010a) was committed before outcomes; this is verified by git history (commit
   order of `RESEARCH_PROTOCOL.md` vs `step36_perturbmap.json`).
6. **Real-data interpretation (medium).** A minimum q of 0.57 and intervals crossing 0 mean
   "not detected at this power", not absence. The current text mostly says this; the Fig. 5
   panel title "none at q < 0.1" is fine.
7. **Benchmark fairness (medium).** MixedLM without small-sample df (D12); SpaGCN subsampled
   (D15); CR1 is a liberal baseline (D13); WCB with 299 draws.
8. **Reproducibility gaps (blocking for submission).** Author metadata, repository URL,
   licence and DOI are missing (D18). Seeds are fixed in every script. The environment is
   recorded in `reports/v5_environment.txt` and `reports/v5_environment_freeze.txt`.

## 4. Reproduction of the v4.0 (= v0.0) baseline (P0.3)

Re-runs use the committed scripts and seeds. Outputs are compared with the committed files
in git. See section 4.1 (filled after the runs complete).

### 4.1 Results

Event during the audit: at 21:11 on 2026-10-22 the project's git repository was
re-initialised outside this session (empty history, remote
`github.com/SaiKawserKazi/SpatialPert`). The committed baseline used for diffing was
therefore unavailable. Four re-runs (steps 27, 28, 31, 37) had already overwritten their
outputs. They were compared with values recorded before the re-run: this session's earlier
console output, the v0.0 manuscript text, and `recheck.json` written at 20:52. All four
reproduce exactly. A snapshot of all results was then taken (`v5/baseline_v0.0/`), and later
re-runs are diffed against it. The author now manages git; no commits are made by the audit.

| Analysis (script) | Manuscript use | Result | Runtime (this PC) |
|---|---|---|---|
| Design-effect formula (step27) | Fig. 3b, Results | Exact (power 0/0.02/0.367/0.827; RMSE 0.325 → 0.269) | < 1 min |
| Xenium small-cluster + permutation (step28) | Fig. 6 | Exact (CIITA cell p 0.003, lesion p 0.4306, CR2 q 0.123) | ~2 min |
| Calibration evaluation (step31) | Fig. 4b, Results | Exact (24/24 raw rates; calibrated 0.043–0.065; brain power ≤ 0.113) | < 1 min |
| Perturb-map integrity (step37) | Methods | Exact (9/9) | ~6 min |
| Coverage, tumour (step33) | Fig. 4c | Exact (48/48 values) | ~5 min |
| Encoder benchmark seed 0 (step35) | Fig. 4d | Running | — |
| Perturb-map analysis (step36) | Fig. 5 | Running | ~60 min |
| Stage-1b benchmark, L = 10 (step19) | Figs 2–3 | Running (row-level comparison with `step19_stage1b.csv`) | ~1.5–2 h |
| Pro-Code counts (step25) | Fig. 3d, Results | **Changed on purpose** (audit D3 fix): 35 codes, median 59, 2,180 unambiguous + 12 ambiguous polygons | < 1 min |

Not re-run (cost; summaries are verified from per-draw CSVs by the recheck): steps 19 (other
L), 22, 24, 26, 29, 32, 34. Estimated full re-run ≈ 20 CPU-hours on this PC.

## 5. Evidence inventory

The machine-readable inventory is `v5/v5_claims.csv`: one row per numeric claim, with
location, value, source file/key, script, figure, recheck item and reproduction status.
It is written after the reproduction runs.
