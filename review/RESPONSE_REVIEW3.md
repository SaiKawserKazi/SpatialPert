# Response to "Can SpatialPert be published in Nature Methods?" (current manuscript v0.0)

All work ran on the authors' laptop (Intel i5-7200U, 2 cores / 4 threads, 8.5 GB RAM, no
CUDA GPU, Windows 10). Every number below is read from `reports/*.json`. The independent
recheck (`scripts/recheck.py`) passes **489/489** checks of the numbers in the current version
(checks of superseded versions were removed with those versions). Manuscript: `reports/SpatialPert_manuscript_v0.0.docx` / `.pdf` (Analysis
format; 139-word abstract; main text 2,858 words; Online Methods 3,425 words; 6 main and
6 Extended Data figures).

**Status of claims.** The brain diagnosis, the calibration, the effect intervals and the
context encoder were developed **after** the brain failure was known. They are registered as
exploratory (AMENDMENT-009). The Perturb-map analysis plan (AMENDMENT-010, 010a) was committed
to git **before** any knockout outcome was computed.

## Priority 1 — the statistical reliability problem

**Why tests fail in brain** (`step29`, `step29b`, `step34`; 150–200 null draws per setting).

| Candidate cause | Test | Type-I at 20 lesions |
|---|---|---|
| Baseline (as in v3) | — | 0.345 |
| Pseudo-knockout in the control fit | refit without it | 0.355 |
| Leave-out relabelling | ~1% of bins relabelled | — |
| Circular labels | circular | 0.350 |
| Lesions straddling niches | drop minority-niche cells (2%) | 0.333 |
| Straddling (alternative) | residual vs majority niche | 0.360 |
| Spatial baselines | local-control residuals | 0.17–0.215 |
| Unequal variance | Welch studentised | 0.49–0.55 |
| Both | local residuals + Welch | 0.089–0.098 |
| Small-cluster test | CR2 two-niche | 0.15–0.23 (all L) |

Diagnosis: lesion means are centred in every niche (|t| ≤ 0.4), but in brain:
- their variance is 7.5 × the sampling variance (tumour ~1);
- their spread differs 2.6-fold between niches (tumour 1.1-fold);
- they stay correlated up to ~500 µm.

So lesion means are **not exchangeable** between niches. With few lesions per niche, no
permutation, Welch or CR2 test can fix this. This is a property of the tissue, not a bug.

**Remedy 1, matched-null calibration** (`spatialpert.calibration`, `step31`, `step32`).
- **How it works:** the knockout's p-value is ranked among control pseudo-knockouts with the
  same lesion sizes in the same tissue, with ties broken at random.
- **Assumption (stated, not testable):** knockout lesions are placed like random control positions.
- **Error rate:** with cross-fitting, lesion permutation holds 0.043–0.065 in both tissues.
- **Power:** brain power stays ≤ 0.11, and the paper says plainly that niche effects of this
  size are undetectable there.
- **Diagnostic:** a lesion-exchangeability check flags such tissues before testing. Variance
  ratio 1.18 in tumour versus 6.9 in brain.
- **A bug caught along the way:** a fixed-width string truncated the pseudo-knockout label, so
  every pseudo-knockout was untestable. Fixed, and the unit test now requires real outputs.

**Remedy 2, effect estimation, not only tests** (`smallcluster.niche_contrast`, `step33`).
- Ψ = Δ_b − Δ_a with CR2 lesion-clustered intervals (Imbens-Kolesár df).
- Bias ≤ 0.012.
- Coverage: 0.95–0.96 in tumour and 0.92–0.95 in brain. Naive cell-level intervals cover
  only 0.46–0.64 in brain.

**Not done:** Kenward-Roger / Satterthwaite mixed models and dream. They need R, and the R
installer needs Windows administrator approval, which was cancelled at the prompt. The
authors can approve it and these can be added.

## Priority 2 — a new spatial-context model

Learned context encoder (`spatialpert.context_gnn`, optional `torch` extra, CPU).
- **Design:** one-hop, two-ring aggregation of control neighbours only. The cell itself is
  always excluded, and pseudo-knockout cells are removed before embedding.
- **A leak caught:** my first two-hop version let a cell's own features return through its
  neighbours, so I replaced it.
- **Benchmark design (`step35`, three seeds):** the same downstream test for every
  representation. The injected effect follows the local host-cell fraction, so no partition
  is favoured by construction.

| Representation | Max type-I | Power, 10 lesions | Power, 20 lesions |
|---|---|---|---|
| Learned context | 0.05 | 0.11 | 0.17 |
| Fixed kernel (v3) | 0.08 | 0.06 | 0.08 |
| Coordinates | 0.05 | 0.08 | 0.10 |
| Own expression | 0.20 ✗ | ~0.02 | ~0.02 |

The learned encoder roughly doubles power while staying calibrated. It is the only
representation never above 0.05.

**Honest limits:**
- The gain is modest, and absolute power stays low.
- It is shown for one effect model in one tissue.
- Histology input, uncertainty in niche assignment, and stability across sections or
  animals were **not** done.

## Priority 3 — stronger benchmark

| Scenario requested | Status |
|---|---|
| No effect / same effect everywhere / niche-specific effect | Done (v2–v3) |
| Different baseline by niche | Done (niche means in the residual model) |
| Spatially correlated residuals | Done implicitly by the real brain and tumour tissue (step29b); no synthetic spatial-noise model |
| Unequal lesion sizes, missing niches | Done (step26 realism; untestable designs reported) |
| Lesions nested in sections and animals | Done on real data: Perturb-map within-mouse permutation, mouse-clustered CR2 |
| Composition changes | Partly: the effect tied to host-cell fraction (step35); no composition-only effect scenario |
| Bias, RMSE, coverage | Done (step33) |
| FDR | Not done as a separate simulation (BH used on real data) |
| Robustness to niche misclassification | Partly (circular vs leave-out labels; four callers) |

## Priority 4 — fair comparisons

- The representation stage was compared with the downstream estimator held fixed (`step35`).
- The full pipelines were compared across eight tests (v3).
- SpaGCN, BANKSY and CellCharter niches were tested (v3; SpaGCN subsampled, stated).

## Priority 5 — independent biological data

**Perturb-map imaging** (BioImage Archive S-BIAD267), Lymphoid panel:
- 41 lobes from 8 mice; 1,377 curated lesions; 12 eligible knockouts.
- **Integrity:** the raw tables of 9 lobes were kept with SHA-256 checksums, and the extracted
  tables equal them exactly (9/9). The other 32 raw tables were deleted after extraction;
  this is stated in Methods.
- **Pre-registered result (negative):** no knockout shows niche-dependent CD8, CD4 or B220
  infiltration (minimum q = 0.57 over 36 tests).
- **The paper's argument on real data:** Tgfbr2's within-mouse permutation p = 0.022 becomes
  calibrated p = 0.060, so an uncalibrated analysis would have reported it.
- **Sanity check (exploratory):** B2m −0.32 (p = 0.0003), Tgfbr2 −0.46 (p = 0.0013) and Socs1
  +0.46 (p = 0.0029) change CD8 infiltration consistently across 8 mice. Tgfbr2 T-cell
  exclusion and Socs1 T-cell infiltration both reproduce the original study (Dhainaut et al.,
  Cell 2022).

**Honest limits:**
- No non-targeting controls exist; knockouts are compared with each other.
- The readout is a protein panel, not a transcriptome.
- Only 11% of lesions sit in immune-rich surroundings, so only large niche differences are
  detectable.
- Pixel size is assumed (0.5 µm).
- No **positive** niche-dependent biological finding was obtained. Perturb-FISH and Spatial
  Perturb-seq data were not obtained; the request e-mails are drafted but still need to be
  sent by the authors.

## Priority 6 — usable tool

New public functions, all unit-tested:
- `calibrate_p`, `matched_null_p`, `lesion_exchangeability`
- `niche_contrast`
- `ContextEncoder`

Still author-only: licence, GitHub, Zenodo and PyPI release (`submission/CHECKLIST.md`).

## Section 6 — Nature Methods readiness, honestly

| Criterion | Status |
|---|---|
| Calibration problem resolved or delimited | **Yes**: explained, calibrated remedy, diagnostic, stated assumption |
| Demonstrably novel method | **Partly**: calibration + Ψ intervals + context encoder; the encoder gain is modest |
| Multiple independent datasets | **Partly**: tumour Xenium, brain Xenium, Perturb-map imaging (8 mice) |
| Fair comparisons | Yes |
| Important biological information | **No positive niche finding**: an honest negative, with replication of known main effects |

Our assessment: the current version is a solid methods and benchmark paper with a real-data application. It
is competitive for Nature Methods as an Analysis only if the editors value the calibration
result and the negative Perturb-map finding. The safer route is Genome Biology or Nature
Communications. A positive niche-dependent finding would most likely need a screen designed
with more immune-rich lesions, or the Perturb-FISH data.
