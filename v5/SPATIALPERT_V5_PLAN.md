# SpatialPert V5 plan

Each item lists the experiment, the outputs and a go/no-go rule. Compute estimates are for
the authors' laptop (Intel i5-7200U, 2 cores / 4 threads, 8.5 GB RAM, no GPU). Version
control: the author manages git/GitHub; this plan writes files only.

## 0. Status of P0 (audit)

- **P0.1** Evidence inventory: provenance mechanism documented (`SPATIALPERT_V5_AUDIT.md` §1);
  claim table `v5_claims.csv` written after reproduction.
- **P0.2** Discrepancies: 20 logged (D1–D20). D1 resolved; D2, D3 and D6 are high severity.
- **P0.3** Reproduction: steps 27, 28, 31 and 37 reproduced exactly. Steps 33 (tumour), 35,
  36 and step19 (L = 10) are running. Snapshot of the v0.0 results: `v5/baseline_v0.0/`.
- **Go/no-go:** no claim-expanding work is reported until P0.3 completes. If a central
  result fails to reproduce, the failure is documented and the claim withdrawn.

## 1. Validity contract (P1.1), to appear in the manuscript

- **Unit.** The operational *spatial lesion*: a connected component of knockout cells at
  100 µm (Xenium), or a curated polygon (Perturb-map). It is not a verified clone.
- **Estimand.** For knockout k and niches a, b: Ψ_k,g = Δ_k,b,g − Δ_k,a,g. Here Δ_k,n,g is the
  mean residual of knockout cells in niche n relative to niche-matched controls (Xenium) or
  to other knockouts' lesions (Perturb-map). It is an **association** among the lesions that
  were observed. It is a causal perturbation-by-niche effect only if knockout assignment,
  survival and placement are independent of niche-specific baselines.
- **Null.** H0: Ψ_k,g = 0 for all g.
- **Test statistic.** Between-niche sum of squares of lesion-aggregated residuals.
- **Calibration claim.** Under H0, the knockout's statistic is exchangeable with the
  statistics of pseudo-knockouts produced by the null-generating mechanism. This requires:
  - **A1 (placement):** knockout lesions are placed like the pseudo-knockout lesions,
    including their niche occupancy (relaxed to "within niche" in the niche-matched variant).
  - **A2 (residual distribution):** under H0 the knockout changes lesion residuals at most
    by a niche-constant shift. A constant main effect cancels in the between-niche statistic.
  - **A3 (definition):** knockout and pseudo-knockout lesions are defined by the same rule.
- **Conditioning.** Calibration is conditional on the observed lesion sizes (and, in the
  matched variant, on the lesions' niche occupancy). It is not conditional on unobserved
  selection.
- **Not identifiable from control data.** Whether a knockout changes survival or placement
  in a way that correlates with niche-specific baselines (violations of A1/A2) cannot be
  identified. Diagnostics can flag observable departures, such as niche-occupancy or size
  mismatch, but cannot prove A1/A2.

## 2. Null-scenario benchmark (P1.2–P1.4)

Script: `scripts/step40_scenarios.py` (new).
- **Tissues:** tumour (Xenium) and brain (binned Xenium).
- **Lesions:** L ∈ {10, 20}, 400 cells per knockout.
- **Tests:** cell-level permutation, lesion permutation, CR2, and lesion permutation
  calibrated two ways (random-placement pool, niche-matched pool).
- **Pool:** 2,000 pseudo-knockouts per tissue × L from the scenario-A generator, as a user
  would build it.
- **Test draws:** 1,000 per scenario cell (Monte Carlo SE ≈ 0.007 at 0.05); 2,000 for A and
  the B mid-level.

| Scenario | Mechanism (all with Ψ = 0) | Type |
|---|---|---|
| A | Random placement, same generator as the pool | Favourable (sanity) |
| B1–B3 | Lesion seeds enriched in one niche, odds 2, 5, 20 | Plausible violation |
| B4 | Seeds enriched near niche boundaries | Plausible violation |
| C1 | Lesion size depends on niche (×2 in one niche) | Perturbation-dependent growth (sensitivity) |
| C2 | Lesions in one niche lost with probability 0.5 (survival selection) | Sensitivity |
| D | Spatially structured baselines: brain vs tumour (real tissue) | Empirically grounded |
| E | Heavy-tailed lesion sizes (real Xenium lesion-size pool) | Empirically grounded |
| F | Hierarchical: pseudo-knockouts from Perturb-map lesions with imbalanced mice; within-mouse permutation and calibration | Real-data null |
| G1–G2 | Lesion-definition error: seeds 150 µm apart so 100 µm linkage merges lesions; linkage at 50 µm splits lesions | Adversarial |
| H1 | Knockout main effect, no interaction (all niches) | Null with an effect |
| H2 | Lesion-specific random effects (τ = 1 × residual SD) | Plausible |

- **Reported per cell:** type-I at α = 0.05 with Wilson CI; testable (non-abstention) rate;
  number of draws and seeds.
- **Go/no-go:**
  - The calibrated method claims validity in a scenario only if type-I ≤ 0.06 there with
    the upper CI bound < 0.075.
  - Where it fails (expected in B and C), the claim is narrowed, and the niche-matched
    variant and/or diagnostics with abstention are tested.
  - The paper will state the failure scenarios explicitly.

**Compute:** ≈ 30,000 draws × ~0.4 s ≈ 3.5 h on 2 cores.

## 3. Method changes to test (P2)

- **P2.1 Niche-matched pseudo-knockouts.** Each pseudo-lesion is seeded in the same niche as
  the corresponding knockout lesion, conditioning on the observed niche occupancy (not on
  outcomes). Compared with the random pool in every scenario. Ablation: size matching
  on/off. Report pool support (number of testable draws).
- **P2.2 Diagnostics and abstention.** Four diagnostics:
  1. niche-occupancy mismatch between the knockout and the pool (total-variation distance);
  2. lesion-size mismatch (KS);
  3. lesion exchangeability (variance ratio, niche SD ratio);
  4. pool support (fraction testable).

  The abstention rule returns "insufficient support" when support < 50%, and "assumptions
  not assessed / mismatch" when occupancy TV > t, with t tuned on development seeds and
  evaluated on held-out seeds. Report false-alarm (scenario A) and detection (B, C) rates.
- **P2.3 Alignment.** Ψ̂ and its CR2 interval are computed on the same lesion set and niche
  pair as the test. The paper explains that the interval is model-based (CR2) and the
  p-value is permutation-calibrated. Coverage is already measured (step33); it is extended
  to scenarios B and H.
- **P2.4 Encoder.** Spatial split: the encoder is trained on control cells in the left half
  of the section, and pseudo-knockouts are drawn from the right half; 3 seeds.
  Compared with fixed kernel, coordinates and no encoder. The encoder leaves the headline
  unless the gain holds out of sample.

## 4. Power, usable-result probability and MDE (P3)

- **Grid (tumour):**
  - L ∈ {5, 10, 20, 40}; λ ∈ {0.25, 0.5, 1, 2};
  - niche balance: as observed, and forced imbalance 90:10;
  - 300 effect draws per cell, calibrated test with abstention.
- **Brain:** L ∈ {10, 20, 40}, λ ∈ {0.5, 2}.
- **Reported:**
  - P(returned);
  - power | returned;
  - P(returned and significant);
  - MDE at 80% overall power where estimable.
- **Compute:** ≈ 2 h.

## 5. Real-data re-audit (P4)

- **Xenium:** report all four knockouts with testability reasons (no new data).
- **Perturb-map:**
  - (i) attrition table (D1);
  - (ii) fix the step25 code count (D3);
  - (iii) declare the ring filter and run a sensitivity without it (D2);
  - (iv) BH across all 36 tests as a sensitivity (D5);
  - (v) a pool excluding knockouts with significant main effects (D6);
  - (vi) the scenario-F null on the real lesion table.
- **Independent dataset (P4.2):** no further public spatial-perturbation dataset with
  lesion-level replication is available to us. Perturb-FISH and Spatial Perturb-seq data were
  not obtained. Stated explicitly; no "independent validation" claim.
- **Positive control:** none independently established for niche dependence. The
  Tgfbr2/Socs1 main effects are reproductions of published main effects, not niche
  controls.

## 6. Manuscript v5.0 (after 1–5)

- **File:** separate `reports/SpatialPert_manuscript_v5.0.docx` / `.pdf`; v0.0 kept until the
  author decides.
- **Title without "clonal":** working title "Lesion-level replication and calibration for
  niche-dependent perturbation effects in spatial CRISPR screens".
- **Structure (brief §10):**
  1. problem;
  2. formulation with the validity contract;
  3. method (lesions, pseudo-knockouts, calibration, diagnostics, abstention);
  4. benchmark with favourable and adversarial scenarios;
  5. power and design;
  6. real data;
  7. discussion.
- **Target journal:** PLOS Computational Biology or Bioinformatics Advances, decided by the
  results; Nature Methods only if the criteria in §12 of the brief are met.

## 7. Deliverables checklist

- [x] `SPATIALPERT_V5_AUDIT.md`
- [x] `SPATIALPERT_V5_PLAN.md`
- [ ] `v5_experiment_manifest.csv` (every run: script, args, seed, output, runtime)
- [ ] `v5_claims.csv` (claim → result → figure → script → limitation)
- [ ] `V5_REPRODUCIBILITY_CHECKLIST.md`
- [ ] Manuscript v5.0

## 8. Author-only items (cannot be done here)

Author list and affiliations; repository URL; licence; release tag/DOI; Zenodo deposit;
sending the data-request e-mails; journal choice.
