# Presubmission enquiry — Nature Methods (Analysis)

Submit through the Nature Methods online system as "Presubmission enquiry". Attach the
abstract below and Fig. 1 + Fig. 3 (`reports/figures/v3/`). Fill in the bracketed fields.

---

Dear Editor,

We would like to ask whether Nature Methods would consider an Analysis entitled
**"Clonal structure confounds tests of location-dependent knockout effects in spatial
CRISPR screens"**.

Spatial CRISPR screens (Perturb-map, Perturb-FISH, Perturb-DBiT and spatial Perturb-seq,
2022-2026) let researchers ask whether a knockout acts differently in different tissue
niches. Knockout cells, however, grow as clonal lesions, so the cells being compared are
not independent. No study has benchmarked how tests of niche dependence behave under this
structure, and published analyses use cell-level tests.

Using semi-synthetic knockouts built from public Xenium tumour and brain sections, we compared
eight tests, including small-cluster corrections from econometrics, under a pre-registered
protocol, and applied the methods to Perturb-map imaging of 1,377 lesions from eight mice.
Our main findings are:

1. Cell-level tests detect niche effects in up to 62% of null draws, and letting knockout
   cells contribute to their own niche definition inflates even lesion-level tests — a
   pitfall for any simulation of spatial perturbations.
2. Lesion permutation and CR2 are calibrated in tumour, but every test fails in brain, where
   lesion baselines are large and unequal across niches. Calibrating against matched control
   pseudo-knockouts restores error control in both tissues (type-I 0.043-0.065), and
   lesion-clustered intervals cover the knockout-by-niche contrast (92-96%). A learned
   context encoder roughly doubles power (8% to 17% with 20 lesions); power still needs about
   20 lesions per knockout.
3. In Perturb-map, B2m, Tgfbr2 and Socs1 change T-cell infiltration consistently across mice
   (reproducing the original study for Tgfbr2 and Socs1), but no knockout's effect differs
   detectably between immune-rich and immune-poor niches. Tgfbr2 shows how an uncalibrated
   test misleads: p = 0.022 before and 0.060 after calibration.

We provide SpatialPert, an open-source Python toolkit that applies the calibrated test,
reports testability per knockout and plans screens from the simulated power curves. We
believe the design rule — count lesions, not cells, and aim for about 20 per knockout — is
directly useful to the many groups now building spatial screens.

The main text has about 2,900 words, six figures and six Extended Data figures. A preprint
will be posted on bioRxiv [date]. We are happy to send the full manuscript.

Yours sincerely,
[Corresponding author, affiliation, e-mail, ORCID]

---

**Abstract (139 words)** — copy from the current manuscript (`reports/SpatialPert_manuscript_v0.0.pdf`), "Abstract".
