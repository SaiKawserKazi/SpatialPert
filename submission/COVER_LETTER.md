# Cover letter — full submission (Nature Methods, Analysis)

[Date]

Dear [Editor name],

Please consider our manuscript **"Clonal structure confounds tests of location-dependent
knockout effects in spatial CRISPR screens"** for publication as an Analysis in Nature
Methods. [If applicable: This follows our presubmission enquiry of (date), reference (no.).]

**The problem.** Several platforms now read CRISPR perturbations in intact tissue
(Perturb-map 2022; Perturb-FISH 2025; Perturb-DBiT and spatial Perturb-seq 2026). Their
most distinctive promise is to find knockouts whose effect depends on the microenvironment.
But knockout cells grow as clones, so a knockout is observed in a few lesions rather than in
many independent cells, and lesions and niches are confounded. Whether existing tests can
answer the question at all had not been examined.

**What we did.** Under a pre-registered protocol with dated amendments, we benchmarked eight
tests on semi-synthetic knockouts built from a public Xenium section, with 1-20 lesions per
knockout, leave-out niche labels, size-adjusted power with confidence intervals, a second
tissue, four niche callers, count-level effects with real lesion-size distributions, and
real knockouts.

**Main findings.**
- Cell-level tests are badly miscalibrated (type-I error up to 0.62), and simulations in
  which perturbed cells define their own niches overstate the reliability of lesion-level
  tests.
- In tumour, lesion permutation, CR2 and a within-lesion sign-flip test are calibrated; the
  wild cluster bootstrap and mixed models are not. Power needs about 20 lesions per knockout and
  is lower for realistic effects; a standard design-effect calculation does not predict it.
- Every test fails in a second tissue (brain), where lesion baselines are large and unequal
  across niches. Calibrating against matched control pseudo-knockouts restores error control,
  lesion-clustered intervals estimate the knockout-by-niche contrast, and a learned context
  encoder roughly doubles power.
- In Perturb-map imaging of 1,377 lesions from eight mice, knockouts change T-cell
  infiltration (reproducing known Tgfbr2 and Socs1 effects), but none does so detectably
  differently between niches; an uncalibrated test would have reported Tgfbr2.

**Why Nature Methods.** The paper gives a practical design rule and an open toolkit
(SpatialPert: calibrated test, testability report, screen planner) at the moment many
groups are designing spatial screens. It fits the Analysis format: a systematic comparison
of established methods with clear recommendations for users.

**Availability.** Code: [GitHub URL], [licence], archived at Zenodo [DOI], installable from
PyPI. All data are public (GEO GSE319123, GSE319277, GSE193460; BioImage Archive S-BIAD267;
Zenodo 11209134); benchmark outputs are deposited at Zenodo [DOI]. A preprint is posted at
bioRxiv [DOI].

**Suggested referees** are listed in the submission form (see `SUGGESTED_REFEREES.md`).
[Excluded referees, if any, with reason.]

The manuscript is not under consideration elsewhere. All authors approved the submission.
The authors declare no competing interests.

Sincerely,
[Corresponding author]
[Affiliation · e-mail · ORCID]
