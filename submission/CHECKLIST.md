# Submission checklist (current manuscript v0.0)

Done in the repository:

- [x] Manuscript v0.0, Analysis format (`reports/SpatialPert_manuscript_v0.0.docx` / `.pdf`)
- [x] Six figures + six Extended Data figures, vector PDF, 183 mm (`reports/figures/v3/`, `reports/figures/v4/`)
- [x] Source data (`reports/figures/v3/source_data/`, `reports/figures/v4/source_data/`)
- [x] Supplementary Note 1 (pre-registration history), inside the manuscript
- [x] Point-by-point response to the latest review (`review/RESPONSE_REVIEW3.md`)
- [x] Presubmission enquiry, cover letter, Reporting Summary answers, referee candidates (this folder)
- [x] Independent recheck 489/489 of the current version's numbers (`reports/recheck.json`)
- [x] Data-request e-mails drafted (`docs/SUBMISSION_PACKAGE.md` §1)

Only the authors can do:

- [ ] Author list, affiliations, ORCIDs, corresponding author; CRediT contributions; acknowledgements
- [ ] Choose a licence (MIT or BSD-3-Clause) and add `LICENSE`
- [ ] Push the repository to GitHub; create a release; archive it on Zenodo (DOI)
- [ ] Deposit benchmark per-draw outputs (`reports/step19_*.csv`, `step22`, `step26`, labels in `data/interim/a6/`) on Zenodo
- [ ] `uv build` and upload to PyPI (name `spatialpert`, check availability)
- [ ] Send the Perturb-FISH and Spatial Perturb-seq data requests; record replies
- [ ] Vet referee candidates (`SUGGESTED_REFEREES.md`)
- [ ] Post the preprint on bioRxiv; replace placeholders ([GitHub URL], [DOI]) in manuscript and letters
- [ ] Send the presubmission enquiry; then the full submission

Open science items that would strengthen the full submission (see the response, "Not done"):

- [ ] Perturb-map lesion-level IFN-pathway test (needs the 15.7 GB segmentation tables and the Pro-Code map)
- [ ] R-based small-cluster methods (lme4/lmerTest with Satterthwaite or Kenward-Roger df; dream)
- [ ] A third tissue with segmented cells
