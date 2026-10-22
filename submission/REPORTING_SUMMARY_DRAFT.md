# Nature Portfolio Reporting Summary — draft answers

Copy into the official form (download from the Nature Portfolio website at submission).
Answers describe this computational study; sections on animals, human participants,
antibodies, cell lines, ChIP-seq, flow cytometry and MRI are "not applicable" because no new
experiments were done.

## Statistics
- Sample size: the unit of analysis in benchmarks is a simulated draw. Null settings: 300
  draws per lesion number (stage 1b), 200 (niche callers, count-level realism); effect
  settings: 150 (stage 1b), 100 (callers, realism). At 300 draws a Wilson 95% interval
  around a type-I error of 0.05 has half-width about 0.025. Mixed model: 40 null /
  20 effect draws (registered deviation; runtime).
- Tests: permutation tests (250 permutations in benchmarks, 1,000 in real data), CR2 t-tests
  with Imbens-Kolesár df, wild cluster bootstrap (199-999 draws), sign-flip test (exact or
  2,000 flips), MixedLM Wald tests. All two-sided where applicable.
- Multiple comparisons: Benjamini-Hochberg across genes (CR2, mixed model) and knockouts.
- Uncertainty: 95% Wilson intervals for rejection rates; IQR for ICC.
- Exclusions: cells with < 20 counts; cells with two or more guides; cells without a
  control neighbour within 100 µm (no niche label). Pre-specified in the protocol.
- Replication: the semi-synthetic benchmark was repeated in a second tissue and with four
  niche callers; an independent recheck script recomputes all reported rates (489/489).
- Randomisation/blinding: not applicable (no experiments); analysis plan pre-registered
  before outcomes were examined (git tag protocol-v1.0; amendments dated).

## Software and code
- Data collection: none (public data).
- Data analysis: Python 3.12; numpy, scipy, pandas, scikit-learn, statsmodels, h5py,
  matplotlib; pybanksy 1.3.5, SpaGCN 1.2.7, CellCharter 0.3.8; SpatialPert 0.2.0
  ([GitHub URL], Zenodo [DOI]).

## Data
- GEO GSE319123, GSE319277, GSE193460; BioImage Archive S-BIAD267; Zenodo 11209134
  (CC BY 4.0). Source data provided with figures; benchmark outputs at Zenodo [DOI].
- No restrictions on data availability.

## Field-specific: life sciences
- Study description: computational benchmark of statistical tests on public spatial
  transcriptomics data, with semi-synthetic perturbations.
- Research sample: one Xenium section of HT-29 lung metastases in mouse (Perturb-DBiT study),
  Xenium mouse brain subset, Perturb-map lesion annotations (8 mice).
- Sampling strategy, data collection, timing: not applicable (public data).
