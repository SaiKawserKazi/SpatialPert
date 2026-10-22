# Items only the authors can complete

Nothing here has been sent or published. Each item is ready to use.

## 1. Data requests (review M3; Stage 3)

**To:** Loïc Binan and Samouil L. Farhi (Perturb-FISH, *Cell* 2025)
**Subject:** Perturb-FISH tumour xenograft cell tables

> Dear Dr Binan and Dr Farhi,
>
> We are developing SpatialPert, a clone-aware test of whether knockout effects depend on
> tissue location, and would like to apply it to your A375 xenograft Perturb-FISH screen. We
> could not find the tumour cell tables in the SSPsyGene or Brain Image Library deposits.
>
> Could you share the per-cell tables used in `tumoranalysis/clusterCells.ipynb`
> (`merfishgoodCells.csv`, `merfishgoodCellsWithGeneName.csv`, `perturbationIDx.csv`,
> `positionsgoodcells.csv`), and say how many sections and animals they cover? We will cite your
> work and share our code and results with you before submission.
>
> With thanks, [name, affiliation]

**To:** Kok Hao Chen, Shyam Prabhakar and Wei Leong Chew (Spatial Perturb-seq, *Nat. Commun.* 2026)
**Subject:** Spatial Perturb-seq cell-level data

> Dear colleagues,
>
> We would like to test SpatialPert, a clone-aware method for location-dependent knockout
> effects, on your Spatial Perturb-seq mouse brain data. Could you point us to, or share, the
> cell-level expression, coordinates and guide assignments? We will cite your work and share
> results before submission.
>
> With thanks, [name, affiliation]

## 2. Pre-registration timestamp (review §5, last point)

1. Create a project on https://osf.io and choose "Registrations → Open-ended registration".
2. Upload `RESEARCH_PROTOCOL.md` as it was at tag `protocol-v1.0`
   (`git show protocol-v1.0:RESEARCH_PROTOCOL.md > protocol_v1.0.md`).
3. In the description, give the commit hash: `git rev-list -n 1 protocol-v1.0`.
4. The OSF date will post-date the git tag. State both honestly in Methods: "frozen in git on
   2026-10-20 (hash …), deposited on OSF on [date]".

## 3. Code archive with DOI

1. Push the repository to GitHub (a public repository is needed for Zenodo).
2. Enable the repository at https://zenodo.org/account/settings/github/.
3. Create a GitHub release `v0.2.0`; Zenodo then mints a DOI automatically.
4. Put the DOI in `CITATION.cff`, the README and the manuscript Code availability statement.

## 4. PyPI (optional)

`uv build` already produces `dist/spatialpert-0.2.0-py3-none-any.whl`. Upload with
`uv publish` using your own PyPI token. Check that the name `spatialpert` is free first.

## 5. Manuscript fields

Still needed from the authors: author list, affiliations, corresponding author, funding,
competing interests, author contributions, licence choice (for the code), and the
repository/DOI URL.
