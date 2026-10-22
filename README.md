# SpatialPert

**Does tissue location change what a gene knockout does?**

SpatialPert tests whether the transcriptional effect of a knockout differs between spatial
domains of a tissue. It works on spatial CRISPR screens (Xenium/MERFISH single cells,
DBiT-seq/Visium pixels).

Three design choices make the test hard to fool:

1. **Domains come from control cells only.** A knockout cannot define its own context.
2. **Effects are residuals against domain-specific controls,** with covariates such as
   depth, host fraction and local density.
3. **Significance comes from permuting domain labels among clone patches (or tiles), not
   cells.** Spatially clustered clones would otherwise inflate false positives: cell-level
   tests reach a type-I error of 0.57 in our benchmark.

## Install

```
uv sync
```

## Use

```python
from spatialpert.api import location_dependence

res = location_dependence(
    y,            # cells x genes, log-normalised (sparse or dense)
    coords,       # cells x 2, µm
    labels,       # per cell: knockout gene or "control"
    domains,      # per cell: integer domain (< 0 = unassigned)
    covariates=cov,   # optional cells x c
    genes=gene_names,
    unit="patch",     # or "tile"
)
print(res.table)        # knockout, units, cells, domains, T, p, q
print(res.pooled_p)     # pooled test across knockouts
res.effects["DIO3"]     # domains x genes mean residual effect
```

Helpers: `spatialpert.domains` (normalisation, tumour calling, control-only neighbourhood
features, domain selection), `spatialpert.xenium` / `spatialpert.dbit` / `spatialpert.visium`
(readers), `spatialpert.benchmark` (semi-synthetic calibration and power).

## Study

The pre-registered protocol and its amendments are in `RESEARCH_PROTOCOL.md` (tags
`protocol-v1.0`, `protocol-v1.1`).

| Step | Content | Result |
|---|---|---|
| Gates 1-3 | Data audits (Perturb-map, Perturb-FISH, Perturb-DBiT) | `docs/GATE*_*.md` |
| 4 | Tumour cells + control-defined domains (Xenium) | 4 testable knockouts → pooled branch |
| 5 | Controls | PC1 Tgfbr2 r = 0.50 (p = 0.001); PC2 FPR = 0.050 |
| 6 | H1-H5 | H1, H2, H4, H5 not supported; H3 p = 0.040 |
| 7 | Semi-synthetic benchmark | patch test calibrated to n = 200; naive tests up to 0.57 |
| Review 1 | Critique + exploratory revisions | `review/`, `reports/step8_revision1.json` |

## Reproduce

```
uv run python scripts/fetch_gse193460.py
uv run python scripts/run_gate1.py; uv run python scripts/run_gate2.py; uv run python scripts/run_gate3.py
uv run python scripts/step4_domains.py
uv run python scripts/step5_controls.py
uv run python scripts/step6_hypotheses.py
uv run python scripts/step7_benchmark.py
uv run python scripts/step8_revision1.py
uv run python scripts/make_figures.py
uv run python scripts/recheck.py        # independent recomputation, exits 1 on mismatch
powershell -File status.ps1             # ruff + mypy + pytest
```

GSE319277 / GSE319123 files are downloaded with curl from GEO (see `docs/GATE3_DBIT_AUDIT.md`).
# SpatialPert
