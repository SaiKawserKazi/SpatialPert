# Gate 3: Perturb-DBiT (GSE319277) and its Xenium validation section (GSE319123)

**Verdict: PASS on the upper bound for the Xenium single-cell section (12 knockouts) and the
paired DBiT lung sections (13). Location-only replication across 7 DBiT lung sections: 17
knockouts.** These replace Perturb-FISH, whose tumour tables are not public (Gate 2).

Run: `uv run python scripts/run_gate3.py` → `reports/gate3_dbit_audit.json` (2026-10-20).
Only guide counts were examined. No expression outcome was looked at before the protocol freeze.

## Data

| Accession | Content | Used for |
|---|---|---|
| GSE319123 | Xenium, HT-29 lung metastasis (NSG mouse), 352,913 cells, 5,001-gene panel + 27 guide probes | Discovery (single cell) |
| GSE319277 | DBiT-seq, 15 sections. Brunello HT-29 lung: 7 with guides, 3 also with expression (LC.6, LC.8, LC.1.33.3) | Cross-platform replication; clone location |
| GSE193460 | Perturb-map Visium | Tgfbr2/Ifngr2 positive control |

LC.4 and LC.5 have expression but no guide data, so they are excluded. All files pass
`gzip -t`, and their sizes match GEO.

## Xenium guide detection

Negative-control probes are positive in 1.08e-5 of cells, which is about 3.8 chance cells per
probe. Knockouts with ≥ 60 positive cells (≥ 1 UMI):

| Target | Cells | ≥ 2 UMI | Above chance |
|---|---|---|---|
| DIO3 | 1,781 | 1,232 | 469× |
| PRSS23 | 783 | 432 | 206× |
| ESAM | 617 | 409 | 162× |
| CIITA | 341 | 144 | 90× |
| PTH1R | 281 | 190 | 74× |
| GRIN2C | 240 | 161 | 63× |
| ALDH5A1 | 112 | 32 | 29× |
| SUMF2 | 100 | 58 | 26× |
| TXNDC16 | 86 | 29 | 23× |
| EHF | 74 | 47 | 19× |
| VAMP5 | 72 | 42 | 19× |
| M6PR | 67 | 28 | 18× |

TEX13A and SLC4A2 probes detected 0 cells. TEX13A is common in DBiT LC.1.33.3, so its Xenium
probe probably failed. It is excluded.

## Caveats carried into the protocol

- Xenium is one section. Replication comes from DBiT sections on a second platform, for
  overlapping knockouts (PRSS23, SUMF2, C10orf53, GPR84, CIITA), and from within-section
  spatial blocks.
- The 5K panel is human; host mouse lung cells cross-hybridise. Tumour cells must be called
  before any effect estimate (Step 4).
- DBiT pixels (20–50 µm) mix clones and cell types, so DBiT effects are pixel-level.
- The full gate (≥ 30 cells in each of ≥ 2 domains) is applied after domain calling, as the
  first check of Step 4.
