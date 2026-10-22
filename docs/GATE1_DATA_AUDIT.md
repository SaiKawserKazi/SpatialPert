# Gate 1: data audit of Perturb-map Visium (GSE193460)

**Verdict: FAIL.** Only 2 of the 35 library knockouts can be tested (Tgfbr2, Ifngr2). The
pre-set rule needs 10. Perturb-map Visium cannot be the discovery dataset.

Run: `uv run python scripts/run_gate1.py` → `reports/gate1_data_audit.json` (2026-10-20).

## Rule (fixed in the proposal before data were opened)

At least 10 knockouts with at least 30 perturbed spots in at least 2 spatial domains. Domains
were not called, so the gate used an upper bound that no domain calling can beat: a knockout
needs at least 60 labelled spots in total. A stricter count (30 spots in each of 2 sections)
is reported beside it.

## Files are consistent

All 4 sections were read: 32,289 features each, and every annotated barcode is in the matrix and
on tissue. KP_4 has 2,546 on-tissue spots but 1,355 annotated ones, because the authors removed
spots they marked as off-tissue. File sizes match GEO's `filelist.txt` exactly.

## Knockout labels exist for 3 genes only

| Knockout | KP_1 | KP_2 | KP_3 | KP_4 | Total spots | Testable |
|---|---|---|---|---|---|---|
| Tgfbr2 | 75 | 62 | - | 87 (2 lesions) | 224 | yes (3 sections) |
| Ifngr2 | - | 41 | 64 | - | 105 | yes (2 sections) |
| Jak2 | 21 | - | - | - | 21 | no |
| other 32 genes | - | - | - | - | 0 | no |

Most tumour spots carry no knockout identity. Of 2,616 k-means tumour spots, 350 are knockout
labelled, 1,061 are in unassigned `KP_n-m` lesions and 595 are `periphery`.

## Why the labels cannot be recovered from the transcripts

- The authors' notebook ([srose89/PERTURB-map](https://github.com/srose89/PERTURB-map),
  `pmap_visium_processing_20211127.Rmd`) shows how the labels were made. Tumour spots are called
  by WPRE UMIs (4 or more). Leiden clusters were then **named by hand** from imaging mass
  cytometry on a serial section. There is no per-spot decoding.
- The matrix has barcode-related features, but they are almost empty. Total UMIs per section:
  `Procode` 0-2, `mCherry_Nter` 2-22, `linker` 0. WPRE has 6,374-8,591, but WPRE is shared by
  every construct, so it identifies tumour, not knockout.
- The paper ([PMC8992964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8992964/)) says few lesions
  per perturbation were sampled by Visium. Full imaging data are "available on request". The
  segmented imaging is in BioImage Archive **S-BIAD267**, but it is protein imaging, not
  transcriptome.

## What remains usable

- **Positive control.** Tgfbr2 has 224 spots in 3 sections from 2 mice, and Ifngr2 has 105 spots
  in 2 sections. This is enough for the H1 positive-control check the proposal requires.
- **Domain calling and composition.** All 6,363 annotated spots support domain calling and
  deconvolution.

## Pivot options (decision for the user)

1. **Recommended.** Make Perturb-FISH (Binan et al., *Cell* 2025) the discovery set. It has
   35 NF-kB genes, single cells in vivo, and a 500-gene panel, which gives far more cells per
   knockout. Keep Perturb-map for the Tgfbr2/Ifngr2 positive control. Run the same gate on
   Perturb-FISH next.
2. Wait for Perturb-DBiT accessions (pending in the preprint), then re-gate.
3. Narrow the study to the 2 testable knockouts. This is a case study, not a general claim.
4. Request the Pro-Code imaging from the authors and register it to the Visium sections.
   This is slow and depends on them.
