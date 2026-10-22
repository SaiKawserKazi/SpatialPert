"""Gate 3: the same week-1 rule on Perturb-DBiT (GSE319277), per screen.

A section counts only if it has BOTH an sgRNA table and an expression matrix. A pixel counts
for a knockout if it carries >= 1 UMI of any guide against that gene. Writes
reports/gate3_dbit_audit.json.
"""

import json
from pathlib import Path

import pandas as pd

from spatialpert.dbit import pixel_knockouts, read_sg_counts
from spatialpert.gate import evaluate_gate
from spatialpert.xenium import guide_cell_table, negative_probe_rate, read_xenium_h5

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE319277"
XENIUM = ROOT / "data" / "raw" / "GSE319123" / "GSM9511110_Xenium_lung_cell_feature_matrix.h5"
OUT = ROOT / "reports" / "gate3_dbit_audit.json"

# section -> (sgRNA file, expression file or None). Bru.LC.6/8 have two sg libraries with the
# same guides; the dedicated gRNA library is used.
SECTIONS: dict[str, tuple[str, str | None]] = {
    "Bru.LC.1": ("GSM9514747_sgCounts_Bru.LC.1.txt.gz", None),
    "Bru.LC.2": ("GSM9514748_sgCounts_Bru.LC.2.txt.gz", None),
    "Bru.LC.3": ("GSM9514749_sgCounts_Bru.LC.3.txt.gz", None),
    "Bru.LC.6": ("GSM9514752_sgCounts_Bru.LC.6_grnalib.txt.gz", "GSM9514752_Bru.LC.6_allRNA.tsv.gz"),
    "Bru.LC.7": ("GSM9514755_sgCounts_Bru.LC.7.txt.gz", None),
    "Bru.LC.8": ("GSM9514756_sgCounts_Bru.LC.8_grnalib.txt.gz", "GSM9514756_Bru.LC.8.tsv.gz"),
    "Bru.LC.1.33.3": ("GSM9514758_sgCounts_Bru.LC.1.33.3.txt.gz", "GSM9514758_Bru.LC.1.33.3.tsv.gz"),
    "mTSG.LV.2": ("GSM9514753_sgCounts_mTSG.LV.2.txt.gz", "GSM9514753_mTSG.LV.2.tsv.gz"),
    "mTSG.LV.3": ("GSM9514754_sgCounts_mTSG.LV.3.txt.gz", "GSM9514754_mTSG.LV.3.updated.tsv.gz"),
}
SCREENS = {"Brunello lung (HT-29, NSG)": "Bru.", "mTSG liver (CRISPR-GEMM)": "mTSG."}


def ko_table(sections: list[str]) -> pd.DataFrame:
    rows = []
    for s in sections:
        ko = pixel_knockouts(read_sg_counts(RAW / SECTIONS[s][0]))
        ko = ko[ko["gene"] != "NTC"]
        for gene, n in ko.groupby("gene")["pixel"].nunique().items():
            rows.append({"section": s, "knockout": gene, "n_spots": int(n), "n_lesions": 1})
    return pd.DataFrame(rows, columns=["section", "knockout", "n_spots", "n_lesions"])


def main() -> None:
    report: dict[str, object] = {"rule": "pass if >= 10 knockouts have >= 30 pixels in >= 2 domains"}
    for screen, prefix in SCREENS.items():
        secs = [s for s in SECTIONS if s.startswith(prefix)]
        paired = [s for s in secs if SECTIONS[s][1] is not None]
        entry = {}
        for label, use in (("paired_sections", paired), ("all_sg_sections", secs)):
            table = ko_table(use)
            res = evaluate_gate(table, library_size=0)
            entry[label] = {
                "sections": use,
                "testable_upper_bound": list(res.testable_upper_bound),
                "testable_replicated": list(res.testable_replicated),
                "passed_upper_bound": res.passed,
                "passed_replicated": len(res.testable_replicated) >= 10,
            }
        report[screen] = entry

    x = read_xenium_h5(XENIUM)
    one = guide_cell_table(x)
    two = guide_cell_table(x, min_umis=2).set_index("probe")["cells"]
    one["cells_2umi"] = one["probe"].map(two)
    above = one[one["cells"] >= 60]
    report["Xenium HT-29 lung (GSE319123, single cell)"] = {
        "cells": x.matrix.counts.shape[0],
        "negative_probe_rate_per_cell": negative_probe_rate(x),
        "probes": one.to_dict(orient="records"),
        "testable_upper_bound": above["target"].tolist(),
        "passed_upper_bound": len(above) >= 10,
        "min_enrichment_among_testable": float(above["enrichment"].min()),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
