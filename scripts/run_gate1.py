"""Week-1 data gate on Perturb-map (GSE193460). Writes reports/gate1_data_audit.json."""

import json
from pathlib import Path

import numpy as np

from spatialpert.gate import (
    evaluate_gate,
    knockout_spot_table,
    non_ensembl_features,
    spot_category_counts,
)
from spatialpert.visium import read_10x_h5, read_positions, read_spot_annotation

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE193460"
OUT = ROOT / "reports" / "gate1_data_audit.json"
SECTIONS = {"KP_1": "GSM5808054", "KP_2": "GSM5808055", "KP_3": "GSM5808056", "KP_4": "GSM5808057"}
LIBRARY_SIZE = 35  # Pro-Code/CRISPR library, GEO series design


def main() -> None:
    annotations = {}
    integrity = {}
    transgenes: dict[str, dict[str, int]] = {}
    for section, gsm in SECTIONS.items():
        stem = f"{gsm}_{section}_"
        m = read_10x_h5(RAW / f"{stem}filtered_feature_bc_matrix.h5")
        pos = read_positions(RAW / f"{stem}tissue_positions_list.csv.gz")
        ann = read_spot_annotation(RAW / f"{stem}spot_annotation.csv.gz")
        annotations[section] = ann

        in_tissue = set(pos.loc[pos["in_tissue"] == 1, "barcode"])
        bc = set(m.barcodes)
        integrity[section] = {
            "matrix_spots": m.counts.shape[0],
            "matrix_features": m.counts.shape[1],
            "in_tissue_spots": len(in_tissue),
            "annotated_spots": len(ann),
            "annotated_in_matrix": int(ann["barcode"].isin(bc).sum()),
            "annotated_in_tissue": int(ann["barcode"].isin(in_tissue).sum()),
        }

        extra = non_ensembl_features(m.feature_ids, m.feature_names)
        col = {n: i for i, n in enumerate(m.feature_names)}
        transgenes[section] = {
            n: int(np.asarray(m.counts[:, col[n]].sum())) for n in extra
        }

    table = knockout_spot_table(annotations)
    categories = spot_category_counts(annotations)
    result = evaluate_gate(table, library_size=LIBRARY_SIZE)

    report = {
        "dataset": "GSE193460 (Perturb-map, Dhainaut et al. Cell 2022)",
        "rule": "pass if >= 10 knockouts have >= 30 spots in >= 2 domains; "
        "upper bound used: >= 60 labelled spots in total",
        "integrity": integrity,
        "non_ensembl_features_total_umis": transgenes,
        "spot_categories": categories.to_dict(orient="records"),
        "knockout_spots": table.to_dict(orient="records"),
        "knockout_totals": table.groupby("knockout")["n_spots"].sum().to_dict(),
        "gate": {
            "library_size": result.library_size,
            "labelled_knockouts": list(result.labelled),
            "testable_upper_bound": list(result.testable_upper_bound),
            "testable_replicated": list(result.testable_replicated),
            "min_testable": result.min_testable,
            "passed": result.passed,
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
