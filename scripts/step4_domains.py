"""Step 4 (protocol §2-3): tumour cells, KO calls, domains, covariates, full testability gate.

Writes data/interim/xenium_cells.parquet and reports/step4_domains.json.
"""

import gzip
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from spatialpert.domains import (
    CONTROL,
    MULTI,
    call_tumour,
    choose_domains,
    knockout_calls,
    marker_score,
    neighbour_counts,
    neighbourhood_mean,
    normalise_log1p,
)
from spatialpert.xenium import guide_probe_indices, probe_target, read_xenium_h5

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE319123"
INTERIM = ROOT / "data" / "interim"
REPORT = ROOT / "reports" / "step4_domains.json"
MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDX2", "CEACAM5"]
MIN_COUNTS = 20
EXCLUDED_PROBES = {"TEX13A", "SLC4A2"}  # 0 cells in Gate 3 (failed probes)


def main() -> None:
    x = read_xenium_h5(RAW / "GSM9511110_Xenium_lung_cell_feature_matrix.h5")
    with gzip.open(RAW / "GSM9511110_Xenium_lung_cells.parquet.gz") as fh:
        cells = pd.read_parquet(io.BytesIO(fh.read()))
    m = x.matrix
    assert cells["cell_id"].astype(str).tolist() == m.barcodes, "cell order mismatch"

    gidx = guide_probe_indices(x)
    eidx = [
        i
        for i, (fid, t) in enumerate(zip(m.feature_ids, x.feature_types, strict=True))
        if t == "Gene Expression" and fid.startswith("ENS")
    ]
    genes = [m.feature_names[i] for i in eidx]
    expr = m.counts[:, eidx].tocsr()
    guides = m.counts[:, gidx].tocsr()
    targets = [probe_target(m.feature_names[i]) for i in gidx]

    totals = np.asarray(expr.sum(axis=1)).ravel()
    keep = totals >= MIN_COUNTS
    expr, guides, cells = expr[keep], guides[keep], cells.loc[keep].reset_index(drop=True)
    coords = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)
    x_log = normalise_log1p(expr)

    markers_used = [g for g in MARKERS if g in genes]
    score = marker_score(x_log, [genes.index(g) for g in markers_used])
    tumour = call_tumour(score)

    labels = knockout_calls(guides, targets)
    labels[np.isin(labels, list(EXCLUDED_PROBES))] = MULTI  # never analysed

    ref = tumour & (labels == CONTROL)
    norm_ref = x_log[ref].copy()
    norm_ref.data = np.expm1(norm_ref.data)  # back to normalised (non-log) scale
    mean_ref = np.asarray(norm_ref.mean(axis=0)).ravel()
    analysis = np.where(mean_ref >= 0.1)[0]  # protocol §2.4: mean normalised expression
    pca = PCA(n_components=30, svd_solver="arpack", random_state=0)
    pcs = pca.fit_transform(x_log[ref][:, analysis])

    feats, _ = neighbourhood_mean(coords, pcs, ref, radius=100.0, length=50.0)
    ok = ~np.isnan(feats).any(axis=1)
    km, sil = choose_domains(feats[ref & ok])
    domain = np.full(len(cells), -1)
    assignable = tumour & ok
    domain[assignable] = km.predict(feats[assignable])

    density, host_frac = neighbour_counts(coords, tumour, radius=50.0)

    out = pd.DataFrame(
        {
            "cell_id": cells["cell_id"].astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
            "log_total": np.log(np.asarray(expr.sum(axis=1)).ravel()),
            "tumour_score": score,
            "tumour": tumour,
            "label": labels.astype(str),
            "domain": domain,
            "density50": density,
            "host_frac50": host_frac,
        }
    )
    INTERIM.mkdir(parents=True, exist_ok=True)
    out.to_parquet(INTERIM / "xenium_cells.parquet")
    pd.Series([genes[i] for i in analysis]).to_csv(
        INTERIM / "analysis_genes.csv", index=False, header=False
    )

    tum = out[out["tumour"] & (out["domain"] >= 0)]
    ko = tum[~tum["label"].isin([CONTROL, MULTI])]
    per = ko.groupby(["label", "domain"]).size().unstack(fill_value=0)
    testable = sorted(per.index[(per >= 30).sum(axis=1) >= 2].tolist())
    report = {
        "cells_total": int(keep.size),
        "cells_kept_min_counts": int(keep.sum()),
        "markers_used": markers_used,
        "tumour_cells": int(tumour.sum()),
        "tumour_fraction": float(tumour.mean()),
        "label_counts_tumour": out.loc[out["tumour"], "label"].value_counts().to_dict(),
        "analysis_genes": len(analysis),
        "pca_explained_var_30": float(pca.explained_variance_ratio_.sum()),
        "silhouette_by_k": sil,
        "chosen_k": int(km.n_clusters),
        "tumour_cells_per_domain": tum["domain"].value_counts().sort_index().to_dict(),
        "control_cells_per_domain": tum[tum["label"] == CONTROL]["domain"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "ko_cells_by_domain": {
            k: {int(d): int(v) for d, v in row.items()} for k, row in per.iterrows()
        },
        "testable_knockouts": testable,
        "n_testable": len(testable),
        "gate_branch": (
            "full" if len(testable) >= 10 else "limited" if len(testable) >= 5 else "pooled"
        ),
        "median_total_counts_kept": float(np.median(np.exp(out["log_total"]))),
    }
    REPORT.write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
