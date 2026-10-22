"""V5 audit P0.2: reconcile Perturb-map lesion counts (2,192 polygons vs 1,377 analysed lesions).

Counts, stage by stage, how curated Lymphoid polygons become analysed lesions under the
AMENDMENT-010/010a rules. Writes reports/step38_lesion_attrition.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as Poly
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step36_perturbmap_niche import ANN, D, RING_EXCL, RING_IN, RING_OUT, polygons  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    files = sorted(D.glob("*_Lymphoid.parquet"))
    pairs = pd.concat([pd.read_parquet(f, columns=["procode", "Gene"]).dropna() for f in files])
    code2gene = pairs.groupby("procode").Gene.agg(lambda v: v.value_counts().index[0]).to_dict()
    n = {"annotation_files": len(list(ANN.glob("*_Lymphoid_*curated.fix.json"))), "lobes_with_tables": len(files),
         "polygons_named": 0, "polygons_with_gene": 0, "polygons_ge30_labelled": 0, "polygons_ring_ok": 0}
    unmapped: dict[str, int] = {}
    for f in files:
        df = pd.read_parquet(f, columns=["mouse", "lobe", "X", "Y", "Gene"])
        lab = df[df.Gene.notna()][["X", "Y"]].to_numpy()
        axy = df[["X", "Y"]].to_numpy()
        near = cKDTree(lab).query(axy, k=1, distance_upper_bound=RING_EXCL)[0] < np.inf if len(lab) else np.zeros(len(df), bool)
        atree = cKDTree(axy)
        for code, poly in polygons(str(df.mouse.iloc[0]), str(df.lobe.iloc[0])):
            n["polygons_named"] += 1
            if code not in code2gene:
                unmapped[code] = unmapped.get(code, 0) + 1
                continue
            n["polygons_with_gene"] += 1
            inside = Poly(poly).contains_points(lab) if len(lab) else np.zeros(0, bool)
            if inside.sum() < 30:
                continue
            n["polygons_ge30_labelled"] += 1
            pts = lab[inside]
            cen, rad = pts.mean(0), np.max(np.hypot(*(pts - pts.mean(0)).T))
            ring = np.asarray(atree.query_ball_point(cen, r=rad + RING_OUT), dtype=int)
            ring = ring[(np.hypot(*(axy[ring] - cen).T) >= rad + RING_IN) & ~near[ring]]
            if ring.size >= 20:
                n["polygons_ring_ok"] += 1
    ann_total = 0
    for fj in ANN.glob("*_Lymphoid_*curated.fix.json"):
        d = json.loads(fj.read_text())
        ann_total += sum(1 for x in (d if isinstance(d, list) else d.get("features", [])) if isinstance(x, dict) and (x.get("properties") or {}).get("name"))
    n["all_named_annotations_in_files"] = ann_total
    n["unmapped_codes"] = unmapped
    (ROOT / "reports/step38_lesion_attrition.json").write_text(json.dumps(n, indent=2))
    print(json.dumps(n, indent=2))


if __name__ == "__main__":
    main()
