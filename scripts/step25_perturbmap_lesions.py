"""A5/A4 context: lesions per knockout (Pro-Code) in the Perturb-map imaging deposit S-BIAD267.

Reads the curated lesion annotations (QuPath GeoJSON, lymphoid panel; one polygon per lesion,
named by its Pro-Code epitope combination). Writes reports/step25_perturbmap_lesions.json.

V5 fix (audit D3): 12 polygons carry an `_M`/`_M<digits>` suffix after the code (meaning not
documented). v4 counted them as 11 extra codes. They are now excluded as ambiguous and
reported separately ("ambiguous_polygons").
"""

import collections
import glob
import json
import re
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/S-BIAD267"


def main() -> None:
    rows = []
    ambiguous: list[str] = []
    for fn in sorted(glob.glob(str(RAW / "*_Lymphoid_*curated.fix.json"))):
        name = Path(fn).name
        mouse = name.split("_")[0]
        m = re.search(r"lobe(\d+)", name)
        lobe = m.group(1) if m else "?"
        d = json.loads(Path(fn).read_text())
        feats = d if isinstance(d, list) else d.get("features", [])
        for x in feats:
            if isinstance(x, dict) and x.get("properties", {}).get("name"):
                code = re.sub(r"_\d+$", "", x["properties"]["name"])
                if re.search(r"_M\d*$", code):
                    ambiguous.append(code)
                    continue
                rows.append((mouse, lobe, code))
    codes = collections.Counter(c for _, _, c in rows)
    per_mouse: dict[str, set[str]] = collections.defaultdict(set)
    for m_, _, c in rows:
        per_mouse[c].add(m_)
    v = sorted(codes.values(), reverse=True)
    out = {
        "source": "BioImage Archive S-BIAD267 (Dhainaut et al. 2022), *_Lymphoid_*curated.fix.json",
        "lesions": len(rows), "ambiguous_polygons": len(ambiguous), "polygons_total": len(rows) + len(ambiguous),
        "mice": sorted({m_ for m_, _, _ in rows}), "lobes": len({(m_, l_) for m_, l_, _ in rows}),
        "distinct_procodes": len(codes), "lesions_per_code_median": statistics.median(v), "lesions_per_code_max": v[0],
        "codes_ge10": sum(x >= 10 for x in v), "codes_ge20": sum(x >= 20 for x in v),
        "codes_in_ge2_mice": sum(len(s) >= 2 for s in per_mouse.values()),
        "per_code": {c: {"lesions": n, "mice": len(per_mouse[c])} for c, n in codes.most_common()},
    }
    (ROOT / "reports/step25_perturbmap_lesions.json").write_text(json.dumps(out, indent=2))
    print({k: v for k, v in out.items() if k != "per_code"})


if __name__ == "__main__":
    main()
