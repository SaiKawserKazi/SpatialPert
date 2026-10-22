"""AMENDMENT-008/008b adoption: literal and corrected size-adjusted power from per-draw p-values.

Writes reports/step19_adoption.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "reports/step19_stage1b.csv")
TESTS = ["loo_cell", "loo_patch30", "loo_lesion", "loo_stud_lesion", "loo_wcb", "loo_cr2", "loo_signflip", "loo_mixedlm",
         "circ_cell", "circ_patch30", "circ_lesion"]


def wilson(k: int, n: int) -> list[float]:
    z = 1.959964
    p = k / n
    den = 1 + z**2 / n
    mid = (p + z**2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return [float(mid - half), float(mid + half)]


out: dict[str, object] = {"per_test": {}}
for t in TESTS:
    rows = []
    for L in (1, 3, 5, 10, 20):
        n = df[(df.lesions == L) & (df["lambda"] == 0)][t].dropna()
        a = df[(df.lesions == L) & (df["lambda"] == 0.5)][t].dropna()
        if len(n) == 0:
            rows.append({"lesions": L, "applicable": 0.0})
            continue
        q05 = float(np.quantile(n, 0.05))
        crit = min(q05, 0.05)
        k = int((n < 0.05).sum())
        ka = int((a <= crit).sum()) if len(a) else 0
        rows.append({
            "lesions": L, "applicable": round(len(n) / 300, 3), "type1": k / len(n), "type1_ci": wilson(k, len(n)),
            "power_raw": float((a < 0.05).mean()) if len(a) else None,
            "power_size_adjusted_literal": float((a <= q05).mean()) if len(a) else None,
            "power_size_adjusted": ka / len(a) if len(a) else None, "power_ci": wilson(ka, len(a)) if len(a) else None,
            "m_null": int(len(n)), "m_alt": int(len(a)),
        })
    out["per_test"][t] = rows  # type: ignore[index]

per = out["per_test"]  # type: ignore[assignment]
worst = {t: max((r["type1"] for r in per[t] if "type1" in r), default=None) for t in TESTS}
qual = [t for t in TESTS if t.startswith("loo_") and worst[t] is not None and worst[t] <= 0.06]
mean_pow = {t: float(np.mean([r["power_size_adjusted"] for r in per[t] if r.get("power_size_adjusted") is not None])) for t in qual}
mean_lit = {t: float(np.mean([r["power_size_adjusted_literal"] for r in per[t] if r.get("power_size_adjusted_literal") is not None])) for t in qual}
out.update({"max_type1": worst, "qualifying": qual, "mean_power_corrected": mean_pow, "mean_power_literal": mean_lit,
            "adopted": max(mean_pow, key=lambda k: mean_pow[k]), "adopted_under_literal_rule": max(mean_lit, key=lambda k: mean_lit[k]),
            "frac_relabelled_by_leave_out": {int(L): float(df[(df.lesions == L) & (df["lambda"] == 0)]["loo_frac_relabelled"].mean()) for L in (1, 3, 5, 10, 20)}})
(ROOT / "reports/step19_adoption.json").write_text(json.dumps(out, indent=2))
print({k: out[k] for k in ("qualifying", "mean_power_corrected", "adopted", "adopted_under_literal_rule")})
