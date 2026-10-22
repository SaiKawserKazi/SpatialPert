"""Priority 1 remedy: calibrate lesion-test p-values against matched control pseudo-knockouts.

p_cal = (#{null p_j < p_obs} + U * (1 + #{null p_j = p_obs})) / (1 + M), U ~ Uniform(0, 1), with M null p-values from control pseudo-knockouts
of the same lesion number in the same tissue (exact under the assumption that knockout lesions
are placed like random control positions). Evaluated by two-fold cross-fitting of the existing
null draws (calibration pool and test draws never overlap) and on the effect draws.
Writes reports/step31_calibrated_p.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"


def calibrate(p_obs: np.ndarray, pool: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Randomised rank of p_obs in the pool (ties broken uniformly), exact under exchangeability."""
    pool = np.sort(pool)
    below = np.searchsorted(pool, p_obs, side="left")
    ties = np.searchsorted(pool, p_obs, side="right") - below
    return (below + rng.uniform(size=p_obs.size) * (ties + 1)) / (1 + pool.size)


def wilson(k: int, n: int) -> list[float]:
    from statsmodels.stats.proportion import proportion_confint  # noqa: PLC0415

    lo, hi = proportion_confint(k, n, method="wilson")
    return [float(lo), float(hi)]


def evaluate(nul: np.ndarray, alt: np.ndarray, rng: np.random.Generator) -> dict[str, object]:
    nul = nul[~np.isnan(nul)]
    alt = alt[~np.isnan(alt)] if alt is not None else np.array([])
    perm = rng.permutation(nul.size)
    half = nul.size // 2
    a, b = nul[perm[:half]], nul[perm[half:]]
    pc = np.r_[calibrate(a, b, rng), calibrate(b, a, rng)]  # cross-fitted
    k = int((pc <= 0.05).sum())
    out: dict[str, object] = {"m_null": int(nul.size), "raw_type1": float((nul < 0.05).mean()),
                              "calibrated_type1": k / pc.size, "calibrated_type1_ci": wilson(k, pc.size)}
    if alt.size:
        pa = calibrate(alt, nul, rng)
        out |= {"m_alt": int(alt.size), "raw_power": float((alt < 0.05).mean()), "calibrated_power": float((pa <= 0.05).mean())}
    return out


def main() -> None:
    rng = np.random.default_rng(31)
    res: dict[str, object] = {}
    br = pd.read_csv(R / "step22_brain.csv")
    res["brain_step22"] = {t: {int(L): evaluate(br[(br.lesions == L) & (br["lambda"] == 0)][t].to_numpy(),
                                                br[(br.lesions == L) & (br["lambda"] == 0.5)][t].to_numpy(), rng)
                               for L in (3, 5, 10, 20)} for t in ("lesion", "patch30", "cell")}
    tu = pd.read_csv(R / "step19_stage1b.csv")
    res["tumour_step19"] = {t: {int(L): evaluate(tu[(tu.lesions == L) & (tu["lambda"] == 0)][t].to_numpy(),
                                                 tu[(tu.lesions == L) & (tu["lambda"] == 0.5)][t].to_numpy(), rng)
                                for L in (3, 5, 10, 20)} for t in ("loo_lesion", "loo_patch30", "loo_cell")}
    (R / "step31_calibrated_p.json").write_text(json.dumps(res, indent=2))
    for k, v in res.items():
        for t, rows in v.items():  # type: ignore[union-attr]
            print(k, t, {L: (round(r["raw_type1"], 3), round(r["calibrated_type1"], 3), round(r.get("calibrated_power", float("nan")), 2)) for L, r in rows.items()})


if __name__ == "__main__":
    main()
