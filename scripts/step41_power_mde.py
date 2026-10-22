"""V5 P3: power, usable-result probability and minimum detectable effect (MDE).

Pseudo-knockouts (400 cells, L lesions) get a niche-specific effect lambda * v (v a random
unit-RMS direction, residual-SD units) in a target niche. Every draw records whether a lesion
test is returned (at least two niches among lesions), the lesion permutation p, the p
calibrated against a random-placement pool and against a niche-matched pool, the lesion
occupancy (for abstention diagnostics), and the CR2 minimum q.
Balance "obs": random seeds. Balance "rare": seeds depleted 10-fold in the target niche, so
the effect sits in a rarely visited niche.

Usage: step41_power_mde.py TISSUE L LAMBDA[,LAMBDA...] BALANCE N_DRAWS
Writes reports/v5/step41_<tissue>_L<L>_<balance>.csv (appended per draw).
"""

import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step40_scenarios import OUT, Pools, Tissue, grow, lesion_majority, sizes_equal, ss_perm_p  # noqa: E402

from spatialpert.calibration import calibrate_p
from spatialpert.smallcluster import cr2_two_niche_min_q


def effect_draw(t: Tissue, L: int, lam: float, balance: str, rng: np.random.Generator) -> dict[str, float]:
    target = int(t.dom[rng.integers(t.k)])
    w = np.where(t.niche == target, 0.1, 1.0) if balance == "rare" else None
    idx, les = grow(t, sizes_equal(L), rng, weights=w)
    d = t.leave_out(idx)
    keep = d >= 0
    idx, d, les = idx[keep], d[keep], np.unique(les[keep], return_inverse=True)[1]
    r = t.resid(idx, d)
    v = rng.normal(size=r.shape[1])
    v /= np.sqrt((v**2).mean())
    r[d == target] += lam * v
    um = lesion_majority(les, d)
    return {"target": target, "target_lesions": int((um == target).sum()), "n_units": int(les.max() + 1),
            "occupancy": ",".join(str(int((um == k).sum())) for k in t.dom),
            "p_lesion": ss_perm_p(r, les, um, rng), "q_cr2": cr2_two_niche_min_q(r, d, les)}


def main() -> None:
    tissue, L, lams, balance, n = sys.argv[1], int(sys.argv[2]), [float(x) for x in sys.argv[3].split(",")], sys.argv[4], int(sys.argv[5])
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    t = Tissue(tissue)
    pools = Pools(t, L, np.random.default_rng(9100 + L + (0 if tissue == "tumour" else 500)))
    f = OUT / f"step41_{tissue}_L{L}_{balance}.csv"
    have = pd.read_csv(f) if f.exists() else pd.DataFrame(columns=["lambda"])
    for lam in lams:
        done = int((have["lambda"] == lam).sum())
        rng = np.random.default_rng(zlib.crc32(f"{tissue}-{L}-{balance}-{lam}".encode()) + done)
        for i in range(done, n):
            rec = effect_draw(t, L, lam, balance, rng)
            p = rec["p_lesion"]
            rec["p_cal_random"] = calibrate_p(p, pools.random_pool("equal"), seed=i)
            rec["p_cal_matched"] = calibrate_p(p, pools.matched_pool("equal", rec["occupancy"], sizes_equal(L)), seed=i)
            rec |= {"tissue": tissue, "L": L, "balance": balance, "lambda": lam, "draw": i}
            pd.DataFrame([rec]).to_csv(f, mode="a", header=not f.exists(), index=False)
            if i % 50 == 0:
                print(f"lam={lam} {i}/{n} {time.time() - t0:.0f}s matched_pools={len(pools.matched)}", flush=True)
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
