"""Review 3, Priorities 2 and 4: does a learned context representation improve inference?

One downstream estimator (lesion permutation on residuals against the representation's niche
means) and four niche representations, all K = 4 k-means on control cells:
  kernel   fixed kernel-weighted neighbourhood PCs (current SpatialPert, NicheModel)
  learned  ContextEncoder (one-hop, two rings, self-masked, trained on controls)
  expr     the cell's own expression PCs (no context)
  coord    x-y coordinates only
Pseudo-knockout cells are left out of their own context (kernel, learned). The niche-dependent
effect follows a continuous covariate, the host-cell fraction within 50 µm (standardised), so
no partition is favoured by construction: effect = delta*sd*u + lambda*delta*sd*v*m_i.
Writes reports/step35_representation.csv / .json.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_lesions import grow_lesions  # noqa: E402

from spatialpert.benchmark import unit_rms
from spatialpert.context_gnn import ContextEncoder
from spatialpert.effects import fit_control_model, patch_domains, patch_permutation_test, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
DELTA = 0.5
N_PERM = 250
REPS = {0.0: 200, 0.5: 150}
LS = (10, 20)
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0  # encoder training seed (replication runs)
TAG = f"_seed{SEED}" if SEED else ""


def perm_p(r: np.ndarray, les: np.ndarray, niche: np.ndarray, seed: int) -> float:
    keep = niche >= 0
    les = np.unique(les[keep], return_inverse=True)[1]
    ud = patch_domains(les, niche[keep])
    if np.unique(ud).size < 2:
        return float("nan")
    return patch_permutation_test(r[keep], les, ud, n_perm=N_PERM, seed=seed).p_value


def main() -> None:
    t0 = time.time()
    a = load_xenium_analysis(ROOT)
    nm = NicheModel.fit(a)
    ok = nm.labels >= 0
    ci = np.where(a.is_control & ok)[0]
    xy, pcs = nm.coords, nm.pcs[:, :30]
    enc = ContextEncoder.fit(xy[ci], pcs[ci], seed=SEED)
    km = {
        "learned": KMeans(4, n_init=10, random_state=0).fit(enc.embed_reference()),
        "expr": KMeans(4, n_init=10, random_state=0).fit(pcs[ci]),
        "coord": KMeans(4, n_init=10, random_state=0).fit(xy[ci]),
    }
    ref_lab = {"kernel": nm.labels[ci], "learned": km["learned"].labels_, "expr": km["expr"].labels_, "coord": km["coord"].labels_}
    models = {k: fit_control_model(a.y[ci], v, a.cov[ci]) for k, v in ref_lab.items()}
    sd = residuals(models["kernel"], a.y[ci[:20000]], nm.labels[ci[:20000]], a.cov[ci[:20000]]).std(0)
    host = a.cov[:, 1]
    m_std = (host - host[ci].mean()) / host[ci].std()
    pos = {c: i for i, c in enumerate(ci)}

    def labels(rep: str, idx: np.ndarray) -> np.ndarray:
        if rep == "kernel":
            return nm.leave_out_labels(idx, idx)
        if rep == "learned":
            excl = np.zeros(ci.size, dtype=bool)
            excl[[pos[i] for i in idx]] = True
            return km["learned"].predict(enc.embed(xy[idx], exclude=excl)).astype(int)
        if rep == "expr":
            return km["expr"].predict(pcs[idx]).astype(int)
        return km["coord"].predict(xy[idx]).astype(int)

    stats = {k: {"largest_share": float(np.bincount(v).max() / v.size)} for k, v in ref_lab.items()}
    print("fitted", round(time.time() - t0), stats, flush=True)
    tree = cKDTree(xy[ci])
    rng = np.random.default_rng(35 + SEED)
    rows = []
    for n_les in LS:
        for lam, reps in REPS.items():
            for rep in range(reps):
                il, les = grow_lesions(tree, xy[ci], 400, n_les, rng)
                idx = ci[il]
                u = unit_rms(rng.normal(size=a.y.shape[1]))
                v = unit_rms(rng.normal(size=a.y.shape[1]))
                seed = int(rng.integers(1 << 31))
                row: dict[str, float] = {"lesions": n_les, "lambda": lam, "rep": rep}
                for name, model in models.items():
                    d = labels(name, idx)
                    dd = np.where(d >= 0, d, 0)
                    r = residuals(model, a.y[idx], dd, a.cov[idx])
                    r = r + DELTA * sd * u + lam * DELTA * sd * v * m_std[idx, None]
                    row[name] = perm_p(r, les, d, seed)
                rows.append(row)
            print(f"L={n_les} lam={lam} {time.time() - t0:.0f}s", flush=True)
            pd.DataFrame(rows).to_csv(R / f"step35_representation{TAG}.csv", index=False)
    df = pd.DataFrame(rows)
    out: dict[str, object] = {"partition": stats, "summary": {}}
    for L in LS:
        for name in models:
            nul = df[(df.lesions == L) & (df["lambda"] == 0)][name].dropna()
            alt = df[(df.lesions == L) & (df["lambda"] == 0.5)][name].dropna()
            crit = min(float(np.quantile(nul, 0.05)), 0.05)
            out["summary"][f"L{L}_{name}"] = {  # type: ignore[index]
                "type1": float((nul < 0.05).mean()), "raw_power": float((alt < 0.05).mean()),
                "size_adjusted_power": float((alt <= crit).mean()), "m_null": int(len(nul)), "m_alt": int(len(alt)),
            }
    (R / f"step35_representation{TAG}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
