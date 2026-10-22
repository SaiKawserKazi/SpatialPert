"""Review Stage 2 (AMENDMENT-003): lesion number as a benchmark factor.

Pseudo-knockouts are L contiguous blobs ('lesions') of real control tumour cells, n cells in
total. The niche effect acts on each cell's own niche. Tests:
  patch30  - current default: permute niche labels among 30 µm single-linkage patches
  lesion   - permute niche labels among lesions (100 µm single linkage), majority niche
  cell     - naive cell-level permutation
An extra arm aligns the niche-specific effect with a real programme (hallmark hypoxia genes).
Writes reports/step12_lesions.csv / .json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from spatialpert.benchmark import cell_permutation_p, unit_rms
from spatialpert.effects import (
    clone_patches,
    fit_control_model,
    patch_domains,
    patch_permutation_test,
    residuals,
)
from spatialpert.pipeline import load_xenium_analysis

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
LS = [1, 3, 5, 10, 20]
NS = [400, 800]
LAMBDAS = [0.0, 0.5, 1.0]
REPS = {0.0: 300, 0.5: 150, 1.0: 150}
N_PERM = 250
DELTA = 0.5
MIN_SEP = 400.0  # µm between lesion seeds


def grow_lesions(tree: cKDTree, coords: np.ndarray, n: int, n_les: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    used = np.zeros(coords.shape[0], dtype=bool)
    seeds: list[int] = []
    tries = 0
    while len(seeds) < n_les and tries < 5000:
        tries += 1
        s = int(rng.integers(coords.shape[0]))
        if all(np.hypot(*(coords[s] - coords[t])) > MIN_SEP for t in seeds):
            seeds.append(s)
    per = max(1, n // n_les)
    idx, les = [], []
    for j, s in enumerate(seeds):
        _, nn = tree.query(coords[s], k=min(per * 3 + 10, coords.shape[0]))
        nn = np.atleast_1d(nn)
        nn = nn[~used[nn]][:per]
        used[nn] = True
        idx.append(nn)
        les.append(np.full(nn.size, j))
    return np.concatenate(idx), np.concatenate(les)


def main() -> None:
    a = load_xenium_analysis(ROOT)
    ctrl = a.is_control
    dom = a.cells["domain"].to_numpy()
    coords = a.cells[["x", "y"]].to_numpy()
    model = fit_control_model(a.y[ctrl], dom[ctrl], a.cov[ctrl])
    ci = np.where(ctrl)[0]
    sub = ci[:20000]
    gene_sd = residuals(model, a.y[sub], dom[sub], a.cov[sub]).std(axis=0)
    hall = (ROOT / "data/external/MSigDB_Hallmark_2020.txt").read_text().splitlines()
    hyp = next(set(x for x in ln.split("\t")[2:] if x) for ln in hall if ln.split("\t")[0] == "Hypoxia")
    hyp_mask = np.array([g in hyp for g in a.genes], dtype=float)
    tree = cKDTree(coords[ci])

    rows = []
    rng = np.random.default_rng(12)
    for n in NS:
        for n_les in LS:
            for lam in LAMBDAS:
                for rep in range(REPS[lam]):
                    for direction in (("random", "hypoxia") if lam > 0 and n == 400 else ("random",)):
                        il, les = grow_lesions(tree, coords[ci], n, n_les, rng)
                        idx = ci[il]
                        d = dom[idx]
                        r = residuals(model, a.y[idx], d, a.cov[idx])
                        if lam > 0:
                            u = unit_rms(rng.normal(size=r.shape[1]))
                            v = unit_rms(hyp_mask + 1e-9) if direction == "hypoxia" else unit_rms(rng.normal(size=r.shape[1]))
                            target = int(rng.choice(np.unique(d)))
                            r = r + DELTA * gene_sd * u
                            r[d == target] += lam * DELTA * gene_sd * v
                        xy = coords[idx]
                        p30 = clone_patches(xy, link=30.0)
                        pd30 = patch_domains(p30, d)
                        l100 = clone_patches(xy, link=100.0)
                        pdl = patch_domains(l100, d)
                        seed = int(rng.integers(1 << 31))
                        row = {"n": n, "lesions": n_les, "lambda": lam, "direction": direction, "rep": rep,
                               "n_patch30": int(pd30.size), "n_lesion_units": int(pdl.size)}
                        row["patch30"] = bool(patch_permutation_test(r, p30, pd30, n_perm=N_PERM, seed=seed).p_value < 0.05) if np.unique(pd30).size >= 2 else None
                        row["lesion"] = bool(patch_permutation_test(r, l100, pdl, n_perm=N_PERM, seed=seed).p_value < 0.05) if np.unique(pdl).size >= 2 else None
                        row["cell"] = bool(cell_permutation_p(r, d, n_perm=N_PERM, seed=seed) < 0.05)
                        rows.append(row)
                print("n", n, "L", n_les, "lambda", lam, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(R / "step12_lesions.csv", index=False)

    def rate(s: pd.Series) -> dict[str, float]:
        v = s.dropna().astype(bool)
        return {"rate": float(v.mean()) if len(v) else float("nan"), "m": int(len(v)), "untestable": int(s.isna().sum())}

    summ = [
        {"n": int(n), "lesions": int(nl), "lambda": float(lam), "direction": dirn, **{t: rate(s[t]) for t in ("patch30", "lesion", "cell")}}
        for (n, nl, lam, dirn), s in df.groupby(["n", "lesions", "lambda", "direction"])
    ]
    (R / "step12_lesions.json").write_text(json.dumps({"settings": {"ls": LS, "ns": NS, "reps": {str(k): v for k, v in REPS.items()}, "n_perm": N_PERM}, "summary": summ}, indent=2))
    print("done", len(df))


if __name__ == "__main__":
    main()
