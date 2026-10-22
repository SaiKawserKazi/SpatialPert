"""Review minor point: H3 against random 8-gene sets (10,000 draws), raw and variance-scaled.

Writes reports/step15_h3_random.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
h6 = json.loads((R / "step6_hypotheses.json").read_text())
eff = pd.read_csv(R / "step6_effects.csv")
genes = [c for c in eff.columns if c not in ("knockout", "domain", "cells")]
energy = np.zeros(len(genes))
for _, s in eff.groupby("knockout"):
    e = s[genes].to_numpy()
    n = s["cells"].to_numpy(float)
    eb = (e * n[:, None]).sum(0) / n.sum()
    energy += (n[:, None] * (e - eb) ** 2).sum(0)
hall = np.isin(genes, h6["H3"]["hallmark_genes"])
k = int(hall.sum())
rng = np.random.default_rng(15)
obs = float(np.mean(energy[hall]))
null = np.array([np.mean(energy[rng.choice(len(genes), k, replace=False)]) for _ in range(10000)])
other = float(np.mean(energy[~hall]))
out = {"k_hallmark": k, "mean_energy_hallmark": obs, "mean_energy_other": other, "ratio": obs / other,
       "p_random_sets": float((1 + np.sum(null >= obs)) / 10001)}
(R / "step15_h3_random.json").write_text(json.dumps(out, indent=2))
print(out)
