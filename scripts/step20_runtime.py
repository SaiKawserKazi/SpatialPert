"""A9: runtime and peak memory on the full Xenium section (352,913 cells), hardware stated.

Run: uv run --with psutil python scripts/step20_runtime.py
"""

import json
import os
import platform
import time
import tracemalloc
from pathlib import Path

import numpy as np
import psutil

from spatialpert.api import location_dependence
from spatialpert.effects import fit_control_model, residuals
from spatialpert.pipeline import NicheModel, load_xenium_analysis
from spatialpert.smallcluster import wild_cluster_bootstrap_p, within_lesion_signflip_p

ROOT = Path(__file__).resolve().parents[1]
proc = psutil.Process(os.getpid())


def timed(fn, *args, **kw):  # noqa: ANN001, ANN002, ANN003, ANN201
    tracemalloc.start()
    t = time.perf_counter()
    out = fn(*args, **kw)
    dt = time.perf_counter() - t
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return out, dt, peak / 1e6


def main() -> None:
    res: dict[str, object] = {
        "hardware": {"platform": platform.platform(), "processor": platform.processor(), "logical_cpus": psutil.cpu_count(),
                     "ram_gb": round(psutil.virtual_memory().total / 1e9, 1), "python": platform.python_version(), "single_threaded_numpy": False},
    }
    a, dt, mb = timed(load_xenium_analysis, ROOT)
    res["load_and_normalise_full_section"] = {"cells_in_file": 352913, "analysis_cells": int(a.y.shape[0]), "seconds": dt, "peak_python_alloc_mb": mb}
    nm, dt, mb = timed(NicheModel.fit, a)
    res["niche_model_fit"] = {"seconds": dt, "peak_python_alloc_mb": mb}
    dom = a.cells["domain"].to_numpy()
    xy = a.cells[["x", "y"]].to_numpy()
    for unit in ("patch", "lesion"):
        out, dt, mb = timed(location_dependence, a.y, xy, a.cells["label"].to_numpy(), dom, covariates=a.cov,
                            knockouts=["CIITA", "DIO3", "ESAM", "PRSS23"], unit=unit, n_perm=1000)
        res[f"location_dependence_{unit}_1000perm_4KO"] = {"seconds": dt, "seconds_per_knockout": dt / 4, "peak_python_alloc_mb": mb}
    ctrl = a.is_control
    model = fit_control_model(a.y[ctrl], dom[ctrl], a.cov[ctrl])
    import sys  # noqa: PLC0415

    from scipy.spatial import cKDTree  # noqa: PLC0415

    sys.path.insert(0, str(ROOT / "scripts"))
    from step12_lesions import grow_lesions  # noqa: PLC0415

    ci = np.where(ctrl)[0]
    il, les = grow_lesions(cKDTree(xy[ci]), xy[ci], 400, 10, np.random.default_rng(0))
    idx = ci[il]
    r = residuals(model, a.y[idx], dom[idx], a.cov[idx])
    p, dt, mb = timed(wild_cluster_bootstrap_p, r, dom[idx], les, n_boot=999)
    res["wild_cluster_bootstrap_999_pseudoKO_10lesions_400cells"] = {"p": p, "seconds": dt, "peak_python_alloc_mb": mb}
    p, dt, mb = timed(within_lesion_signflip_p, r, dom[idx], les, n_flip=2000)
    res["signflip_2000_pseudoKO_10lesions_400cells"] = {"p": p, "seconds": dt, "peak_python_alloc_mb": mb}
    res["process_peak_rss_gb"] = round(proc.memory_info().peak_wset / 1e9, 2) if hasattr(proc.memory_info(), "peak_wset") else None
    (ROOT / "reports/step20_runtime.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
