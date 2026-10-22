"""Independent recomputation of headline numbers via separate code paths.

Each check re-derives a reported number from raw or intermediate files without using the
function that produced it, and compares. Writes reports/recheck.json; exits 1 on any mismatch.
"""

import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reports"
checks: list[dict[str, object]] = []


def check(name: str, reported: object, recomputed: object, tol: float = 1e-6) -> None:
    if isinstance(reported, float) or isinstance(recomputed, float):
        ok = abs(float(reported) - float(recomputed)) <= tol * max(1.0, abs(float(reported)))  # type: ignore[arg-type]
    else:
        ok = reported == recomputed
    checks.append({"check": name, "reported": reported, "recomputed": recomputed, "ok": bool(ok)})


def gate1() -> None:
    rep = json.loads((R / "gate1_data_audit.json").read_text())["knockout_totals"]
    cnt: Counter[str] = Counter()
    for f in (ROOT / "data/raw/GSE193460").glob("*spot_annotation.csv.gz"):
        with gzip.open(f, "rt") as fh:
            for row in csv.DictReader(fh):
                lab = row["phenotypes"]
                for g in ("Tgfbr2", "Ifngr2", "Jak2"):
                    if lab.startswith(g + "_"):
                        cnt[g] += 1
    for g in ("Tgfbr2", "Ifngr2", "Jak2"):
        check(f"gate1 {g} spots", rep[g], cnt[g])


def gate3_xenium() -> None:
    rep = json.loads((R / "gate3_dbit_audit.json").read_text())
    probes = {p["probe"]: p["cells"] for p in rep["Xenium HT-29 lung (GSE319123, single cell)"]["probes"]}
    with h5py.File(ROOT / "data/raw/GSE319123/GSM9511110_Xenium_lung_cell_feature_matrix.h5") as fh:
        ids = [x.decode() for x in fh["matrix/features/id"][:]]
        idx = fh["matrix/indices"][:]  # feature index of each nonzero (CSC by cell)
        data = fh["matrix/data"][:]
    for name in ("DIO3-4", "PRSS23-2", "ESAM-4", "CIITA-1"):
        col = ids.index(name)
        # one nonzero per (cell, feature): count cells with data>0 at this feature
        check(f"xenium {name} positive cells", probes[name], int(np.sum((idx == col) & (data > 0))))


def step4() -> None:
    rep = json.loads((R / "step4_domains.json").read_text())
    cells = pd.read_parquet(ROOT / "data/interim/xenium_cells.parquet")
    check("step4 tumour cells", rep["tumour_cells"], int(cells["tumour"].sum()))
    t = cells[cells["tumour"] & (cells["domain"] >= 0)]
    tab = pd.crosstab(t["label"], t["domain"])
    tab = tab.drop(index=[x for x in ("control", "multi") if x in tab.index])
    testable = sorted(tab.index[(tab >= 30).sum(axis=1) >= 2])
    check("step4 testable knockouts", rep["testable_knockouts"], testable)


def h1_naive() -> None:
    """Recompute observed T_k with explicit loops from the saved domain effects."""
    rep = json.loads((R / "step6_hypotheses.json").read_text())
    eff = pd.read_csv(R / "step6_effects.csv")
    genes = [c for c in eff.columns if c not in ("knockout", "domain", "cells")]
    total = 0.0
    for ko, sub in eff.groupby("knockout"):
        n = sub["cells"].to_numpy(float)
        e = sub[genes].to_numpy()
        grand = [sum(n[d] * e[d, g] for d in range(len(n))) / n.sum() for g in range(len(genes))]
        t = 0.0
        for d in range(len(n)):
            for g in range(len(genes)):
                t += n[d] * (e[d, g] - grand[g]) ** 2
        check(f"H1 T[{ko}] naive", rep["H1"]["per_knockout_exploratory"][ko]["T"], t, tol=1e-6)
        total += t
    check("H1 T pooled naive", rep["H1"]["T_pooled"], total, tol=1e-6)


def external_review() -> None:
    """Stage-1 benchmark (independent clones) rates recomputed from the per-replicate CSV."""
    d10 = json.loads((R / "step10_stage1.json").read_text())
    c10 = pd.read_csv(R / "step10_stage1.csv")
    for r in d10["summary"]:
        s = c10[(c10.n_target == r["n"]) & (c10["lambda"] == r["lambda"])]
        for m in ("patch", "studentized", "pseudobulk", "sandwich", "tile", "cell", "anova"):
            v = s[m].dropna()
            if len(v):
                check(f"stage1 n={r['n']} lam={r['lambda']} {m}", r[m]["rate"], float(v.astype(bool).mean()))



def _wilson(k: int, n: int) -> tuple[float, float]:
    from statsmodels.stats.proportion import proportion_confint  # noqa: PLC0415

    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return float(lo), float(hi)


def review2() -> None:
    """Second external review: stage-1b adoption, realism, brain, Perturb-map and real-data lesions."""
    ad = json.loads((R / "step19_adoption.json").read_text())
    df = pd.read_csv(R / "step19_stage1b.csv")
    for t, rows in ad["per_test"].items():
        for r in rows:
            if "type1" not in r:
                continue
            nul = df.loc[(df.lesions == r["lesions"]) & (df["lambda"] == 0.0), t].dropna().to_numpy()
            alt = df.loc[(df.lesions == r["lesions"]) & (df["lambda"] == 0.5), t].dropna().to_numpy()
            k = int(np.count_nonzero(nul < 0.05))
            check(f"step19 {t} L={r['lesions']} type1", r["type1"], k / nul.size)
            check(f"step19 {t} L={r['lesions']} type1 CI low (statsmodels)", r["type1_ci"][0], _wilson(k, nul.size)[0])
            if alt.size:
                crit = min(float(np.sort(nul)[int(np.floor(0.05 * (nul.size - 1)))] if nul.size == 1 else np.percentile(nul, 5)), 0.05)
                check(f"step19 {t} L={r['lesions']} size-adjusted power", r["power_size_adjusted"], float(np.count_nonzero(alt <= crit)) / alt.size)
    qual = sorted(t for t, v in ad["max_type1"].items() if t.startswith("loo_") and v is not None and v <= 0.06)
    check("step19 qualifying tests", sorted(ad["qualifying"]), qual)
    check("step19 adopted test", ad["adopted"], max(qual, key=lambda t: ad["mean_power_corrected"][t]))
    re_ = json.loads((R / "step26_realism.json").read_text())
    dr = pd.read_csv(R / "step26_realism.csv")
    for r in re_["summary"]:
        v = dr.loc[(dr.lesions == r["lesions"]) & (dr.delta == 0), "lesion"].dropna().to_numpy()
        check(f"step26 L={r['lesions']} lesion type1", r["lesion"]["type1"], float(np.count_nonzero(v < 0.05)) / v.size)
    br = json.loads((R / "step22_brain.json").read_text())
    db = pd.read_csv(R / "step22_brain.csv")
    for r in br["summary"]:
        for t in ("cell", "patch30", "lesion", "wcb"):
            v = db.loc[(db.lesions == r["lesions"]) & (db["lambda"] == 0), t].dropna().to_numpy()
            if r[t]["type1"] is not None:
                check(f"step22 L={r['lesions']} {t} type1", r[t]["type1"], float(np.count_nonzero(v < 0.05)) / v.size)
    pm = json.loads((R / "step25_perturbmap_lesions.json").read_text())
    names: list[str] = []
    for fn in sorted((ROOT / "data/raw/S-BIAD267").glob("*_Lymphoid_*curated.fix.json")):
        d = json.loads(fn.read_text())
        for x in d if isinstance(d, list) else d.get("features", []):
            nm = (x.get("properties") or {}).get("name") if isinstance(x, dict) else None
            if nm:
                head, _, tail = nm.rpartition("_")
                code = head if tail.isdigit() and head else nm
                base, _, suffix = code.rpartition("_")
                if base and re.fullmatch(r"M\d*", suffix):
                    continue  # ambiguous "_M" annotations, excluded (audit D3)
                names.append(code)
    cnt = Counter(names)
    check("perturbmap lesions", pm["lesions"], len(names))
    check("perturbmap codes", pm["distinct_procodes"], len(cnt))
    check("perturbmap median lesions/code", pm["lesions_per_code_median"], float(np.median(list(cnt.values()))))
    check("perturbmap codes >= 20", pm["codes_ge20"], sum(v >= 20 for v in cnt.values()))
    from scipy.cluster.hierarchy import fcluster, linkage  # noqa: PLC0415

    rl = json.loads((R / "step28_real_smallcluster.json").read_text())
    cells = pd.read_parquet(ROOT / "data/interim/xenium_cells.parquet")
    cells = cells[cells["tumour"]]
    from scipy.spatial import cKDTree  # noqa: PLC0415

    ref = cKDTree(cells.loc[(cells["label"] == "control") & (cells["domain"] >= 0), ["x", "y"]].to_numpy())
    for ko, v in rl.items():
        xy = cells.loc[cells["label"] == ko, ["x", "y"]].to_numpy()
        # niche needs at least one control neighbour within 100 um (else unlabelled, excluded)
        xy = xy[np.array([len(h) > 0 for h in ref.query_ball_point(xy, r=100.0)])]
        check(f"real cells[{ko}]", v["cells"], len(xy))
        n = int(fcluster(linkage(xy, "single"), t=100.0, criterion="distance").max()) if len(xy) > 1 else len(xy)
        check(f"real lesions[{ko}] via scipy linkage", v["lesions"], n)


def review3() -> None:
    """Review 3: calibration rates, representation benchmark, Perturb-map Psi and main effects."""
    cal = json.loads((R / "step31_calibrated_p.json").read_text())
    for key, csv_, tests in (("brain_step22", "step22_brain.csv", ("lesion", "patch30", "cell")),
                             ("tumour_step19", "step19_stage1b.csv", ("loo_lesion", "loo_patch30", "loo_cell"))):
        df = pd.read_csv(R / csv_)
        for t in tests:
            for L, rec in cal[key][t].items():
                v = df.loc[(df.lesions == int(L)) & (df["lambda"] == 0), t].dropna().to_numpy()
                check(f"step31 {key} {t} L={L} raw type1", rec["raw_type1"], float(np.count_nonzero(v < 0.05)) / v.size)
    for tag in ("", "_seed1", "_seed2"):
        js = json.loads((R / f"step35_representation{tag}.json").read_text())
        df = pd.read_csv(R / f"step35_representation{tag}.csv")
        for k, rec in js["summary"].items():
            L, rep = k[1:].split("_", 1)
            nul = df.loc[(df.lesions == int(L)) & (df["lambda"] == 0), rep].dropna().to_numpy()
            alt = df.loc[(df.lesions == int(L)) & (df["lambda"] == 0.5), rep].dropna().to_numpy()
            crit = min(float(np.percentile(nul, 5)), 0.05)
            check(f"step35{tag} {k} type1", rec["type1"], float(np.count_nonzero(nul < 0.05)) / nul.size)
            check(f"step35{tag} {k} size-adjusted power", rec["size_adjusted_power"], float(np.count_nonzero(alt <= crit)) / alt.size)
    sys.path.insert(0, str(ROOT / "scripts"))
    from step36_perturbmap_niche import build  # noqa: PLC0415 - lesion table only; contrasts recomputed below

    les = build()
    pm = json.loads((R / "step36_perturbmap.json").read_text())
    check("perturbmap lesions", pm["design"]["lesions"], len(les))
    check("perturbmap mice", pm["design"]["mice"], int(les.mouse.nunique()))
    main = json.loads((R / "step36b_main_effects.json").read_text())["y_CD8"]
    for y, tab in pm["results"].items():
        for g, rec in tab.items():
            other = les[les.gene != g]
            ref = other.groupby(["mouse", "niche"], as_index=False)[y].mean().rename(columns={y: "ref"})
            k = les[les.gene == g].merge(ref, on=["mouse", "niche"], how="left").dropna(subset=["ref"])
            e = k[y] - k["ref"]
            psi = float(e[k.niche == 1].mean() - e[k.niche == 0].mean())
            check(f"perturbmap {y} {g} psi (independent merge)", rec["psi"], psi, tol=1e-9)
            check(f"perturbmap {y} {g} lesions", rec["lesions"], int(len(k)))
            if y == "y_CD8":
                mm = k.assign(e=e).groupby("mouse").e.mean()
                check(f"perturbmap CD8 {g} main effect", main[g]["mean_of_mouse_means"], float(mm.mean()), tol=1e-9)
    ver = json.loads((R / "step37_verify_perturbmap.json").read_text())
    check("perturbmap raw extraction verified", ver["failed"], 0)


def main() -> None:
    for fn in (gate1, gate3_xenium, step4, h1_naive, external_review, review2, review3):
        try:
            fn()
        except FileNotFoundError as exc:
            checks.append({"check": fn.__name__, "ok": False, "error": f"missing {exc.filename}"})
    bad = [c for c in checks if not c["ok"]]
    out = {"n_checks": len(checks), "n_failed": len(bad), "failed": bad, "checks": checks}
    (R / "recheck.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"recheck: {len(checks) - len(bad)}/{len(checks)} OK")
    for c in bad:
        print("FAIL", c)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
