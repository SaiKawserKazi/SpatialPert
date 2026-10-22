"""Download the Perturb-map Visium files needed for the week-1 gate (GSE193460)."""

from pathlib import Path
from urllib.request import urlretrieve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "GSE193460"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5808nnn"
SAMPLES = {"GSM5808054": "KP_1", "GSM5808055": "KP_2", "GSM5808056": "KP_3", "GSM5808057": "KP_4"}
FILES = [
    "filtered_feature_bc_matrix.h5",
    "spot_annotation.csv.gz",
    "tissue_positions_list.csv.gz",
    "scalefactors_json.json.gz",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for gsm, section in SAMPLES.items():
        for suffix in FILES:
            name = f"{gsm}_{section}_{suffix}"
            dest = OUT / name
            if dest.exists():
                continue
            url = f"{BASE}/{gsm}/suppl/{name}"
            print(f"fetch {name}")
            urlretrieve(url, dest)  # noqa: S310 - fixed https GEO URL


if __name__ == "__main__":
    main()
