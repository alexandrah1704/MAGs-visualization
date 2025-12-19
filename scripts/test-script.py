from __future__ import annotations

import sys
import subprocess
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]

    # Ordner: test-data enthält bee-use-case Daten
    use_case_folder = repo_root / "test-data"  # z.B. test-data/coverm.tsv etc.
    out_dir = repo_root / "test-plots"

    main_py = repo_root / "scripts" / "main.py"

    # Basis-Command (entspricht deinem PowerShell Beispiel)
    cmd = [
        sys.executable,
        str(main_py),

        "--coverm", str(use_case_folder / "coverm.tsv"),
        "--checkm", str(use_case_folder / "checkm.tsv"),
        "--checkm2", str(use_case_folder / "checkm2.tsv"),
        "--gtdb", str(use_case_folder / "gtdb.tsv"),
        "--drep", str(use_case_folder / "drep.csv"),
        "--quast", str(use_case_folder / "quast.tsv"),
        "--bakta", str(use_case_folder / "bakta.tsv"),
        "--metadata", str(use_case_folder / "metadata.tsv"),

        # metadata columns (mehrere Werte!)
        "--meta_cols",
        "Infection by Nosema ceranae",
        "Chronic exposure to neonicotinoid",
        "Treatment with probiotic",

        "--color_by", "tax",
        "--tax_level", "phylum",

        # in deinem aktuellen main.py: -n/--top_n_counts (nicht --top_n)
        "--top_n", "30",

        "--top_bar_spacer", "-0.5",
        "--spacer_meta", "2.5",

        "-o", str(out_dir),
    ]

    print("[TEST] Running:")
    print(" ".join(f'"{c}"' if " " in c else c for c in cmd))
    print()

    out_dir.mkdir(parents=True, exist_ok=True)

    # cwd auf repo_root setzen, damit relative Imports/Pfade stabil sind
    result = subprocess.run(cmd, cwd=str(repo_root))
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
