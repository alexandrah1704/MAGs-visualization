from __future__ import annotations

import sys
import shutil
import subprocess
from pathlib import Path

def run_and_check(cmd: list[str], cwd: Path, expected_files: list[Path]) -> bool:
    print("[TEST] Running:")
    print(" ".join(f'"{c}"' if " " in c else c for c in cmd))
    print()

    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        print(f"[ERROR] Command failed with exit code {result.returncode}")
        return False

    ok = True
    for path in expected_files:
        if not path.exists():
            print(f"[ERROR] Expected output file not found: {path}")
            ok = False
        elif path.stat().st_size == 0:
            print(f"[ERROR] Output file is empty: {path}")
            ok = False
        else:
            print(f"[OK] Output created: {path}")

    print()
    return ok


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    test_data = repo_root / "test-data"
    out_root = repo_root / "test-plots"

    if shutil.which("mags-visualization") is None:
        print("[ERROR] mags-visualization not found. Did you run 'pip install .'?")
        return 1

    out_root.mkdir(parents=True, exist_ok=True)

    python_cmd = [sys.executable, "-m", "mags_visualization.main"]

    tests = [
        {
            "name": "sample-heatmap",
            "cmd": python_cmd + [
                "sample-heatmap",
                "--coverm", str(test_data / "coverm.tsv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--metadata", str(test_data / "metadata.tsv"),
                "--meta_cols",
                "Infection by Nosema ceranae",
                "Chronic exposure to neonicotinoid",
                "Treatment with probiotic",
                "--tax_level", "phylum",
                "--output", str(out_root / "sample-heatmap"),
            ],
            "expected": [
                out_root / "sample-heatmap" / "heatmap_with_bars_species_['Infection by Nosema ceranae', 'Chronic exposure to neonicotinoid', 'Treatment with probiotic'].png",
            ],
        },

        {
            "name": "comp-conta",
            "cmd": python_cmd + [
                "comp-conta",
                "--checkm", str(test_data / "checkm.tsv"),
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--color_by", "tax",
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--tax_level", "phylum",
                "--output", str(out_root / "comp-conta"),
            ],
            "expected": [
                out_root / "comp-conta" / "comp_conta_marginals_checkm2.png",
            ],
        },

        {
            "name": "taxa-sankey",
            "cmd": python_cmd + [
                "taxa-sankey",
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--output", str(out_root / "taxa-sankey"),
            ],
            "expected": [
                out_root / "taxa-sankey" / "sankey_plot.html",
            ],
        },

        {
            "name": "drep-cluster-annot",
            "cmd": python_cmd + [
                "drep-cluster",
                "--drep", str(test_data / "drep.csv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--top_n", "30",
                "--output", str(out_root / "drep-cluster"),
            ],
            "expected": [
                out_root / "drep-cluster" / "drep_cluster_top_30_phylum.png",
            ],
        },

        {
            "name": "drep-cluster-func",
            "cmd": python_cmd + [
                "functional-annotation",
                "--drep", str(test_data / "drep.csv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--kegg_pathway_completeness", str(test_data / "kegg_pathway_completeness.tsv"),
                "--top_n", "30",
                "--output", str(out_root / "functional-annotation"),
            ],
            "expected": [out_root / "functional-annotation" / "drep_cluster_functional_core_top30.png"]
        },

        {
            "name": "all",
            "cmd": python_cmd + [
                "all",
                "--coverm", str(test_data / "coverm.tsv"),
                "--checkm", str(test_data / "checkm.tsv"),
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--drep", str(test_data / "drep.csv"),
                "--metadata", str(test_data / "metadata.tsv"),
                "--kegg_pathway_completeness", str(test_data / "kegg_pathway_completeness.tsv"),
                "--meta_cols",
                "Infection by Nosema ceranae",
                "Chronic exposure to neonicotinoid",
                "Treatment with probiotic",
                "--color_by", "tax",
                "--tax_level", "phylum",
                "--top_n", "30",
                "--output", str(out_root / "all"),
            ],
            "expected": [
                out_root / "all" / "heatmap_with_bars_species_['Infection by Nosema ceranae', 'Chronic exposure to neonicotinoid', 'Treatment with probiotic'].png",
                out_root / "all" / "comp_conta_marginals_checkm2.png",
                out_root / "all" / "sankey_plot.html",
                out_root / "all" / "drep_cluster_top_30_phylum.png",
                out_root / "all" / "drep_cluster_functional_core_top30.png",
            ],
        },
    ]

    all_ok = True
    for test in tests:
        print(f"# ---- Running test: {test['name']} ---- #")
        if not run_and_check(test["cmd"], repo_root, test["expected"]):
            all_ok = False

    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())