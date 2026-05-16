from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

def run_and_check(cmd: list[str], cwd: Path, expected_files: list[Path]) -> bool:
    print("[TEST] Running:", flush=True)
    print(" ".join(f'"{c}"' if " " in c else c for c in cmd), flush=True)
    print(flush=True)

    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        print(f"[ERROR] Command failed with exit code {result.returncode}")
        if result.stdout:
            print("[STDOUT]")
            print(result.stdout)
        if result.stderr:
            print("[STDERR]")
            print(result.stderr)
        print()
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
    
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    cli_md = ["mags-visualization"]

    tests = [
        {
            "name": "sample-heatmap",
            "cmd": cli_md + [
                "sample-heatmap",
                "--coverm", str(test_data / "coverm.tsv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--metadata", str(test_data / "metadata.tsv"),
                "--meta_cols",
                "Infection by Nosema ceranae",
                "Chronic exposure to neonicotinoid",
                "Treatment with probiotic",
                "--tax_level", "phylum",
                "--output", str(out_root),
            ],
            "expected": [
                out_root / "heatmap_with_bars_phylum.png",
            ],
        },

        {
            "name": "comp-conta",
            "cmd": cli_md + [
                "comp-conta",
                "--tools", "checkm2",
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--mode", "tax",
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--tax_level", "phylum",
                "--output", str(out_root),
            ],
            "expected": [
                out_root / "comp_conta_by_rank_marginals_CheckM2.png",
            ],
        },

        {
            "name": "taxa-sankey",
            "cmd": cli_md + [
                "taxa-sankey",
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--format", "html",
                "--output", str(out_root),
            ],
            "expected": [
                out_root / "sankey_plot.html",
            ],
        },

        {
            "name": "drep-cluster-annot",
            "cmd": cli_md + [
                "drep-cluster-annot",
                "--drep", str(test_data / "drep.csv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--quast", str(test_data / "quast.tsv"),
                "--bakta", str(test_data / "bakta.tsv"),
                "--top_n", "30",
                "--output", str(out_root),
            ],
            "expected": [
                out_root / "drep_cluster_top30_phylum-genus.png",
            ],
        },

        {
            "name": "drep-cluster-func",
            "cmd": cli_md + [
                "drep-cluster-func",
                "--drep", str(test_data / "drep.csv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--pathways", str(test_data / "kegg_pathway_completeness.tsv"),
                "--top_n", "30",
                "--output", str(out_root),
            ],
            "expected": [out_root / "drep_cluster_functional_core_top30.png"]
        },

        {
            "name": "all",
            "cmd": cli_md + [
                "all",
                "--coverm", str(test_data / "coverm.tsv"),
                "--checkm", str(test_data / "checkm.tsv"),
                "--checkm2", str(test_data / "checkm2.tsv"),
                "--gtdb", str(test_data / "gtdb.tsv"),
                "--drep", str(test_data / "drep.csv"),
                "--quast", str(test_data / "quast.tsv"),
                "--metadata", str(test_data / "metadata.tsv"),
                "--pathways", str(test_data / "kegg_pathway_completeness.tsv"),
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
                out_root / "all" / "heatmap_with_bars_phylum.png",
                out_root / "all" / "comp_conta_by_rank_marginals_CheckM2.png",
                out_root / "all" / "sankey_plot.html",
                out_root / "all" / "drep_cluster_top30_phylum-genus.png",
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