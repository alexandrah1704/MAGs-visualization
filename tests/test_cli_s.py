from __future__ import annotations
import subprocess
import sys
from pathlib import Path
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = REPO_ROOT / "test-data"


def run_cli(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-m", "mags_visualization.main"] + args
    return subprocess.run(
        cmd,
        cwd=str(cwd or REPO_ROOT),
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    "args",
    [
        ["--help"],
        ["sample-heatmap", "--help"],
        ["comp-conta", "--help"],
        ["taxa-sankey", "--help"],
        ["drep-cluster-annot", "--help"],
        ["drep-cluster-func", "--help"],
        ["all", "--help"],
    ],
)
def test_help_commands_work(args: list[str]) -> None:
    result = run_cli(args)
    assert result.returncode == 0, result.stderr


def test_sample_heatmap_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "sample-heatmap"
    result = run_cli(
        [
            "sample-heatmap",
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--metadata", str(TEST_DATA / "metadata.tsv"),
            "--meta_cols",
            "Infection by Nosema ceranae",
            "Chronic exposure to neonicotinoid",
            "Treatment with probiotic",
            "--tax_level", "phylum",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "heatmap_with_bars_phylum_['Infection by Nosema ceranae', 'Chronic exposure to neonicotinoid', 'Treatment with probiotic'].png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_comp_conta_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "comp-conta"
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--tax_level", "phylum",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "comp_conta_marginals_checkm2.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_taxa_sankey_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "taxa-sankey"
    result = run_cli(
        [
            "taxa-sankey",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "sankey_plot.html"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_drep_cluster_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "drep-cluster-annot"
    result = run_cli(
        [
            "drep-cluster-annot",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--top_n", "30",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "drep_cluster_top30_phylum-genus.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_functional_annotation_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "drep-cluster-func"
    result = run_cli(
        [
            "drep-cluster-func",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--pathways", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--top_n", "30",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "drep_cluster_functional_core_top30.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_all_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "all"
    result = run_cli(
        [
            "all",
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--drep", str(TEST_DATA / "drep.csv"),
            "--metadata", str(TEST_DATA / "metadata.tsv"),
            "--pathways", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--meta_cols",
            "Infection by Nosema ceranae",
            "Chronic exposure to neonicotinoid",
            "Treatment with probiotic",
            "--tax_level", "phylum",
            "--top_n", "30",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr

    expected = [
        outdir / "heatmap_with_bars_phylum_['Infection by Nosema ceranae', 'Chronic exposure to neonicotinoid', 'Treatment with probiotic'].png",
        outdir / "comp_conta_marginals_checkm2.png",
        outdir / "sankey_plot.html",
        outdir / "drep_cluster_top30_phylum-genus.png",
        outdir / "drep_cluster_functional_core_top30.png",
    ]
    for path in expected:
        assert path.exists(), f"Missing expected output: {path}"
        assert path.stat().st_size > 0, f"Empty output file: {path}"