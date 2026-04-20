from __future__ import annotations
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = REPO_ROOT / "test-data"


def run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-m", "mags_visualization.main"] + args
    return subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


def test_sample_heatmap_requires_coverm() -> None:
    result = run_cli(
        [
            "sample-heatmap",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_sample_heatmap_requires_gtdb() -> None:
    result = run_cli(
        [
            "sample-heatmap",
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_gtdb_when_coloring_by_tax() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "tax",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_metadata_when_coloring_by_meta() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "meta",
            "--meta_col", "Infection by Nosema ceranae",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_meta_col_when_coloring_by_meta() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "meta",
            "--metadata", str(TEST_DATA / "metadata.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_drep_cluster_requires_drep() -> None:
    result = run_cli(
        [
            "drep-cluster-annot",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_taxa_sankey_requires_gtdb() -> None:
    result = run_cli(
        [
            "taxa-sankey",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_functional_annotation_requires_kegg_pathway_completeness() -> None:
    result = run_cli(
        [
            "drep-cluster-func",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_functional_annotation_requires_drep() -> None:
    result = run_cli(
        [
            "drep-cluster-func",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--kegg_pathway_completeness", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_all_requires_required_inputs() -> None:
    result = run_cli(
        [
            "all",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--kegg_pathway_completeness", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0