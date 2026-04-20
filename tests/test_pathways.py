from __future__ import annotations
import pandas as pd
import pytest

from mags_visualization.drep_cluster_func import prepare_pathway_matrix, select_top_modules


def make_pathway_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "contig_normalized": ["g1_fasta", "g1_fasta", "g2_fasta", "g2_fasta", "g3_fasta", "g3_fasta"],
            "module_accession": ["M1", "M2", "M1", "M2", "M1", "M2"],
            "completeness": [100, 10, 90, 80, 95, 20],
        }
    )


def make_full_pathway_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "contig": ["g1", "g1", "g2", "g2"],
            "module_accession": ["M1", "M2", "M1", "M2"],
            "completeness": [100, 20, 80, 50],
            "pathway_name": ["Path1", "Path2", "Path1", "Path2"],
            "pathway_class": [
                "Metabolism; Carbohydrate metabolism",
                "Metabolism; Energy metabolism",
                "Metabolism; Carbohydrate metabolism",
                "Metabolism; Energy metabolism",
            ],
        }
    )


def test_select_top_modules_none_returns_all() -> None:
    pw = make_pathway_df()
    reps = ["g1_fasta", "g2_fasta", "g3_fasta"]
    result = select_top_modules(pw, reps, top_modules=None)
    assert len(result) == len(pw)


def test_select_top_modules_mean() -> None:
    pw = make_pathway_df()
    reps = ["g1_fasta", "g2_fasta", "g3_fasta"]
    result = select_top_modules(pw, reps, top_modules=1, mode="mean")
    assert set(result["module_accession"]) == {"M1"}


def test_select_top_modules_variance() -> None:
    pw = pd.DataFrame(
        {
            "contig_normalized": ["g1_fasta", "g1_fasta", "g2_fasta", "g2_fasta", "g3_fasta", "g3_fasta"],
            "module_accession": ["M1", "M2", "M1", "M2", "M1", "M2"],
            "completeness": [50, 0, 50, 100, 50, 50],
        }
    )
    reps = ["g1_fasta", "g2_fasta", "g3_fasta"]
    result = select_top_modules(pw, reps, top_modules=1, mode="variance")
    assert set(result["module_accession"]) == {"M2"}


def test_select_top_modules_no_mode() -> None:
    pw = make_pathway_df()
    reps = ["g1_fasta", "g2_fasta", "g3_fasta"]
    with pytest.raises(ValueError, match="mode must be"):
        select_top_modules(pw, reps, top_modules=2, mode="invalid")


def test_prepare_pathway_matrix_missing_columns() -> None:
    df = pd.DataFrame({"contig": ["g1"]})
    with pytest.raises(ValueError, match="missing required columns"):
        prepare_pathway_matrix(df, representatives=["g1_fasta"])


def test_prepare_pathway_matrix_no_matching_rows() -> None:
    df = make_full_pathway_df()
    with pytest.raises(ValueError, match="No pathway rows matched"):
        prepare_pathway_matrix(df, representatives=["not_present_fasta"])


def test_prepare_pathway_matrix_right() -> None:
    df = make_full_pathway_df()
    heatmap_df, module_meta, class_blocks = prepare_pathway_matrix(
        df,
        representatives=["g1_fasta", "g2_fasta"],
        top_modules=None,
    )

    assert heatmap_df.shape[0] == 2
    assert list(heatmap_df.index) == ["g1_fasta", "g2_fasta"]
    assert len(module_meta) == 2
    assert len(class_blocks) >= 1
    assert "module_accession" in module_meta.columns
    assert "pathway_class_simple" in module_meta.columns