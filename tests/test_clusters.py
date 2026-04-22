from __future__ import annotations
import pandas as pd
import pytest

from mags_visualization.drep_cluster_func import prepare_cluster_data


def test_prepare_cluster_data_missing_drep_columns() -> None:
    drep = pd.DataFrame({"wrong": [1]})
    gtdb = pd.DataFrame(
        {
            "user_genome": ["g1"],
            "classification": ["d__Bacteria; p__Firmicutes"],
        }
    )

    with pytest.raises(ValueError, match="dRep must contain"):
        prepare_cluster_data(drep, gtdb, tax_levels=["phylum"], top_n=5)


def test_prepare_cluster_data_missing_gtdb_columns() -> None:
    drep = pd.DataFrame(
        {
            "genome": ["g1"],
            "secondary_cluster": ["c1"],
        }
    )
    gtdb = pd.DataFrame({"wrong": [1]})

    with pytest.raises(ValueError, match="GTDB file must contain"):
        prepare_cluster_data(drep, gtdb, tax_levels=["phylum"], top_n=5)


def test_prepare_cluster_data_right() -> None:
    drep = pd.DataFrame(
        {
            "genome": ["g1", "g2", "g3"],
            "secondary_cluster": ["c1", "c1", "c2"],
        }
    )
    gtdb = pd.DataFrame(
        {
            "user_genome": ["g1", "g2", "g3"],
            "classification": [
                "d__Bacteria; p__Firmicutes",
                "d__Bacteria; p__Pseudomonadota",
                "d__Bacteria; p__Actinobacteriota",
            ],
        }
    )

    cluster_data, drep_out, total_clusters = prepare_cluster_data(
        drep,
        gtdb,
        tax_levels=["phylum"],
        top_n=2,
    )

    assert total_clusters == 2
    assert len(cluster_data) == 2
    assert "representative" in cluster_data.columns
    assert "representative_display" in cluster_data.columns
    assert "phylum" in cluster_data.columns
    assert isinstance(drep_out, pd.DataFrame)


def test_prepare_cluster_data_respects_top_n() -> None:
    drep = pd.DataFrame(
        {
            "genome": ["g1", "g2", "g3", "g4", "g5"],
            "secondary_cluster": ["c1", "c1", "c2", "c3", "c3"],
        }
    )
    gtdb = pd.DataFrame(
        {
            "user_genome": ["g1", "g2", "g3", "g4", "g5"],
            "classification": [
                "d__Bacteria; p__Firmicutes",
                "d__Bacteria; p__Firmicutes",
                "d__Bacteria; p__Pseudomonadota",
                "d__Bacteria; p__Actinobacteriota",
                "d__Bacteria; p__Actinobacteriota",
            ],
        }
    )

    cluster_data, _, _ = prepare_cluster_data(
        drep,
        gtdb,
        tax_levels=["phylum"],
        top_n=2,
    )

    assert len(cluster_data) == 2


def test_prepare_cluster_data_with_pathway_repr_in_pathway() -> None:
    drep = pd.DataFrame(
        {
            "genome": ["g1", "g2"],
            "secondary_cluster": ["c1", "c1"],
        }
    )
    gtdb = pd.DataFrame(
        {
            "user_genome": ["g1", "g2"],
            "classification": [
                "d__Bacteria; p__Firmicutes",
                "d__Bacteria; p__Pseudomonadota",
            ],
        }
    )
    pathway_df = pd.DataFrame(
        {
            "contig": ["g2"],
        }
    )

    cluster_data, _, _ = prepare_cluster_data(
        drep,
        gtdb,
        tax_levels=["phylum"],
        top_n=1,
        pathway_df=pathway_df,
    )

    assert cluster_data.loc[0, "representative"] == "g2_fasta"