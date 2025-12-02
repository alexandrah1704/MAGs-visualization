import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from id_normalizer import normalize_genome_series
from heatmap import extract_rank

def drep_taxa_plot(drep_df: pd.DataFrame, gtdb_df: pd.DataFrame, 
                   output_dir: str, rank: str = "phylum", top_n: int = 15,
                   fig_size=(8, 4), fmt: str = "png",
                   ):
    """
    Two plots that shows dRep clusters by taxonomy.

    Plot 1: 
    - x-axis = tax_level
    - y-axis = number of MAGs (cluster members)

    Plot 2:
    - x-axis = number of MAGs per cluster
    - y-axis = number of clusters
    - stacked bar colored by majority taxon in each cluster
    """

    os.makedirs(output_dir, exist_ok=True)

    # ---- dRep ----
    df = drep_df.copy()

    if "genome" not in df.columns:
        if df.index.name == "genome":
            df = df.reset_index()
        else:
            raise ValueError("dRep table has no column 'genome'.")

    if "primary_cluster" not in df.columns:
        raise ValueError("dRep table has no column 'primary_cluster'.")

    df["genome_norm"] = normalize_genome_series(df["genome"])
    df = df.drop_duplicates(subset=["genome_norm", "primary_cluster"])

    cluster_sizes = (df.groupby("primary_cluster")["genome_norm"]
                     .nunique().rename("cluster_size").reset_index())

    # ---- gtdb ----
    gtdb = gtdb_df.copy()

    if "user_genome" in gtdb.columns:
        gtdb = gtdb.rename(columns={"user_genome": "genome_raw"})
    else:
        gtdb = gtdb.rename_axis("genome_raw").reset_index()

    gtdb["genome_norm"] = normalize_genome_series(gtdb["genome_raw"])

    if "classification" in gtdb.columns:
        gtdb[rank] = gtdb["classification"].apply(lambda tx: extract_rank(tx, rank))
        tax_col = rank
    else:
        tax_col = None
        for cand in [rank, rank.capitalize(), rank.lower()]:
            if cand in gtdb.columns:
                tax_col = cand
                break
        if tax_col is None:
            raise ValueError(
                f"GTDB has no column 'classification' and no column '{rank}'."
            )

    gtdb = gtdb[["genome_norm", tax_col]].copy()

    df = df.merge(gtdb, on="genome_norm", how="left")
    df[tax_col] = df[tax_col].fillna("Unclassified")

    # ---- Plot 1 ----
    sns.set_style("whitegrid")

    genome_counts = (
        df.groupby(tax_col)["genome_norm"]
        .nunique()
        .reset_index()
        .rename(columns={tax_col: "taxon", "genome_norm": "n_genomes"})
        .sort_values("n_genomes", ascending=False)
    )

    # Top_n Taxa + "Others"
    if len(genome_counts) > top_n:
        top_taxa = genome_counts.head(top_n)
        others_count = genome_counts.iloc[top_n:]["n_genomes"].sum()
        others_row = pd.DataFrame([{"taxon": "Others", "n_genomes": others_count}])
        genome_counts = pd.concat([top_taxa, others_row], ignore_index=True)

    fig1, ax1 = plt.subplots(figsize=fig_size)

    pal = sns.color_palette("tab20", n_colors=len(genome_counts))
    colors = pal[: len(genome_counts)]

    bars = ax1.bar(
        range(len(genome_counts)),
        genome_counts["n_genomes"],
        color=colors,
        edgecolor="black",
        linewidth=0.3,
        alpha=0.85,
    )

    for bar, count in zip(bars, genome_counts["n_genomes"]):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            str(count),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax1.set_xlabel(rank.capitalize(), fontsize=12)
    ax1.set_ylabel("Number of genomes (cluster members)", fontsize=12)
    ax1.set_title(f"Total cluster members by {rank}", fontweight="bold", fontsize=14)

    ax1.set_xticks(range(len(genome_counts)))
    ax1.set_xticklabels(
        genome_counts["taxon"], rotation=45, ha="right", fontsize=10
    )

    ax1.grid(axis="y", linestyle="--", alpha=0.3)
    ax1.grid(axis="x", linestyle="--", alpha=0.3)
    ax1.set_axisbelow(True)

    plt.tight_layout()
    safe_rank = rank.lower().replace(" ", "_")
    out_path1 = os.path.join(output_dir, f"drep_total_members_by_{safe_rank}.{fmt}")
    plt.savefig(out_path1, dpi=300, bbox_inches="tight")
    plt.close(fig1)
    print(f"[OK] Saved Plot 1 (Total members): {out_path1}")

    # ---- Plot 2 ----

    # compute majority of a cluster
    cluster_diversity = (df.groupby("primary_cluster")[tax_col]
                         .agg(lambda x: x.nunique())
                         .reset_index()
                         .rename(columns={tax_col: "n_taxa"}))

    # Warning for heterogen clusters
    # heterogeneous = cluster_diversity[cluster_diversity["n_taxa"] > 1]
    # for cluster_id in heterogeneous["primary_cluster"].head(5):  # Top 5 anschauen
    #     cluster_genomes = df[df["primary_cluster"] == cluster_id]
    #     print(f"\n=== Cluster: {cluster_id} ===")
    #     print(cluster_genomes[[tax_col, "genome_norm"]].value_counts())
    
    # majority taxon
    cluster_taxa = (
        df.groupby("primary_cluster")[tax_col]
        .agg([
            lambda x: x.mode()[0] if len(x.mode()) > 0 else "Unclassified",  # Mehrheit
            lambda x: (x == x.mode()[0]).sum() / len(x) if len(x.mode()) > 0 else 0  # Anteil
        ])
        .reset_index()
    )
    cluster_taxa.columns = ["primary_cluster", "cluster_taxon", "consensus_fraction"]
    cluster_data = cluster_sizes.merge(cluster_taxa, on="primary_cluster")

    bins = [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 100, np.inf]
    labels = ["1", "2", "3-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-49", "50-99", "100+"]

    cluster_data["size_bin"] = pd.cut(cluster_data["cluster_size"], bins=bins,
                                      labels=labels, right=False, include_lowest=True)

    top_taxa_names = genome_counts["taxon"].tolist()
    if "Others" in top_taxa_names:
        top_taxa_names.remove("Others")
        
    cluster_data["cluster_taxon_grouped"] = cluster_data["cluster_taxon"].apply(
        lambda x: x if x in top_taxa_names else "Others"
    )

    # pivot for stacked bar
    pivot_data = (
        cluster_data.groupby(["size_bin", "cluster_taxon_grouped"])
        .size()
        .unstack(fill_value=0)
    )

    # Order tax
    taxa_order = [t for t in top_taxa_names if t in pivot_data.columns]
    if "Others" in pivot_data.columns:
        taxa_order.append("Others")
    pivot_data = pivot_data[taxa_order]

    fig2, ax2 = plt.subplots(figsize=fig_size)

    pivot_data.plot(
        kind="bar",
        stacked=True,
        ax=ax2,
        color=colors[:len(taxa_order)],
        edgecolor="black",
        linewidth=0.3,
        alpha=0.85,
        width=0.8
    )

    ax2.set_xlabel("Cluster size (number of genomes)", fontsize=12)
    ax2.set_ylabel("Number of clusters", fontsize=12)
    ax2.set_title(
        f"Distribution of cluster sizes by {rank}",
        fontweight="bold",
        fontsize=14
    )

    ax2.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
    ax2.legend(title=rank.capitalize(), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax2.grid(axis="y", linestyle="--", alpha=0.3)
    ax2.set_axisbelow(True)

    plt.tight_layout()
    out_path2 = os.path.join(output_dir, f"drep_cluster_sizes_by_{safe_rank}.{fmt}")
    plt.savefig(out_path2, dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print(f"[OK] Saved Plot 2 (Cluster sizes): {out_path2}")