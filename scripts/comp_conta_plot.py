import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
import math
import numpy as np
import os                               # Dateipfade
import pandas as pd                     # Tabellen
import seaborn as sns                   # High-level Plots

prefix_map = {
    'domain': 'd',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
}

def extract_rank(classification, rank):
        prefix = prefix_map.get(rank)
        try:
            return next(x for x in classification.split(';') if x.startswith(f'{prefix}__')).replace(f'{prefix}__', '').strip()
        except StopIteration:
            return None

def completeness_contamination_plot(checkm: pd.DataFrame, output_path: str):
    """Scatter Completeness vs Contamination with marginal histograms"""
    # ---- Data ----
    df = checkm.loc[:, ["Completeness", "Contamination"]].copy()
    x = pd.to_numeric(df["Completeness"], errors="coerce")
    y = pd.to_numeric(df["Contamination"], errors="coerce")
    df = pd.DataFrame({"x": x, "y": y}).dropna()

    col_low = "#86cbd5"
    col_mq  = "#7f7f7f"
    col_hq  = "#b64a4a"

    # ---- Layout ----
    fig = plt.figure(figsize=(9, 8), constrained_layout=True)

    gs = gridspec.GridSpec(
        nrows=2, ncols=2,
        width_ratios=[16, 5],
        height_ratios=[4, 16],
        figure=fig
    )

    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_histx   = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
    ax_histy   = fig.add_subplot(gs[1, 1], sharey=ax_scatter)

    # ---- Limits / Grid ----
    # xmin = max(40, int(np.floor(df["x"].min() / 5) * 5))
    xmin = 40
    xmax = 100
    ymax = max(5.0, min(7.0, np.ceil(df["y"].quantile(0.995))))
    ax_scatter.set_xlim(xmin, xmax)
    ax_scatter.set_ylim(0, ymax)
    ax_scatter.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
    ax_scatter.set_xlabel("Completeness %")
    ax_scatter.set_ylabel("Contamination %")
    ax_scatter.axvline(90, linestyle="--", linewidth=1.2, color=col_hq, alpha=0.9)
    ax_scatter.axhline(5,  linestyle="--", linewidth=1.0, color="#bbbbbb", alpha=0.9)

    # only visible points in scatter plot - needed for histogram
    df_vis = df[(df["x"] >= xmin) & (df["x"] <= xmax) &
                (df["y"] >= 0)    & (df["y"] <= ymax)].copy()
    
    # redefined categories based on df_vis
    cat_low_vis = (df_vis["x"] < 70)
    cat_mq_vis  = (df_vis["x"] >= 70) & (df_vis["x"] < 90)
    cat_hq_vis  = (df_vis["x"] >= 90) & (df_vis["y"] <= 5)

    total_n = len(df_vis)
    
    ax_scatter.scatter(df_vis.loc[cat_low_vis, "x"], df_vis.loc[cat_low_vis, "y"], s=16, alpha=0.85, edgecolors="none", color=col_low)
    ax_scatter.scatter(df_vis.loc[cat_mq_vis,  "x"], df_vis.loc[cat_mq_vis,  "y"],  s=16, alpha=0.85, edgecolors="none", color=col_mq)
    ax_scatter.scatter(df_vis.loc[cat_hq_vis,  "x"], df_vis.loc[cat_hq_vis,  "y"],  s=18, alpha=0.95, edgecolors="none", color=col_hq)

    # ---- Histogram Completeness ----
    bins_x = np.arange(xmin, xmax + 1, 1)
    bin_centers_x = (bins_x[:-1] + bins_x[1:]) / 2

    cnt_low, _ = np.histogram(df_vis.loc[cat_low_vis, "x"], bins=bins_x)
    cnt_mq,  _ = np.histogram(df_vis.loc[cat_mq_vis,  "x"], bins=bins_x)
    cnt_hq,  _ = np.histogram(df_vis.loc[cat_hq_vis,  "x"], bins=bins_x)

    p_low = cnt_low / total_n
    p_mq  = cnt_mq  / total_n
    p_hq  = cnt_hq  / total_n

    # stacked bins
    ax_histx.bar(bin_centers_x, p_low, width=1.0, color=col_low, edgecolor="white", linewidth=0.5)
    ax_histx.bar(bin_centers_x, p_mq,  width=1.0, bottom=p_low, color=col_mq, edgecolor="white", linewidth=0.5)
    ax_histx.bar(bin_centers_x, p_hq,  width=1.0, bottom=p_low+p_mq, color=col_hq, edgecolor="white", linewidth=0.5)

    ax_histx.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax_histx.set_ylabel("Percentage of MAGs")
    ax_histx.set_ylim(0, (p_low+p_mq+p_hq).max() * 1.15)

    # --- axis lines ---
    for side in ("left", "bottom"):
        ax_histx.spines[side].set_visible(True)
    for side in ("top", "right"):
        ax_histx.spines[side].set_visible(False)

    ax_histx.tick_params(axis="both", length=4, labelleft=True, labelbottom=True)

    ax_histx.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))

    # ---- Histogram Contamination ----
    bins_y = np.arange(0, ymax + 0.05, 0.1)
    bin_centers_y = (bins_y[:-1] + bins_y[1:]) / 2

    cnt_low_y, _ = np.histogram(df_vis.loc[cat_low_vis, "y"], bins=bins_y)
    cnt_mq_y,  _ = np.histogram(df_vis.loc[cat_mq_vis,  "y"], bins=bins_y)
    cnt_hq_y,  _ = np.histogram(df_vis.loc[cat_hq_vis,  "y"], bins=bins_y)

    p_low_y = cnt_low_y / total_n
    p_mq_y  = cnt_mq_y  / total_n
    p_hq_y  = cnt_hq_y  / total_n

    # stacked horizontal bins
    ax_histy.barh(bin_centers_y, p_low_y, height=0.1, color=col_low, edgecolor="white", linewidth=0.5)
    ax_histy.barh(bin_centers_y, p_mq_y,  height=0.1, left=p_low_y, color=col_mq, edgecolor="white", linewidth=0.5)
    ax_histy.barh(bin_centers_y, p_hq_y,  height=0.1, left=p_low_y+p_mq_y, color=col_hq, edgecolor="white", linewidth=0.5)

    ax_histy.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax_histy.set_xlim(0, (p_low_y+p_mq_y+p_hq_y).max() * 1.15)

    # --- axis lines ---
    for side in ("left", "bottom"):
        ax_histy.spines[side].set_visible(True)
    for side in ("top", "right"):
        ax_histy.spines[side].set_visible(False)

    ax_histy.tick_params(axis="both", length=4, labelleft=True, labelbottom=True)

    ax_histy.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))

    from matplotlib.patches import Patch
    ax_scatter.legend(
        handles=[Patch(color=col_low, label="Partial MAGs"),
                 Patch(color=col_mq,  label="Medium-quality MAGs"),
                 Patch(color=col_hq,  label="High-quality MAGs")],
        loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False
    )

    out = os.path.join(output_path, "comp_conta_marginals.png")
    fig.savefig(out, dpi=220)
    print(f"[INFO] Saved: {out}")


def rank_completeness_contamination_plot(checkm, gtdb_bac, gtdb_ar, rank, output_path, n):

    checkm.index = checkm.index.str.replace('.', '_', regex=False)

    # connects GTDB-Archaea and Bacteria Tables
    gtdb_df = pd.concat([gtdb_ar, gtdb_bac], ignore_index=False)
    
    gtdb_df[f'{rank.capitalize()}'] = gtdb_df['classification'].apply(lambda x: extract_rank(x, rank))

    # ---- Join GTDB phylum into CheckM based on index ----
    # hängt Rang-Spalte an die CheckM-Tabelle
    merged_df = checkm.join(gtdb_df[[f'{rank.capitalize()}']])

    # ---- Create a "Phylum (n)" column for labeling ----
    rank_counts = merged_df[f'{rank.capitalize()}'].value_counts()
    top_ranks = rank_counts.nlargest(n).index.tolist()

    other_count = rank_counts[~rank_counts.index.isin(top_ranks)].sum()
    
    def map_rank(r):
        if pd.isnull(r):
            return f"Unknown {rank.capitalize()}"
        elif r in top_ranks:
            return f"{r} ({rank_counts[r]})"
        else:
            return f"Other ({other_count})"

    merged_df[f'{rank.capitalize()} (n)'] = merged_df[f'{rank.capitalize()}'].map(map_rank)

    # ---- Plot: Completeness vs Contamination with counts in legend ----
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 6))
    sns.scatterplot(
        data=merged_df,
        y='Contamination',
        x='Completeness',
        hue=f'{rank.capitalize()} (n)',
        palette='tab20',
        s=80,
        alpha=0.8
    )

    plt.title(f'CheckM: Completeness vs Contamination (Colored by {rank.capitalize()})', fontsize=14, weight='bold')
    plt.ylabel('Contamination (%)')
    plt.xlabel('Completeness (%)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=f'{rank.capitalize()} (n)')
    plt.tight_layout()

    plt.savefig(os.path.join(output_path, "comp_conta_by_rank.png"))