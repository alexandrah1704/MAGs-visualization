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

def completeness_contamination_plot(checkm: pd.DataFrame, output_path: str, tag: str = "checkm", title: str | None = None):
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
    ax_scatter.axvline(70, linestyle="--", linewidth=1.2, color=col_mq, alpha=0.9)

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
    ax_histx.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)

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
    ax_histy.set_xlabel("Percentage of MAGs")
    ax_histy.set_xlim(0, (p_low_y+p_mq_y+p_hq_y).max() * 1.15)
    ax_histy.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)

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

    out = os.path.join(output_path, f"comp_conta_marginals_{tag}.png")
    fig.savefig(out, dpi=220)
    # plt.close(fig)
    # print(f"[INFO] Saved: {out}")


def rank_completeness_contamination_plot(checkm, checkm2, gtdb, rank, output_path, n):
    """
    Create two plots (for CheckM and CheckM2):
    Completeness vs Contamination, colored by GTDB <rank> (e.g., phylum),
    incl. stacked marginal histograms. Robust join via normalized IDs.
    """
    import os, re
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as mtick

    def normalize_id(s: str) -> str:
        s = str(s).strip()
        s = s.replace('.', '_').replace('-', '_')
        s = re.sub(r'\s+', '_', s)
        s = re.sub(r'(?i)_?fastq_?', '_', s)
        s = re.sub(r'(?i)_?sub_?', '_', s)
        s = re.sub(r'(?i)(\.fa(sta)?|_fasta)$', '', s)
        s = re.sub(r'_+', '_', s).strip('_')
        return s.lower()

    def _plot_rank(df_checkm, gtdb_df, rank, output_path, n, tag):
        rank_col = rank.capitalize()

        # --- CheckM: normalize key for join ---
        
        # --- DEBUG: check normalized sets BEFORE join ---
        chk_keys = pd.Index(df_checkm.index.map(normalize_id), name="key")
        gtd_keys = pd.Index(gtdb_df.index.map(normalize_id),   name="key")

        print(f"[DBG][{tag}] #checkm bins: {len(chk_keys)} | #gtdb ids: {len(gtd_keys)}")
        inter = set(chk_keys).intersection(set(gtd_keys))
        print(f"[DBG][{tag}] intersection size: {len(inter)}")

    
        # --- Join via canonical key ---
        df_checkm = df_checkm.copy()
        df_checkm["__key__"] = df_checkm.index.map(normalize_id)

        gtdb_df = gtdb_df.copy()
        if "user_genome" in gtdb_df.columns and gtdb_df.index.name != "user_genome":
            gtdb_df.set_index("user_genome", inplace=True)
        gtdb_df["__key__"] = gtdb_df.index.map(normalize_id)

        # rank extrahieren (z. B. phylum)
        rank_col = rank.capitalize()
        gtdb_df[rank_col] = gtdb_df["classification"].apply(lambda s: extract_rank(s, rank))

        left  = df_checkm.set_index("__key__", drop=True)
        right = gtdb_df.set_index("__key__", drop=True)[[rank_col]]
        merged = left.join(right, how="left")

        # --- Quick diagnostics ---
        have_rank = merged[rank_col].notna().sum()
        total_bins = len(merged)
        print(f"[DBG][{tag}] matched ranks: {have_rank}/{total_bins} ({have_rank/total_bins*100:.1f}%)")

        # --- Plot data frame ---
        x = pd.to_numeric(merged["Completeness"],  errors="coerce")
        y = pd.to_numeric(merged["Contamination"], errors="coerce")
        df = pd.DataFrame({"x": x, "y": y, rank_col: merged[rank_col]}).dropna(subset=["x", "y"])

        # --- Top-n ranks; others -> Other; NaN -> Unknown ---
        counts = df[rank_col].value_counts(dropna=False)
        top_ranks = counts.dropna().nlargest(n).index.tolist()
        other_count = counts.dropna()[~counts.dropna().index.isin(top_ranks)].sum()

        def map_rank(v):
            if pd.isna(v):
                return f"Unknown {rank_col}"
            return f"{v} ({counts[v]})" if v in top_ranks else f"Other ({other_count})"

        label_col = f"{rank_col} (n)"
        df[label_col] = df[rank_col].map(map_rank)

        # --- Layout (2x2) ---
        fig = plt.figure(figsize=(10, 8), constrained_layout=True)
        gs  = gridspec.GridSpec(nrows=2, ncols=2,
                                width_ratios=[16, 5], height_ratios=[4, 16],
                                figure=fig)
        ax_scatter = fig.add_subplot(gs[1, 0])
        ax_histx   = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
        ax_histy   = fig.add_subplot(gs[1, 1], sharey=ax_scatter)

        # --- Limits/Guides ---
        xmin, xmax = 40, 100
        ymax = max(5.0, min(7.0, np.ceil(df["y"].quantile(0.995)))) if len(df) else 7.0
        ax_scatter.set_xlim(xmin, xmax)
        ax_scatter.set_ylim(0, ymax)
        ax_scatter.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
        ax_scatter.set_xlabel("Completeness (%)")
        ax_scatter.set_ylabel("Contamination (%)")
        ax_scatter.set_title(f"{tag}: Completeness vs Contamination (Colored by {rank_col})",
                             fontsize=14, weight="bold", pad=8)
        ax_scatter.axvline(90, linestyle="--", linewidth=1.2, color="#b64a4a", alpha=0.9)
        ax_scatter.axvline(70, linestyle="--", linewidth=1.2, color="#7f7f7f", alpha=0.9)
        ax_scatter.axhline(5,  linestyle="--", linewidth=1.0, color="#bbbbbb", alpha=0.9)

        # --- Colors/Labels ---
        sns.set_theme(style="whitegrid")
        labels = list(df[label_col].value_counts().index)
        palette = sns.color_palette("tab20")
        color_map = {lab: palette[i % len(palette)] for i, lab in enumerate(labels)}

        # --- Scatter grouped by label ---
        for lab in labels:
            sub = df[df[label_col] == lab]
            ax_scatter.scatter(sub["x"], sub["y"], s=36, alpha=0.85, edgecolors="none",
                               color=color_map[lab], label=lab)

        ax_scatter.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12),
                          ncol=min(3, len(labels)), frameon=False)

        # --- Visible points for marginals ---
        vis = df[(df["x"] >= xmin) & (df["x"] <= xmax) & (df["y"] >= 0) & (df["y"] <= ymax)].copy()
        total_n = len(vis) if len(vis) > 0 else 1

        # --- Top bar (Completeness) ---
        bins_x = np.arange(xmin, xmax + 1, 1)
        centers_x = (bins_x[:-1] + bins_x[1:]) / 2
        cum_x = np.zeros(len(centers_x), dtype=float)
        for lab in labels:
            vals = vis.loc[vis[label_col] == lab, "x"]
            cnt, _ = np.histogram(vals, bins=bins_x)
            p = cnt / total_n
            ax_histx.bar(centers_x, p, width=1.0, bottom=cum_x,
                         color=color_map[lab], edgecolor="white", linewidth=0.5)
            cum_x += p
        ax_histx.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
        ax_histx.set_ylabel("% of MAGs")
        ax_histx.set_ylim(0, (cum_x.max() * 1.15) if cum_x.size else 1.0)
        ax_histx.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
        for side in ("left", "bottom"): ax_histx.spines[side].set_visible(True)
        for side in ("top", "right"):   ax_histx.spines[side].set_visible(False)
        ax_histx.tick_params(axis="both", length=4, labelleft=True, labelbottom=True)

        # --- Right bar (Contamination) ---
        bins_y = np.arange(0, ymax + 0.05, 0.1)
        centers_y = (bins_y[:-1] + bins_y[1:]) / 2
        cum_y = np.zeros(len(centers_y), dtype=float)
        for lab in labels:
            vals = vis.loc[vis[label_col] == lab, "y"]
            cnt, _ = np.histogram(vals, bins=bins_y)
            p = cnt / total_n
            ax_histy.barh(centers_y, p, height=0.1, left=cum_y,
                          color=color_map[lab], edgecolor="white", linewidth=0.5)
            cum_y += p
        ax_histy.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
        ax_histy.set_xlabel("% of MAGs")
        ax_histy.set_xlim(0, (cum_y.max() * 1.15) if cum_y.size else 1.0)
        ax_histy.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
        for side in ("left", "bottom"): ax_histy.spines[side].set_visible(True)
        for side in ("top", "right"):   ax_histy.spines[side].set_visible(False)
        ax_histy.tick_params(axis="both", length=4, labelleft=True, labelbottom=True)

        # --- Save ---
        out = os.path.join(output_path, f"comp_conta_by_rank_marginals_{tag}.png")
        fig.savefig(out, dpi=220)
        plt.close(fig)
        print(f"[INFO] Saved: {out}")

    # Run for both tools
    _plot_rank(checkm,  gtdb, rank, output_path, n, tag="CheckM")
    _plot_rank(checkm2, gtdb, rank, output_path, n, tag="CheckM2")