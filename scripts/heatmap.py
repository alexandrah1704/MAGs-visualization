import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
import os


def extract_phylum(tax):
    if pd.isna(tax):
        return None
    for part in str(tax).split(";"):
        if part.startswith("p__"):
            return part.replace("p__", "")
    return None

def normalize_id(id_str: str) -> str:
    """Normalize genome IDs to match gdtb IDs"""
    s = id_str.replace(".", "_")
    if not (s.endswith("_fasta") or s.endswith(".fasta")):
        s = s + "_fasta"
    return s

def clean_sample_label(s: str) -> str:
    return s.split()[0].replace(".fastq", "")

def mag_heatmap(coverm_df: pd.DataFrame, gtdb_df: pd.DataFrame, output_path: str,
                present_threshold: float = 0.0,
                top_bar_spacing: float = 0.95,
                top_bar_width: float = 0.90):
    """
    Combined visualization:
    - top: log10(MAGs/Phylum)
    - center: Heatmap showing relative abundance
    - right: MAGs/sample
    """
    
    # GTDB: Phylum-Column
    gtdb = gtdb_df.copy()
    if "classification" not in gtdb.columns:
        raise ValueError("GTDB-DataFrame has no column 'classification'.")
    
    gtdb["phylum"] = gtdb["classification"].apply(extract_phylum)
    gtdb = gtdb[["phylum"]].dropna()
    gtdb.index.name = "user_genome"

    # CoverM: Index->Column, map MAGs to gtdb taxonomy
    cov = coverm_df.copy()
    cov = cov.reset_index().rename(columns={cov.index.name or "index": "Genome", cov.columns[0]: cov.columns[0]})
    cov["user_genome"] = cov["Genome"].apply(normalize_id)

    # Join MAGs to their gtdb phylum assignments
    merged = cov.merge(gtdb, left_on="user_genome", right_index=True, how="left").dropna(subset=["phylum"])

    # Coverm table to long format (sample, abundance, phylum)
    id_cols = ["Genome", "user_genome", "phylum"]
    value_cols = [c for c in merged.columns if c not in id_cols]

    long_df = merged.melt(id_vars=["phylum"], value_vars=value_cols,
                          var_name="sample", value_name="abundance")

    all_samples = value_cols[:]
    clean_map = {s: clean_sample_label(s) for s in all_samples}

    # ---- Heatmap-Matrix: Sample × Phylum ----
    heat = (long_df.groupby(["sample", "phylum"], as_index=False)["abundance"].sum()
                    .pivot(index="sample", columns="phylum", values="abundance")
                    .fillna(0.0))

    heat.index = heat.index.map(lambda s: clean_map.get(s, clean_sample_label(s)))

    # sort phyla by total abundance
    heat = heat.loc[:, heat.sum(axis=0).sort_values(ascending=False).index]
    n_rows, n_cols = heat.shape

    # ---- Top bar chart ----
    mags_per_phylum = (merged.groupby("phylum")["Genome"]
                       .nunique()
                       .reindex(heat.columns)
                       .fillna(0).astype(int))
    top_vals = pd.Series(np.log10(mags_per_phylum.replace(0, np.nan)),
                         index=mags_per_phylum.index)
    # top_vals = mags_per_phylum # wenn ohne log10 


    # ---- Right bar chart ----
    cov_mag_sample = coverm_df.copy()
    cov_mag_sample.rename(columns=clean_map, inplace=True)
    mags_per_sample = (cov_mag_sample > present_threshold).sum(axis=0)
    mags_per_sample = mags_per_sample.reindex(heat.index).fillna(0).astype(int)

    # ---- Colormap abundance ----
    boundaries = [0, 1, 2, 4, 8, 16, 40, 60, 80, 1000]
    colors = [
        "#ffffff", "#e2f5e8", "#bfe6c9", "#88d0a6", "#48b07c",
        "#219c6a", "#ffb67a", "#e0554a", "#7f1d1d"
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    bin_labels = ["0", "1–2", "2–4", "4–8", "8–16", "16–40", "40–60", "60–80", ">80"]

    # ---- Layout ----
    fig = plt.figure(figsize=(max(10, n_cols * 0.6), max(8, n_rows * 0.3)))
    gs = gridspec.GridSpec(
        3, 3, figure=fig,
        height_ratios=[1.0, 0.2, 8.0],     # Top bar, space, heatmap
        width_ratios=[8.0, 0.5, 2.0],      # Heatmap, space, right
        wspace=0.20, hspace=0.18
    )

    ax_top   = fig.add_subplot(gs[0, 0])
    ax_heat  = fig.add_subplot(gs[2, 0])
    ax_right = fig.add_subplot(gs[2, 2])

    # ---- Heatmap + Grid ----
    im = ax_heat.imshow(heat.values, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)

    # Grid
    ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heat.grid(which="minor", color="#d0d0d0", linewidth=0.5)
    ax_heat.tick_params(which="both", length=0)

    # Ticks/Labels
    ax_heat.set_xlim(-0.5, n_cols - 0.5)
    ax_heat.set_xticks(np.arange(n_cols))
    ax_heat.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=9)
    ax_heat.set_yticks(np.arange(n_rows))
    ax_heat.set_yticklabels(heat.index, rotation=0, fontsize=9)
    ax_heat.set_xlabel("Phylum")
    ax_heat.set_ylabel("Samples")

    # Colorbar abundance
    cbar = fig.colorbar(
        im, ax=ax_heat, fraction=0.03, pad=0.02,
        ticks=[(boundaries[i]+boundaries[i+1])/2 for i in range(len(boundaries)-1)]
    )
    cbar.ax.set_yticklabels(bin_labels)
    cbar.set_label("Relative abundance (%)", rotation=90)

    # ---- Top bar ----
    x_pos = np.linspace(0, n_cols - 1, n_cols) * top_bar_spacing
    ax_top.bar(x_pos, top_vals.values, color="#6b6b6b", edgecolor="#444444",
               width=top_bar_width, align="center")
    ax_top.set_xlim(-0.5, n_cols - 0.5)  # gleiche Breite wie Heatmap
    ax_top.set_ylabel("log$_{10}$(MAGs/Phylum)")
    ax_top.set_xticks([])
    ax_top.axhline(0, color="#888888", linewidth=0.8)

    # ---- Right bar ----
    ax_right.barh(np.arange(n_rows), mags_per_sample.values, edgecolor="#444444", height=0.9)
    ax_right.set_ylim(-0.5, n_rows - 0.5)
    ax_right.set_xlabel("MAGs/sample")
    ax_right.set_yticks([])
    ax_right.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.6)

    plt.suptitle("MAG distribution: samples × phyla", y=0.98, fontsize=12)
    os.makedirs(output_path, exist_ok=True)
    out_png = os.path.join(output_path, "heatmap_with_bars.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    return out_png