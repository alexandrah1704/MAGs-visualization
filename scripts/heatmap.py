import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle
from id_normalizer import normalize_genome_id


def extract_rank(tax, rank: str):
    """Extract a GTDB rank from a taxonomy string."""
    if pd.isna(tax):
        return None
    prefix_map = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }
    prefix = prefix_map.get(rank.lower())
    if not prefix:
        raise ValueError(f"Invalid rank '{rank}'. Choose one of {list(prefix_map.keys())}.")
    for part in str(tax).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            name = part[len(prefix):].strip()
            return name if name else "Unnamed"
    return "Unnamed"

def normalize_id(s: str) -> str:
    return normalize_genome_id(s)


def clean_sample_from_coverm_col(col: str) -> str:
    """
    Extract a compact sample label from a CoverM column like:
    """
    m = re.search(r"(SRR\d+)", col)
    if m:
        return m.group(1)
    return col.split()[0].replace(".fastq", "")


def mag_heatmap(coverm_df: pd.DataFrame, gtdb_df: pd.DataFrame, output_path: str,
    present_threshold: float = 0.0, top_bar_spacing: float = 1.0,
    top_bar_width: float = 0.90, rank: str = "phylum",
    metadata_df: pd.DataFrame | None = None, meta_col: str | None = None,
    meta_bin_width: float = 5.0, fmt: str = "png"):

    """
    Combined visualization:
    - top: log10(MAGs/Phylum)
    - center: Heatmap showing relative abundance
    - right: MAGs/sample
    - left: legend Abundance, legend weather
    """

    os.makedirs(output_path, exist_ok=True)

    # ---- GTDB - selected rank ----
    gtdb = gtdb_df.copy()
    if "classification" not in gtdb.columns:
        raise ValueError("GTDB table must contain column 'classification'.")

    if gtdb.index.name != "user_genome":
        gtdb = gtdb.set_index("user_genome")
    gtdb.index = gtdb.index.map(normalize_id)

    gtdb[rank] = gtdb["classification"].apply(lambda tax: extract_rank(tax, rank))
    gtdb = gtdb[[rank]].dropna()

    # ---- CoverM ----
    cov = coverm_df.copy()
    cov = cov.loc[:, ~cov.columns.duplicated(keep="first")]

    if "Genome" not in cov.columns:
        cov = cov.reset_index().rename(columns={cov.columns[0]: "Genome"})

    cov = cov.loc[:, ~cov.columns.duplicated(keep="first")]

    sample_cols_raw = [c for c in cov.columns if "Relative Abundance (%)" in c]

    cov = cov[["Genome"] + sample_cols_raw].copy()
    cov["Genome"] = cov["Genome"].astype(str)
    cov = cov[cov["Genome"].str.lower() != "unmapped"].copy()

    # clean sample labels
    clean_map = {c: clean_sample_from_coverm_col(c) for c in sample_cols_raw}
    cov.rename(columns=clean_map, inplace=True)
    for c in clean_map.values():
        cov[c] = pd.to_numeric(cov[c], errors="coerce").fillna(0.0)
    cov["user_genome"] = cov["Genome"].map(normalize_id)

    # Join MAG
    merged = cov.merge(gtdb, left_on="user_genome", right_index=True, how="left").dropna(subset=[rank])
    value_cols = [c for c in cov.columns if c not in ["Genome", "user_genome", rank]]
    
    # ---- Long form for heatmap ----
    long_df = merged.melt(id_vars=[rank], value_vars=[c for c in merged.columns if c in value_cols],
                          var_name="sample", value_name="abundance")

    # ---- Heat matrix: Sample × Rank ----
    heat = (long_df.groupby(["sample", rank], as_index=False)["abundance"].sum()
            .pivot(index="sample", columns=rank, values="abundance")
            .fillna(0.0))

    # sort taxa columns by total abundance;'Unnamed' to the end
    cols = heat.sum(axis=0).sort_values(ascending=False).index.tolist()
    if "Unnamed" in cols:
        cols = [c for c in cols if c != "Unnamed"] + ["Unnamed"]
    heat = heat.loc[:, cols]

    n_rows, n_cols = heat.shape

    # ---- Top bar ----
    mags_per_rank = (merged.groupby(rank)["Genome"]
                     .nunique()
                     .reindex(heat.columns)
                     .fillna(0).astype(int))
    top_vals = mags_per_rank.replace(0, np.nan).apply(np.log10)
    if top_vals.isna().all():
        top_vals = pd.Series([0] * len(mags_per_rank), index=mags_per_rank.index)

    # ---- Right bar ----
    cov_for_count = cov.set_index("Genome")[
        [c for c in cov.columns if c not in ["Genome", "user_genome", rank]]]
    mags_per_sample = (
        (cov_for_count > present_threshold).sum(axis=0).reindex(heat.index).fillna(0).astype(int))

    # ---- Metadata per sample ----
    meta_category_per_sample = None
    meta_color_per_category = None
    meta_title = None

    if metadata_df is not None and meta_col is not None:
        if meta_col not in metadata_df.columns:
            print(f"[WARN] meta_col '{meta_col}' not found in heatmap metadata. Available columns: {list(metadata_df.columns)}")
        else:
            # align metadata to sample names used in heatmap
            meta_raw = metadata_df[meta_col]
            meta_aligned = meta_raw.reindex(heat.index)

            print(f"[INFO] Using heatmap metadata '{meta_col}': {meta_aligned.notna().sum()}/{len(meta_aligned)} non-NA after alignment")

            if meta_aligned.notna().sum() == 0:
                print("[WARN] All values are NA after aligning heatmap metadata to samples. No metadata colors will be shown.")
            else:
                unknown_label = f"Unknown {meta_col}"
                meta_num = pd.to_numeric(meta_aligned, errors="coerce")
                frac_numeric = meta_num.notna().mean()

                if frac_numeric >= 0.8 and meta_num.notna().sum() > 0:
                    # numeric (temperature)
                    vmin = float(np.nanmin(meta_num.values))
                    vmax = float(np.nanmax(meta_num.values))
                    if vmin == vmax:
                        vmax = vmin + meta_bin_width

                    bins = np.arange(
                        np.floor(vmin / meta_bin_width) * meta_bin_width,
                        np.ceil(vmax / meta_bin_width) * meta_bin_width + meta_bin_width,
                        meta_bin_width
                    )
                    labels_bins = [
                        f"{int(left)}–{int(right)}"
                        for left, right in zip(bins[:-1], bins[1:])
                    ]

                    meta_binned_arr = pd.cut(
                        meta_num.to_numpy(),
                        bins=bins,
                        labels=labels_bins,
                        include_lowest=True,
                    )

                    meta_binned = pd.Series(
                        meta_binned_arr,
                        index=meta_num.index,
                        dtype="object",
                    )
                    meta_binned = meta_binned.where(meta_binned.notna(), other=unknown_label)

                    meta_category_per_sample = meta_binned.to_dict()
                    meta_title = f"{meta_col} (binned, {meta_bin_width:g})"
                else:
                    # categorical (weather)
                    meta_cat = (
                        meta_aligned.astype("string")
                        .fillna(unknown_label)
                        .replace("", unknown_label)
                    )
                    meta_category_per_sample = meta_cat.to_dict()
                    meta_title = meta_col

                categories = list(sorted(set(meta_category_per_sample.values())))
                base_palette = [
                    "#60a5fa", "#f97316", "#22c55e", "#8b5cf6", "#ec4899",
                    "#e5e7eb", "#10b981", "#facc15", "#64748b", "#fb7185",
                    "#0ea5e9", "#a3e635", "#fbbf24", "#f97373", "#4b5563",
                ]
                meta_color_per_category = {
                    cat: base_palette[i % len(base_palette)]
                    for i, cat in enumerate(categories)
                }

    # ---- colormap and bins ----
    boundaries = [0, 1, 2, 4, 8, 16, 40, 60, 80, 1000]
    colors = [
        "#ffffff", "#e2f5e8", "#bfe6c9", "#88d0a6", "#48b07c",
        "#219c6a", "#ffb67a", "#e0554a", "#7f1d1d"
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    bin_labels = ["0", "1–2", "2–4", "4–8", "8–16", "16–40", "40–60", "60–80", ">80"]

    # ---- layout ----
    safe_cols = max(1, n_cols)
    safe_rows = max(1, n_rows)
    fig = plt.figure(figsize=(max(10, safe_cols * 0.62), max(8, safe_rows * 0.32)))

    gs = gridspec.GridSpec(
        3, 5, figure=fig,
        height_ratios=[1.0, 0.0, 8.0],             # Top, spacer, heatmap row
        width_ratios=[3.0, 0.5, 8.0, 0.2, 2.0],    # Legends | spacer | Heatmap | spacer | Right bars
        wspace=0.05, hspace=0.18
    )

    ax_top   = fig.add_subplot(gs[0, 2])
    ax_heat  = fig.add_subplot(gs[2, 2])
    ax_right = fig.add_subplot(gs[2, 4])

    # left column split into abundance legend & weather legend
    left = gs[2, 0].subgridspec(2, 1, height_ratios=[3, 2], hspace=0.35)
    ax_abund_legend  = fig.add_subplot(left[0, 0])
    ax_meta_legend = fig.add_subplot(left[1, 0])

    # ---- heatmap ----
    im = ax_heat.imshow(heat.values, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)

    ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heat.grid(which="minor", color="#d0d0d0", linewidth=0.5)
    ax_heat.grid(which="major", visible=False)

    ax_heat.set_xlim(-0.5, safe_cols - 0.5)
    ax_heat.set_ylim(safe_rows - 0.5, -0.5)
    ax_heat.set_xticks(np.arange(n_cols))
    ax_heat.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=9)
    ax_heat.set_yticks(np.arange(n_rows))
    ax_heat.set_yticklabels(heat.index, rotation=0, fontsize=9)
    ax_heat.set_xlabel(rank.capitalize())
    ax_heat.set_ylabel("")

    # ---- left: abundance legend ----
    ax_abund_legend.set_title("Abundance", fontsize=10, pad=6)
    ax_abund_legend.title.set_x(0.15)
    n_bins = len(bin_labels)
    ax_abund_legend.set_xlim(0, 1)
    ax_abund_legend.set_ylim(0, n_bins)
    for i, (label, color) in enumerate(zip(reversed(bin_labels), reversed(colors))):
        y = i
        ax_abund_legend.add_patch(
            Rectangle((0.05, y), 0.25, 1.0, color=color, ec="#888888", linewidth=0.5)
        )
        ax_abund_legend.text(0.35, y + 0.5, label, va="center", fontsize=9)
    ax_abund_legend.axis("off")

    # ---- Left: metadata legend (if available) ----
    if meta_category_per_sample is not None and meta_color_per_category is not None and meta_title is not None:
        ax_meta_legend.set_title(meta_title, fontsize=10, pad=6)
        ax_meta_legend.title.set_x(0.15)

        cats = list(meta_color_per_category.keys())
        ax_meta_legend.set_xlim(0, 1)
        ax_meta_legend.set_ylim(0, len(cats))

        for i, cat in enumerate(reversed(cats)):
            y = i
            color = meta_color_per_category[cat]
            ax_meta_legend.add_patch(
                Rectangle((0.05, y), 0.25, 1.0, color=color, ec="#888888", linewidth=0.5)
            )
            ax_meta_legend.text(0.35, y + 0.5, str(cat), va="center", fontsize=9)
        ax_meta_legend.axis("off")
    else:
        # no metadata, no second legend
        ax_meta_legend.axis("off")
    
    # ---- top bar ----
    x_pos = np.linspace(0, n_cols - 1, n_cols) * top_bar_spacing
    ax_top.bar(x_pos, top_vals.values, color="#6b6b6b", edgecolor="#444444",
               width=top_bar_width, align="center")
    ax_top.set_xlim(-0.5, n_cols - 0.5)
    ax_top.set_ylabel("log$_{10}$(MAGs/Rank)")
    ax_top.set_xticks([])

    # ---- Dynamic grid/ticks setup ----
    if np.isfinite(top_vals).any():
        ymax = float(np.nanmax(top_vals.values))
        yceil = int(np.ceil(ymax))
        ax_top.set_ylim(0, max(1, yceil))
        
        yticks = np.arange(0, yceil + 1, 1)
        ax_top.set_yticks(yticks)
        
        # Draw horizontal reference line at 0 only
        ax_top.axhline(0, color="#888888", linewidth=0.8)
        if yceil >= 2:
            ax_top.axhline(2, color="#cccccc", linewidth=0.5, linestyle="--")
    else:
        ax_top.set_ylim(0, 1)
        ax_top.set_yticks([0, 1])

    # Remove visual clutter
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.grid(axis='y', color='#888888', linewidth=0.5, zorder=10)
    ax_top.set_axisbelow(False)

    # ---- right bar ----
    ypos = np.arange(n_rows)

    if meta_category_per_sample is not None and meta_color_per_category is not None:
        bar_colors = [
            meta_color_per_category.get(meta_category_per_sample.get(sample, None), "#668dcc")
            for sample in heat.index
        ]
    else:
        bar_colors = "#668dcc"

    ax_right.barh(
        ypos,
        mags_per_sample.values,
        height=0.7,
        facecolor=bar_colors,
        edgecolor="#1f2937",
        linewidth=0.6,
        zorder=1
    )

    max_val = int(mags_per_sample.max())
    nice_max = max(10, int(np.ceil(max_val / 50.0) * 50))
    half = nice_max // 2

    ax_right.set_xlim(0, nice_max)
    ax_right.set_ylim(n_rows - 0.5, -0.5)

    # 50 step ticks
    ticks = np.arange(0, nice_max + 1, 50)
    ax_right.set_xticks(ticks)
    ax_right.set_xticklabels([str(int(t)) for t in ticks])

    # reference lines above bars
    for x in (half, nice_max):
        ax_right.axvline(x, color="#6b7280", linewidth=1.0, zorder=10)

    ax_right.grid(axis='x', color="#717172", linewidth=0.2, zorder=10)
    ax_right.set_axisbelow(False)
    ax_right.set_title("MAGs/sample", pad=6, fontsize=10)
    ax_right.set_yticks([])
    ax_right.tick_params(axis="x", labelsize=9)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    ax_right.spines["bottom"].set_visible(False)

    plt.suptitle(f"MAG distribution: samples × {rank}", y=0.98, fontsize=12)
    out_png = os.path.join(output_path, f"heatmap_with_bars.{fmt}")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    return out_png