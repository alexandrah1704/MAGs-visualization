# Heatmap visualization for MAG:
# - Taxonomic distribution of MAGs
# - Relative abundance visualization
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle
from id_normalizer import normalize_genome_id
from pandas.api.types import is_bool_dtype

# ---- Tax and ID processing ----
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
            return name if name else "Unclassified"
    return "Unclassified"

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

# ---- Metadata ----
def process_metadata_column(meta_aligned, meta_col, meta_bin_width=5.0, palette_index=0):
    """
    Process a single metadata column and return category and color mapping
    """
    # Different color palettes for each metadata column
    palettes = [
        # Palette 0: Blues/Oranges/Greens
        ["#60a5fa", "#f97316", "#22c55e", "#8b5cf6", "#ec4899",
         "#e5e7eb", "#10b981", "#facc15", "#64748b", "#fb7185",
         "#0ea5e9", "#a3e635", "#fbbf24", "#f97373", "#4b5563"],
        # Palette 1: reds/yellows/browns
        ["#fbbf24", "#d97706", "#ea580c", "#dc2626", "#b91c1c",
         "#78350f", "#f59e0b", "#ca8a04", "#92400e", "#f87171",
         "#fb923c", "#fde047", "#ef4444", "#c2410c", "#ea580c",],
        # Palette 2: purples/pinks/teals
        ["#8b5cf6", "#d946ef", "#ec4899", "#06b6d4", "#0891b2",
         "#a855f7", "#e879f9", "#f472b6", "#14b8a6", "#0d9488",
         "#c026d3", "#a21caf", "#be185d", "#6366f1", "#4f46e5"],
        # Palette 3: greens/browns/olives
        ["#16a34a", "#65a30d", "#84cc16", "#a3e635", "#15803d",
         "#4d7c0f", "#166534", "#14532d", "#854d0e", "#713f12",
         "#22c55e", "#86efac", "#bef264", "#365314", "#422006"],
        # Palette 4: Pastels
        ["#bfdbfe", "#fed7aa", "#bbf7d0", "#ddd6fe", "#fbcfe8",
         "#fde68a", "#a5f3fc", "#fecaca", "#c7d2fe", "#fae8ff",
         "#99f6e4", "#fef3c7", "#dbeafe", "#e9d5ff", "#fce7f3"],
    ]
    
    base_palette = palettes[palette_index % len(palettes)]
    
    unknown_label = f"Unknown {meta_col}"

    raw = meta_aligned.dropna()

    bool_like = False    
    meta_num = pd.to_numeric(meta_aligned, errors="coerce")
    frac_numeric = meta_num.notna().mean()

    if is_bool_dtype(raw):
        bool_like = True
    else:
        lowered = {str(v).strip().lower() for v in pd.unique(raw)}
        bool_tokens = {"true", "false", "yes", "no", "0", "1"}
        if lowered and lowered.issubset(bool_tokens):
            bool_like = True

    if bool_like:
        meta_cat = (
            meta_aligned.astype("string")
            .fillna(unknown_label)
            .replace("", unknown_label)
        )
        meta_category_per_sample = meta_cat.to_dict()
        meta_title = meta_col

        categories = list(sorted(set(meta_category_per_sample.values())))
        meta_color_per_category = {
            cat: base_palette[i % len(base_palette)]
            for i, cat in enumerate(categories)
        }
        return meta_category_per_sample, meta_color_per_category, meta_title


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
    meta_color_per_category = {
        cat: base_palette[i % len(base_palette)]
        for i, cat in enumerate(categories)
    }

    return meta_category_per_sample, meta_color_per_category, meta_title


def mag_heatmap(coverm_df: pd.DataFrame, gtdb_df: pd.DataFrame, output_path: str,
    present_threshold: float = 0.0, top_bar_spacing: float = 1.0,
    top_bar_width: float = 0.90, rank: str = "phylum",
    metadata_df: pd.DataFrame | None = None, meta_cols: str | None = None,
    meta_bin_width: float = 5.0, fmt: str = "png",
    # Layout parameters
    top_bar_height: float = 0.7,
    hspace: float = 0.25,
    heatmap_width: float = 11.0,
    spacer_legend: float = 0.3,
    spacer_meta: float = 2.0,
    spacer_heatmap: float = 0.10,
    legend: float = 3.5,
    meta_bar_add: float = 1.5,
    top_bar_spacer: float = 0.0,
    max_col: int = 10,
    log_top: bool = True,
    ):

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
        if cov.index.name == "Genome":
            cov = cov.reset_index()
        else:
            cov = cov.reset_index()
            old_name = cov.columns[0]
            if old_name != "Genome":
                cov = cov.rename(columns={old_name: "Genome"})
    else:
        pass

    cov = cov.loc[:, ~cov.columns.duplicated(keep="first")]

    sample_cols_raw = [c for c in cov.columns if c not in ["Genome", "user_genome"]]

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

    sample_labels = [c for c in cov.columns if c not in ["Genome", "user_genome"]]
    heat = heat.reindex(index=sample_labels, fill_value=0.0)

    # sort taxa columns by total abundance
    cols = heat.sum(axis=0).sort_values(ascending=False).index.tolist()

    # ensure Unclassified is always last (optional)
    has_unclassified = "Unclassified" in cols
    classified = [c for c in cols if c != "Unclassified"]
    unclassified = ["Unclassified"] if has_unclassified else []

    # ---- Top-k + Others ----
    if max_col is not None and max_col > 0 and len(classified) > max_col:
        top = classified[:max_col]
        rest = classified[max_col:]

        heat["Others"] = heat[rest].sum(axis=1)
        heat = heat.drop(columns=rest)

        # column order: top + Others + Unclassified
        cols = top + ["Others"] + unclassified
    else:
        cols = classified + unclassified

    heat = heat.loc[:, cols]

    n_rows, n_cols = heat.shape

    # ---- Top bar ----
    mags_per_rank = (merged.groupby(rank)["Genome"]
                     .nunique()
                     .reindex(heat.columns)
                     .fillna(0).astype(int))
    
    # If we created "Others", compute its MAG count too
    if "Others" in heat.columns:
        top_set = set([c for c in heat.columns if c not in ("Others", "Unclassified")])

        # all ranks present in merged excluding selected top + Unclassified
        merged_ranks = merged[rank].astype(str)
        other_mask = ~merged_ranks.isin(top_set) & (merged_ranks != "Unclassified")

        mags_per_rank.loc["Others"] = merged.loc[other_mask, "Genome"].nunique()

    if log_top:
        top_vals = mags_per_rank.replace(0, np.nan).apply(np.log10)
        y_label = "log$_{10}$(MAGs / rank)"
    else:
        top_vals = mags_per_rank.astype(float)
        y_label = "MAGs / rank"

    if top_vals.isna().all():
        top_vals = pd.Series([0] * len(mags_per_rank), index=mags_per_rank.index)

    # ---- Right bar ----
    cov_for_count = cov.set_index("Genome")[
        [c for c in cov.columns if c not in ["Genome", "user_genome", rank]]]
    mags_per_sample = (
        (cov_for_count > present_threshold).sum(axis=0).reindex(heat.index).fillna(0).astype(int))

    # ---- Process all metadata columns ----
    all_metadata = []
    
    if metadata_df is not None and meta_cols is not None:
        for idx, meta_col in enumerate(meta_cols):
            if meta_col not in metadata_df.columns:
                print(f"[WARN] meta_col '{meta_col}' not found in heatmap metadata. Available columns: {list(metadata_df.columns)}")
                continue
            
            meta_raw = metadata_df[meta_col]
            meta_aligned = meta_raw.reindex(heat.index)

            print(f"[INFO] Using heatmap metadata '{meta_col}': {meta_aligned.notna().sum()}/{len(meta_aligned)} non-NA after alignment")

            if meta_aligned.notna().sum() == 0:
                print(f"[WARN] All values are NA after aligning heatmap metadata to samples for '{meta_col}'. Skipping.")
                continue
            
            category_map, color_map, title = process_metadata_column(meta_aligned, meta_col, meta_bin_width, palette_index=idx)
            all_metadata.append({
                'category_map': category_map,
                'color_map': color_map,
                'title': title
            })

    # First metadata column is used for right bar coloring
    first_meta = all_metadata[0] if all_metadata else None

    # ---- colormap and bins Abundance ----
    boundaries = [0, 1, 2, 4, 8, 16, 40, 60, 80, 1000]
    colors = [
        "#ffffff", "#e2f5e8", "#bfe6c9", "#88d0a6", "#48b07c",
        "#219c6a", "#ffb67a", "#e0554a", "#7f1d1d"
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    bin_labels = ["0-1", "1–2", "2–4", "4–8", "8–16", "16–40", "40–60", "60–80", ">80"]

    # ---- layout ----
    safe_cols = max(1, n_cols)
    safe_rows = max(1, n_rows)

    # Calculate width ratios: legends | spacer | metadata bars | spacer | heatmap | spacer | right bars
    n_meta_cols = len(all_metadata)
    meta_bars_width = n_meta_cols * 0.50 if n_meta_cols > 0 else 0.0
    width_ratios = [legend, spacer_legend, meta_bars_width + meta_bar_add, spacer_meta, heatmap_width, spacer_heatmap, 2.5] if n_meta_cols > 0 else [3.5, 0.6, 8.0, 0.3, 2.5]
    
    fig = plt.figure(figsize=(max(12, safe_cols * 0.80 + meta_bars_width), max(7, safe_rows * 0.35 + 1.5)))

    if n_meta_cols > 0:
        gs = gridspec.GridSpec(
            3, 7, figure=fig,
            height_ratios=[top_bar_height, top_bar_spacer, 8.0],
            width_ratios=width_ratios,
            wspace=0.10, hspace=hspace
        )
        ax_top = fig.add_subplot(gs[0, 4])
        ax_heat = fig.add_subplot(gs[2, 4])
        ax_right = fig.add_subplot(gs[2, 6])
        ax_meta_bars = fig.add_subplot(gs[2, 2])
    else:
        gs = gridspec.GridSpec(
            3, 5, figure=fig,
            height_ratios=[top_bar_height, top_bar_spacer, 8.0],
            width_ratios=width_ratios,
            wspace=0.10, hspace=hspace
        )
        ax_top = fig.add_subplot(gs[0, 2])
        ax_heat = fig.add_subplot(gs[2, 2])
        ax_right = fig.add_subplot(gs[2, 4])
        ax_meta_bars = None


    # left column split into abundance legend & metadata legend
    n_meta_legends = len(all_metadata) if all_metadata else 0

    abund_height = 2.0 if n_meta_legends == 0 else 4.0
    meta_heights = []

    for meta_info in all_metadata:
        num_cats = len(meta_info["color_map"])
        meta_heights.append( max(3.0, num_cats * 0.25) )

    left_height_ratios = [abund_height] + meta_heights

    left = gs[2, 0].subgridspec(1 + len(all_metadata), 1,
                                height_ratios=left_height_ratios,
                                hspace=0.35)

    ax_abund_legend = fig.add_subplot(left[0, 0])
    ax_meta_legends = [fig.add_subplot(left[i+1, 0]) for i in range(n_meta_legends)] if n_meta_legends > 0 else []

    # ---- heatmap ----
    im = ax_heat.imshow(heat.values, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)

    # Add grid
    ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heat.grid(which="minor", color="#d0d0d0", linewidth=0.5)
    ax_heat.grid(which="major", visible=False)

    # Limits
    ax_heat.set_xlim(-0.5, safe_cols - 0.5)
    ax_heat.set_ylim(safe_rows - 0.5, -0.5)

    # Labels
    ax_heat.set_xticks(np.arange(n_cols))
    ax_heat.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=9)
    ax_heat.set_yticks(np.arange(n_rows))
    ax_heat.set_yticklabels(heat.index, rotation=0, fontsize=9)
    ax_heat.set_xlabel(rank.capitalize(), labelpad=10)
    ax_heat.set_ylabel("", labelpad=10)
    ax_heat.tick_params(axis='x', which='major', pad=5)

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


    # ---- left: meta legend ----
    for idx, (ax_meta, meta_info) in enumerate(zip(ax_meta_legends, all_metadata if all_metadata else [])):
            ax_meta.set_title(meta_info['title'], fontsize=10, pad=6)
            ax_meta.title.set_x(0.15)

            row_height = 1.4

            cats = list(meta_info['color_map'].keys())
            ax_meta.set_xlim(0, 1)
            ax_meta.set_ylim(0, len(cats) * row_height)

            for i, cat in enumerate(reversed(cats)):
                y = i * row_height
                color = meta_info['color_map'][cat]
                ax_meta.add_patch(
                    Rectangle((0.05, y), 0.25, row_height, color=color, ec="#888888", linewidth=0.5)
                )
                ax_meta.text(0.35, y + row_height/2, str(cat), va="center", fontsize=9)
            ax_meta.axis("off")
        
    if not all_metadata:
        for ax_meta in ax_meta_legends:
            ax_meta.axis("off")

    # ---- Metadata bars between legends and heatmap ----
    if ax_meta_bars is not None and all_metadata:
        ax_meta_bars.set_xlim(-0.5, n_meta_cols - 0.5)
        ax_meta_bars.set_ylim(n_rows - 0.5, -0.5)

        bar_width = 0.4
        
        for col_idx, meta_info in enumerate(all_metadata):
            category_map = meta_info['category_map']
            color_map = meta_info['color_map']
            
            for row_idx, sample in enumerate(heat.index):
                category = category_map.get(sample, None)
                color = color_map.get(category, "#cccccc")
                
                ax_meta_bars.add_patch(
                    Rectangle((col_idx - bar_width/2, row_idx - 0.45), bar_width, 0.9, 
                             color=color, ec="#888888", linewidth=0.3))
        
        ax_meta_bars.set_xticks(np.arange(n_meta_cols))
        ax_meta_bars.set_xticklabels([m['title'] for m in all_metadata], 
                                      rotation=45, ha="right", fontsize=9)
        ax_meta_bars.set_yticks([])
        ax_meta_bars.tick_params(axis='both', which='both', length=0, labelleft=False)
        ax_meta_bars.xaxis.set_label_position('bottom')
        ax_meta_bars.xaxis.tick_bottom()

        for spine in ax_meta_bars.spines.values():
            spine.set_visible(False)

        ax_meta_bars.yaxis.set_visible(False)
    
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

        if log_top:
            yceil = int(np.ceil(ymax))
            ax_top.set_ylim(0, max(1, yceil))
            ax_top.set_yticks(np.arange(0, yceil + 1, 1))
        else:
            ax_top.set_ylim(0, ymax * 1.15)
            ax_top.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    else:
        ax_top.set_ylim(0, 1)
        ax_top.set_yticks([0, 1])

    ax_top.set_ylabel(y_label)

    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.spines["top"].set_visible(False)
    ax_top.grid(axis='y', color='#888888', linewidth=0.5, zorder=10)
    ax_top.set_axisbelow(False)

    # ---- right bar ----
    ypos = np.arange(n_rows)

    if first_meta is not None:
        bar_colors = [
            first_meta['color_map'].get(first_meta['category_map'].get(sample, None), "#668dcc")
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
    # Dynamic tick spacing based on max value
    if max_val <= 10:
        tick_step = 2
        nice_max = max(10, int(np.ceil(max_val / tick_step) * tick_step))
    elif max_val <= 20:
        tick_step = 5
        nice_max = int(np.ceil(max_val / tick_step) * tick_step)
    elif max_val <= 50:
        tick_step = 10
        nice_max = int(np.ceil(max_val / tick_step) * tick_step)
    elif max_val <= 100:
        tick_step = 20
        nice_max = int(np.ceil(max_val / tick_step) * tick_step)
    else:
        tick_step = 50
        nice_max = int(np.ceil(max_val / tick_step) * tick_step)
    
    half = nice_max // 2

    ax_right.set_xlim(0, nice_max)
    ax_right.set_ylim(n_rows - 0.5, -0.5)

    # Dynamic step ticks
    ticks = np.arange(0, nice_max + 1, tick_step)
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

    plt.suptitle(f"MAG distribution: samples × {rank}", y=0.95, fontsize=12)
    out_png = os.path.join(output_path, f"heatmap_with_bars_{rank}_{meta_cols}.{fmt}")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    return out_png
