# Bakta functional annotation vs completeness (checkm2)
# x = completeness bins, y = bakta features
# Coloring options: quality, tax, metadata
import numpy as np
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import hashlib
from .id_normalizer import normalize_genome_id, normalize_genome_series
from pandas.api.types import is_bool_dtype


def normalize_id_series(s: pd.Series) -> pd.Series:
    return normalize_genome_series(s)

def normalize_id_universal(s: str) -> str:
    return normalize_genome_id(s)

def _first_present_col(df: pd.DataFrame, candidates) -> str | None:
    lowmap = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    return None

def _coerce_bakta_to_wide(bakta_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Brings Bakta into wide form:
    - rows = genomes
    - columns = bakta features
    """
    df = bakta_raw.copy()
    wide = None

    feature_names = {
        "cdss", "trnas", "rrnas", "hypotheticals",
        "crispr arrays", "gaps", "ncrna regions", "ncrnas",
        "orics", "orits", "orivs"
    }
    if any(str(x).strip().lower() in feature_names for x in df.index[:10]):
        drop_mask = df.index.to_series().astype(str).str.strip().str.lower().isin(
            {"#key", "key", "annotation"}
        )
        df2 = df[~drop_mask]
        wide = df2.T

    elif "Annotation" in df.columns:
        wide = df.set_index("Annotation").T

    # Case: 'Annotation'
    elif (df.index.name is not None) and (str(df.index.name).strip().lower() == "annotation"):
        wide = df.T

    # Case: first column = typical feature name
    elif df.shape[1] >= 2 and any(
        str(x).strip().lower() in feature_names
        for x in df.iloc[:, 0].head(10)
    ):
        df = df.rename(columns={df.columns[0]: "Annotation"})
        wide = df.set_index("Annotation").T

    # *_Count
    elif (df.columns.to_series().astype(str)
              .str.contains(r"_count$", case=False, regex=True).mean() > 0.5):
        wide = df.T
    elif (df.index.to_series().astype(str)
              .str.contains(r"_count$", case=False, regex=True).mean() > 0.5):
        wide = df

    # Fallback
    if wide is None:
        print("[WARN] Couldn't detect Bakta orientation. Using raw table.")
        wide = df

    # Genome-IDs normalizing
    norm = normalize_id_series(wide.index.to_series())
    wide.index = pd.Index(norm.values, name="Genome")

    wide.columns = wide.columns.str.strip()
    for c in wide.columns:
        wide[c] = pd.to_numeric(wide[c], errors="coerce")

    return wide


def stable_jitter(keys: pd.Series, scale=0.06) -> np.ndarray:
    """ Deterministic jitter """
    out = np.empty(len(keys), dtype=float)
    for i, k in enumerate(keys.astype(str)):
        h = hashlib.blake2b(k.encode("utf-8"), digest_size=4).digest()
        u = int.from_bytes(h, "little") / 2**32
        out[i] = (2*u - 1) * scale
    return out


def bakta_annotation_plot(dfs: dict, output_dir: str, metrics=None,
                               color_by: str | None = None,
                               tax_rank: str = "phylum",
                               meta_col: str | None = None,
                               fig_size = (5,5), fmt: str ="png",
                               plot_ratio=False,):
    """
    Plot Bakta features vs completeness (checkm2)
    """

    os.makedirs(output_dir, exist_ok=True)

    # --- Normalize Bakta to wide ----
    bakta_raw = dfs["bakta"].copy()
    bakta = _coerce_bakta_to_wide(bakta_raw)
    bakta.columns = bakta.columns.str.strip().str.lower()

    # --- CheckM2 ---
    checkm2 = dfs["checkm2"].copy()
    compl_col = _first_present_col(checkm2, ["Completeness", "completeness", "CheckM2 completeness"])
    cont_col = _first_present_col(checkm2, ["Contamination", "contamination", "CheckM2 contamination"])

    # passende ID-Spalte für CheckM2 suchen (analog wie im dRep-Code)
    if checkm2.index.name and str(checkm2.index.name).strip().lower() not in ("", "index"):
        id_series = checkm2.index.to_series()
    else:
        id_col = _first_present_col(
            checkm2,
            ["Name", "Bin Id", "Bin_Id", "Bin", "Genome", "user_genome"]
        )
        if id_col is None:
            raise ValueError(
                "Could not find genome ID column in CheckM2 table "
                "(tried Name / Bin Id / Bin / Genome / user_genome)."
            )
        id_series = checkm2[id_col]

    checkm2["__id__"] = id_series.astype(str).apply(normalize_id_universal)

    # build normalized join keys
    bakta = bakta.reset_index().rename(columns={"Genome": "__id__"})
    bakta["__id__"] = bakta["__id__"].apply(normalize_id_universal)

    # --- Choose features ---
    if metrics is None:
        metrics = {
            "CDS": "cdss",
            "Hypotheticals": "hypotheticals",
            "rRNAs": "rrnas",
            "tRNAs": "trnas",
        }
        metrics = {disp: col for disp, col in metrics.items() if col in bakta.columns}

    metric_cols = list(metrics.values())

    metric_cols_for_merge = metric_cols.copy()
    if "cdss" in bakta.columns and "cdss" not in metric_cols_for_merge:
        metric_cols_for_merge.append("cdss")

    bakta_sub = bakta[["__id__"] + metric_cols_for_merge]

    # ---- Merge Bakta features with completeness ----
    cols_keep = ["__id__", compl_col]
    if cont_col is not None:
        cols_keep.append(cont_col)
    checkm2_sub = checkm2[cols_keep]

    merged = bakta_sub.merge(checkm2_sub, on="__id__", how="inner")
    
    # ---- calculate Ratios: Metric / cdss ----
    ratio_cols = {}

    if "cdss" in merged.columns:
        for disp_name, col_name in metrics.items():
            # CDS not dividing by themselves
            if col_name == "cdss":
                continue
            if col_name in merged.columns:
                ratio_col_name = f"{col_name}_per_cds"
                merged[ratio_col_name] = (
                    merged[col_name] / merged["cdss"].replace(0, np.nan)
                )
                ratio_cols[disp_name] = ratio_col_name
    else:
        print("[WARN] No 'cdss' column in merged Bakta+CheckM2 table -> no ratios will be plotted.")
        print("Available columns in merged:", list(merged.columns))

    # ---- Bin completeness ----
    bins = [-np.inf, 50, 70, 90, np.inf]
    labels = ["<50%", "50-70%", "70-90%", "≥90%"]
    merged["Completeness_bin"] = pd.Categorical(
        pd.cut(
            pd.to_numeric(merged[compl_col], errors="coerce"),
            bins=bins,
            labels=labels,
        ),
        categories=labels,
        ordered=True,
    )

    # ---- Color logic ----
    color_by_effective = color_by
    quality_col = None
    tax_col = None
    tax_count_map: dict[str, int] | None = None

    # a) color_by: quality - hq, mq, lq
    if color_by_effective == "quality" and cont_col is not None:
        c = pd.to_numeric(merged[compl_col], errors="coerce")
        k = pd.to_numeric(merged[cont_col], errors="coerce")

        q = np.full(len(merged), "Low-quality", dtype=object)
        q[(c >= 70)] = "Medium-quality"
        q[(c >= 90) & (k <= 5)] = "High-quality"
        q[pd.isna(c) | pd.isna(k)] = "Unknown"

        merged["Quality"] = pd.Categorical(
            q,
            categories=["Low-quality", "Medium-quality", "High-quality", "Unknown"],
            ordered=True,
        )
        quality_col = "Quality"

    # b) color_by_ tax
    elif color_by_effective == "tax":
        if "gtdb" not in dfs:
            print("[WARN] No GTDB table found → back to completeness.")
            color_by_effective = None
        else:
            from .heatmap import extract_rank
            gtdb = dfs["gtdb"].copy()
            gtdb = gtdb.rename_axis("orig_id").reset_index()
            gtdb["__id__"] = gtdb["orig_id"].astype(str).apply(normalize_id_universal)

            if "classification" in gtdb.columns:
                gtdb[tax_rank] = gtdb["classification"].apply(
                    lambda tx: extract_rank(tx, tax_rank)
                )
                tax_col = tax_rank
            else:
                tax_col = None
                for c in [tax_rank, tax_rank.capitalize(), tax_rank.lower()]:
                    if c in gtdb.columns:
                        tax_col = c
                        break

            if tax_col is None:
                print(f"[WARN] Could not find taxonomy for rank '{tax_rank}'.")
                color_by_effective = None
            else:
                merged = merged.merge(
                    gtdb[["__id__", tax_col]],
                    on="__id__",
                    how="left"
                )
                merged[tax_col] = merged[tax_col].fillna("Unknown")

                top_n = 10
                value_counts = merged[tax_col].value_counts()
                top_taxa = value_counts.head(top_n).index.tolist()
                merged[tax_col] = merged[tax_col].apply(
                    lambda x: x if x in top_taxa else "Others"
                )
                tax_count_map = merged[tax_col].value_counts().to_dict()

    # c) color_by: metadata
    elif color_by_effective == "meta":
        if "metadata" not in dfs:
            print("[WARN] No metadata table found → falling back to completeness.")
            color_by_effective = None
        else:
            metadata = dfs["metadata"].copy()
            
            if "user_genome" in metadata.columns and metadata.index.name != "user_genome":
                metadata.set_index("user_genome", inplace=True)
            
            meta_index_str = metadata.index.astype(str)
            has_bin = meta_index_str.str.contains("_bin_").any()

            def sample_key(s: str):
                if pd.isna(s):
                    return s
                s_norm = normalize_id_universal(str(s))
                return re.split(r'[_\.]', s_norm, 1)[0]

            if has_bin:
                metadata["__key__"] = metadata.index.to_series().map(normalize_id_universal)
                merged["__key__"]   = merged["__id__"].astype(str).map(normalize_id_universal)
            else:
                metadata["__key__"] = metadata.index.to_series().map(sample_key)
                merged["__key__"]   = merged["__id__"].astype(str).map(sample_key)

            if meta_col not in metadata.columns:
                print(f"[WARN] Column '{meta_col}' not found in metadata → available: {metadata.columns.tolist()}")
                color_by_effective = None
            else:
                merged = merged.merge(
                    metadata[["__key__", meta_col]],
                    on="__key__",
                    how="left"
                )
                
                raw = merged[meta_col]
                raw_non_na = raw.dropna()
                unknown_label = f"Unknown {meta_col}"

                bool_like = False
                if is_bool_dtype(raw_non_na):
                    bool_like = True
                else:
                    lowered = {str(v).strip().lower() for v in pd.unique(raw_non_na)}
                    bool_tokens = {"true", "false", "yes", "no", "0", "1"}
                    if lowered and lowered.issubset(bool_tokens):
                        bool_like = True

                if bool_like:
                    merged[meta_col] = (
                        raw.astype("string")
                           .fillna(unknown_label)
                           .replace("", unknown_label)
                    )
                else:
                    meta_num = pd.to_numeric(raw, errors="coerce")
                    frac_numeric = meta_num.notna().mean()

                    if frac_numeric >= 0.8 and meta_num.notna().sum() > 0:
                        bin_width = 5.0
                        vmin = np.floor(meta_num.min() / bin_width) * bin_width
                        vmax = np.ceil(meta_num.max() / bin_width) * bin_width
                        if vmin == vmax:
                            vmax = vmin + bin_width
                        
                        bins_meta = np.arange(vmin, vmax + bin_width, bin_width)
                        labels_meta = [
                            f"{int(left)}–{int(right)}"
                            for left, right in zip(bins_meta[:-1], bins_meta[1:])
                        ]
                        
                        merged[meta_col] = pd.cut(
                            meta_num,
                            bins=bins_meta,
                            labels=labels_meta,
                            include_lowest=True
                        ).astype("object").fillna(unknown_label)  # type: ignore
                    else:
                        merged[meta_col] = (
                            raw.astype("string")
                               .fillna(unknown_label)
                               .replace("", unknown_label)
                        )

                top_n_meta = 10
                counts = merged[meta_col].value_counts()
                if len(counts) > top_n_meta:
                    top_cats = counts.nlargest(top_n_meta).index.tolist()
                    merged[meta_col] = merged[meta_col].apply(
                        lambda x: x if x in top_cats else "Others"
                    )

    # ---- long format for plotting ----
    display_map = {v: k for k, v in metrics.items()}
    
    id_vars = ["__id__", "Completeness_bin"]
    if quality_col:
        id_vars.append(quality_col)
    if tax_col:
        id_vars.append(tax_col)
    if meta_col and meta_col in merged.columns:
        id_vars.append(meta_col)
    
    long = merged.melt(
        id_vars=id_vars,
        value_vars=metric_cols,
        var_name="Metric", value_name="Value"
    ).dropna(subset=["Value", "Completeness_bin"])
    
    long["Metric"] = long["Metric"].map(lambda col: display_map.get(col, col))
    metrics_display = [dsp for dsp in metrics.keys() if dsp in long["Metric"].unique()]

    # ---- Plot each metric separately ----
    for metric in metrics_display:
        data = long[long["Metric"] == metric].copy()
        if data.empty:
            print(f"[WARN] No data for metric {metric}")
            continue

        # color setup - depending on color_by
        if color_by_effective == "quality" and quality_col:
            label_series = data[quality_col].astype(str)
            categories_color = [
                lab
                for lab in ["Low-quality", "Medium-quality", "High-quality", "Unknown"]
                if lab in label_series.unique()
            ]
            pal = sns.color_palette("Set2", len(categories_color))
            legend_title = "Quality"

        elif color_by_effective == "tax" and tax_col:
            label_series = data[tax_col].astype(str)
            categories_color = sorted(label_series.unique())
            pal = sns.color_palette("tab20", len(categories_color))
            legend_title = tax_rank

        elif color_by_effective == "meta" and meta_col:
            label_series = data[meta_col].astype(str)
            categories_color = sorted(label_series.unique())
            pal = sns.color_palette("tab20", len(categories_color))
            legend_title = meta_col

        else:
            label_series = data["Completeness_bin"].astype(str)
            categories_color = labels
            pal = sns.color_palette("Set2", len(categories_color))
            legend_title = None

        color_map = {lab: pal[i] for i, lab in enumerate(categories_color)}
        point_colors = [color_map.get(v, "#555555") for v in label_series]

        if plot_ratio:
            fig, (ax_top, ax_bottom) = plt.subplots(
                2, 1,
                figsize=(fig_size[0], fig_size[1] * 1.4),
                sharex=True,
                gridspec_kw={'height_ratios': [3, 1]}
            )
        else:
            fig, ax_top = plt.subplots(figsize=fig_size)
            ax_bottom = None

        #  Top Boxplot
        sns.boxplot(
            data=data,
            x="Completeness_bin",
            y="Value",
            order=labels,
            color="#dddddd",
            ax=ax_top,
            showfliers=False,
            dodge=False,
            zorder=1,
        )
        
        # Jitter positions
        cats_x = pd.Categorical(data["Completeness_bin"], categories=labels, ordered=True)
        xpos = cats_x.codes.astype(float)
        keys = data["__id__"].astype(str) + "|" + data["Completeness_bin"].astype(str)
        x_jit = xpos + stable_jitter(keys, scale=0.06)

        # Scatter points on top of boxplot
        ax_top.scatter(
            x_jit,
            data["Value"],
            s=22,
            facecolors=point_colors,
            edgecolors="black",
            linewidths=0.5,
            alpha=0.9,
            zorder=3,
        )

        # Log scale
        if metric in ["CDS", "Hypotheticals"]:
            ax_top.set_yscale("log")

        ax_top.set_ylabel(metric)
        ax_top.set_title(metric, fontweight="bold")
        ax_top.grid(axis="y", linestyle="--", alpha=0.3)

        # Legend
        if legend_title:
            handles = []
            legend_labels = categories_color

            if color_by_effective == "tax" and tax_col and tax_count_map is not None:
                legend_labels = [
                    f"{lab} ({tax_count_map.get(lab, 0)})"
                    for lab in categories_color
                ]

            for lab, disp_lab in zip(categories_color, legend_labels):
                pt = plt.Line2D(
                    [], [], marker="o", linestyle="", color=color_map[lab],
                    label=disp_lab, markersize=6
                )
                handles.append(pt)

            ax_top.legend(
                handles,
                legend_labels,
                title=legend_title,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.15),
                frameon=False,
                ncol=min(3, len(legend_labels)),
                fontsize=7,
                title_fontsize=8,
                markerscale=0.7,
            )
        
        # Bottom plot
        if plot_ratio and ax_bottom is not None:
            ratio_col = ratio_cols.get(metric, None)

            if plot_ratio and ratio_col and ratio_col in merged.columns:
                ratio_stats = (
                    merged
                    .groupby("Completeness_bin")[ratio_col]
                    .median()
                    .reindex(labels)
                )

                ax_bottom.plot(
                    range(len(labels)),
                    ratio_stats.values,
                    marker="o",
                    linestyle="-",
                    linewidth=1.5,
                )

                ax_bottom.set_ylabel(f"{metric} / CDS")
                rmax = np.nanmax(ratio_stats.values)
                if np.isfinite(rmax):
                    ax_bottom.set_ylim(0, rmax * 1.2)
                ax_bottom.grid(axis="y", linestyle="--", alpha=0.3)
            else:
                ax_bottom.axis("off")

        if plot_ratio and ax_bottom is not None:
            ax_x = ax_bottom
        else:
            ax_x = ax_top
            
        ax_x.set_xlabel("Completeness")
        ax_x.set_xticks(range(len(labels)))
        ax_x.set_xticklabels(labels)


        plt.tight_layout()
        safe_metric = metric.replace(" ", "_").replace("/", "_")
        if color_by == "meta":
            out = os.path.join(output_dir, f"bakta_{safe_metric}_{color_by}_{meta_col}.{fmt}")
        elif color_by == "tax":
            out = os.path.join(output_dir, f"bakta_{safe_metric}_{color_by}_{tax_rank}.{fmt}")
        elif color_by == "quality":
            out = os.path.join(output_dir, f"bakta_{safe_metric}_{color_by}.{fmt}")
        else:
            out = os.path.join(output_dir, f"bakta_{safe_metric}.{fmt}")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] Saved: {out}")
