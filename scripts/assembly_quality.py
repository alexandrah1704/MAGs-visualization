# Box + Scatter for assembly quality
# x = completeness bins (Checkm2), y = QUAST metrics
# Coloring options: quality, tax, metadata
import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from heatmap import extract_rank
from id_normalizer import normalize_genome_id
from pandas.api.types import is_bool_dtype


def normalize_id_universal(s: str) -> str:
    return normalize_genome_id(s)

def _normalize_names(seq):
    return [str(s).strip().lower() for s in seq]

def _reshape_quast(quast_df: pd.DataFrame) -> pd.DataFrame:
    """Return QUAST in wide shape: rows = assemblies, columns = metrics."""
    df = quast_df.copy()

    if "Assembly" in df.columns:
        df = df.set_index("Assembly").T
    else:
        cols_norm = set(_normalize_names(df.columns))
        idx_norm = set(_normalize_names(df.index))
        cand_norm = set(
            _normalize_names(
                [
                    "N50",
                    "Total length (>= 0 bp)",
                    "# contigs",
                    "# contigs (>= 0 bp)",
                    "# contigs (>= 1000 bp)",
                    "Largest contig",
                    "GC (%)",
                ]
            )
        )
        if len(idx_norm & cand_norm) > len(cols_norm & cand_norm):
            df = df.T
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _first_present_col(df: pd.DataFrame, candidates) -> str | None:
    lowmap = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    return None


def assembly_quality_plot(dfs: dict, output_dir: str, metrics=None,
                               color_by: str | None = None,
                               tax_rank: str = "phylum",
                               meta_col: str | None = None,
                               fig_size=(5,5), fmt: str = "png",):
    """
    Create box + scatter plots of quast metrics vs checkm2 completeness
    """

    os.makedirs(output_dir, exist_ok=True)

    # ---- Load + reshape quast, checkm2 ----
    quast = _reshape_quast(dfs["quast"]).copy()
    checkm2 = dfs["checkm2"].copy()

    compl_col = _first_present_col(
        checkm2, ["CheckM2 completeness", "Completeness", "completeness"]
    )

    cont_col = _first_present_col(
        checkm2, ["CheckM2 contamination", "Contamination", "contamination"]
    )

    if compl_col is None:
        raise ValueError("No completeness column found in CheckM2")

    # ---- Normalize IDs ----
    quast_ids = quast.index.to_series().apply(normalize_id_universal)
    quast["__id__"] = quast_ids

    checkm2_ids = checkm2.index.to_series().apply(normalize_id_universal)
    checkm2["__id__"] = checkm2_ids

    # ---- Choose metrics to plot (by default) ----
    available_metrics = [c for c in quast.columns if c != "__id__"]

    if metrics is None:
        default_priority = ["N50", "# contigs", "Largest contig", "GC (%)"]
        metrics = [m for m in default_priority if m in available_metrics]
    else:
        metrics = [m for m in metrics if m in available_metrics]
        if not metrics:
            raise ValueError("No valid metrics found in QUAST table.")

    # ---- Merge quast + checkm2 ----
    cols_keep = ["__id__", compl_col]
    if cont_col is not None:
        cols_keep.append(cont_col)

    merged = quast[["__id__"] + metrics].merge(
        checkm2[cols_keep],
        on="__id__",
        how="inner",
    )

    # no overlap
    if merged.empty:
        print("[WARN] No overlap between QUAST and CheckM2 after ID normalization.")
        return

    # ---- Completeness bins ----
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

    # b) color_by: taxonomy
    if color_by_effective == "tax":
        if "gtdb" not in dfs:
            print("[WARN] No GTDB table found → back to completeness.")
            color_by_effective = None
        else:
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
    if color_by_effective == "meta":
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
                """SRR-Level-Schlüssel: normalisieren und Stamm vor erstem '_' oder '.'."""
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

                # ---- Bool / numeric / categorical Erkennung ----
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
                        # numeric (Temperature)
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
                        ).astype("object").fillna(unknown_label)
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

    # ---- convert to long format ----
    id_vars = ["__id__", "Completeness_bin"]
    if quality_col:
        id_vars.append(quality_col)
    if tax_col:
        id_vars.append(tax_col)
    if meta_col and meta_col in merged.columns:
        id_vars.append(meta_col)

    long = merged.melt(
        id_vars=id_vars,
        value_vars=metrics,
        var_name="Metric",
        value_name="Value",
    ).dropna(subset=["Value", "Completeness_bin"])

    # ---- Plot ----
    rng = np.random.default_rng(42)

    for metric in metrics:
        data = long[long["Metric"] == metric].copy()
        if data.empty:
            print(f"[WARN] No data for metric {metric}")
            continue

        fig, ax = plt.subplots(figsize=fig_size)

        # Boxplots
        sns.boxplot(
            data=data,
            x="Completeness_bin",
            y="Value",
            order=labels,
            color="#dddddd",
            ax=ax,
            showfliers=False,
            dodge=False,
            zorder=1,
        )

        # Jitter, so points don't overlap too much
        cats_x = pd.Categorical(data["Completeness_bin"], categories=labels, ordered=True)
        xpos = cats_x.codes.astype(float)
        x_jit = xpos + rng.normal(0, 0.06, len(xpos))

        # ---- color setup - depending on color_by ----
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
            # default = color by completeness
            label_series = data["Completeness_bin"].astype(str)
            categories_color = labels
            pal = sns.color_palette("Set2", len(categories_color))
            legend_title = None

        color_map = {lab: pal[i] for i, lab in enumerate(categories_color)}
        point_colors = [color_map.get(v, "#555555") for v in label_series]

        # Scatter points on top of boxplot
        ax.scatter(
            x_jit,
            data["Value"],
            s=22,
            facecolors=point_colors,
            edgecolors="black",
            linewidths=0.5,
            alpha=0.9,
            zorder=3,
        )

        # log-scale
        if any(k in metric.lower() for k in ["n50", "contig", "length"]) and "gc" not in metric.lower():
            ax.set_yscale("log")

        ax.set_xlabel("Completeness")
        ax.set_ylabel(metric)
        ax.set_title(metric, fontweight="bold")

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)

        ax.grid(axis="y", linestyle="--", alpha=0.3)

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

            ax.legend(
                handles,
                legend_labels,
                title=legend_title,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.15),
                frameon=False,
                ncol=min(3, len(legend_labels)),
                fontsize=8,
                title_fontsize=9,
            )

        plt.tight_layout()
        if color_by == "meta":
            out = os.path.join(output_dir, f"assembly_quality_{metric}_{color_by}_{meta_col}.{fmt}")
        elif color_by == "tax":
            out = os.path.join(output_dir, f"assembly_quality_{metric}_{color_by}_{tax_rank}.{fmt}")
        elif color_by == "quality":
            out = os.path.join(output_dir, f"assembly_quality_{metric}_{color_by}.{fmt}")
        else:
            out = os.path.join(output_dir, f"assembly_quality_{metric}.{fmt}")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] Saved: {out}")
