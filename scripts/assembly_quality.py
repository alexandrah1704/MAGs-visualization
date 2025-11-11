# Box + Scatter dashboard for assembly quality
# x = completeness bins, y = metrix
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


CANDIDATE_METRICS = [
    "N50",
    "Total length (>= 0 bp)",
    "# contigs",
    "# contigs (>= 0 bp)",
    "# contigs (>= 1000 bp)",
    "Largest contig",
    "GC (%)",
]

def _normalize_names(seq):
    return [str(s).strip().lower() for s in seq]

def _reshape_quast(quast_df: pd.DataFrame) -> pd.DataFrame:
    """
    Return Quast in wide shape: rows = assemblies, columns = metrics
    """
    df = quast_df.copy()

    if "Assembly" in df.columns:
        df = df.set_index("Assembly").T
    else:
        cols_norm = set(_normalize_names(df.columns))
        idx_norm  = set(_normalize_names(df.index))
        cand_norm = set(_normalize_names(CANDIDATE_METRICS))
        if len(idx_norm & cand_norm) > len(cols_norm & cand_norm):
            df = df.T

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def _normalize_id_series(s: pd.Series, for_quast: bool = False) -> pd.Series:
    s = s.astype(str)
    s = s.str.replace(r".*[\\/]", "", regex=True)
    if for_quast:
        s = s.str.replace(r"(?i)\.(?:fa|fna|fasta|faa).*$", "", regex=True)
    else:
        s = s.str.replace(r"(?i)\.(?:fa|fna|fasta|faa|gz)$", "", regex=True)
    s = s.str.replace(r"\s+", "_", regex=True)
    s = s.str.replace(r"[^\w\-]+", "_", regex=True)
    s = s.str.replace(r"_+", "_", regex=True).str.strip("_").str.lower()
    return s

def _first_present_col(df: pd.DataFrame, candidates) -> str | None:
    lowmap = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    return None

def _ensure_id_col(df: pd.DataFrame, prefer_cols, for_quast: bool = False) -> pd.Series:
    """ Build normalized IDs """
    col = _first_present_col(df, prefer_cols)
    base = df[col] if col is not None else df.index.to_series()
    return _normalize_id_series(base, for_quast=for_quast)


def assembly_quality_dashboard(dfs: dict, output_dir: str, metrics=None):

    os.makedirs(output_dir, exist_ok=True)

    # ---- Load & normalize ----
    quast = _reshape_quast(dfs["quast"]).copy()
    checkm2 = dfs["checkm2"].copy()

    compl_col = _first_present_col(checkm2, ["Completeness", "completeness", "CheckM2 completeness"])

    quast_ids = _normalize_id_series(quast.index.to_series(), for_quast=True)
    quast = quast.assign(__id__=quast_ids)
    checkm2_ids = _ensure_id_col(checkm2, ["Genome", "user_genome", "Bin Id", "Name", "genome"], for_quast=False)
    checkm2 = checkm2.assign(__id__=checkm2_ids)

    # hue_series = None
    # if use_domain_hue and "gtdb" in dfs:
    #     gtdb = dfs["gtdb"].copy()
    #     gtdb_ids = _ensure_id_col(gtdb, ["Genome", "user_genome", "Bin Id", "Name", "genome"])
    #     gtdb = gtdb.assign(__id__=gtdb_ids)
    #     for cand in ["domain", "Domain", "superkingdom"]:
    #         if cand in gtdb.columns:
    #             hue_series = gtdb.set_index("__id__")[cand]
    #             break

    # ---- Select metrics ----
    if metrics is None:
        priority = ["N50", "# contigs", "Largest contig", "GC (%)"]
        metrics = [m for m in priority if m in quast.columns]

    # ---- Merge on normalized IDs ----
    merged = quast[["__id__"] + metrics].merge(
        checkm2[["__id__", compl_col]], on="__id__", how="inner"
    )

    # ---- Bin completeness ----
    bins = [-np.inf, 50, 70, 90, np.inf]
    labels = ["<50%", "50-70%", "70-90%", "≥90%"]
    merged["Completeness"] = pd.Categorical(
        pd.cut(pd.to_numeric(merged[compl_col], errors="coerce"), bins=bins, labels=labels),
        categories=labels, ordered=True
    )
    order = labels

    # ---- long format for plotting ----
    long = merged.melt(
        id_vars=["__id__", "Completeness"],
        value_vars=metrics,
        var_name="Metric", value_name="Value"
    ).dropna(subset=["Value", "Completeness"])
    
    # ---- Plotting ----
    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(5*n_metrics, 5), sharey=False)
    if n_metrics == 1:
        axes = [axes]
    
    palette = sns.color_palette("Set2", n_colors=len(order))
    color_map = {lab: palette[i] for i, lab in enumerate(order)}
    rng = np.random.default_rng(42)

    for ax, metric in zip(axes, metrics):
        data = long[long["Metric"] == metric].copy()
        if data.empty:
            ax.set_axis_off()
            continue

        # ---- Box Layer, hue = completeness ----
        bp = sns.boxplot(
            data=data, x="Completeness", y="Value",
            hue="Completeness", dodge=False, order=order,
            palette=palette, ax=ax, showfliers=False, legend=False,
            zorder=1
        )
        for patch in ax.artists:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, 0.6))
            patch.set_zorder(1)
        for line in ax.lines:
            line.set_zorder(1)

        # ---- Scatter Layer ----
        cats = pd.Categorical(data["Completeness"], categories=order, ordered=True)
        xpos = cats.codes.astype(float)
        x_jit = xpos + rng.normal(0, 0.06, size=len(xpos))
        point_colors = [color_map[str(cat)] for cat in cats.astype(str)]

        ax.scatter(
            x_jit, data["Value"].values,
            s=22, facecolors=point_colors,
            edgecolors="black", linewidths=0.5, alpha=0.9,
            zorder=3
        )

        # ---- Styling ----
        if any(k in metric.lower() for k in ["length", "n50", "largest", "contig"]) and "gc" not in metric.lower():
            ax.set_yscale("log")

        ax.set_xlabel("Completeness", fontsize=11)
        ax.set_ylabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', rotation=0)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)

    
    plt.tight_layout()
    out_path = os.path.join(output_dir, "assembly_quality_dashboard.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
