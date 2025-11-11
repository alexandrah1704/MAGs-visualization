import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import hashlib

def normalize_id_series(s: pd.Series) -> pd.Series:
    """Normalize IDs"""
    s = s.astype(str)
    s = s.str.replace(r".*[\\/]", "", regex=True)
    s = s.str.replace(r"_Count$", "", regex=True)
    s = s.str.replace(r"_fasta$", "", regex=True)
    s = s.str.replace(r"(?i)\.(?:fa|fna|fasta|faa|gz|gbff)$", "", regex=True)
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

def _coerce_bakta_to_wide(bakta_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Brings Bakta into wide form:
    - rows = genomes
    - columns = features
    """
    df = bakta_raw.copy()

    if "Annotation" in df.columns:
        wide = df.set_index("Annotation").T
    elif (df.index.name is not None) and (str(df.index.name).strip().lower() == "annotation"):
        wide = df.T
    elif df.shape[1] >= 2 and any(
        str(x).strip().lower() in {"cdss", "trnas", "rrnas", "hypotheticals", "crispr arrays", "gaps"}
        for x in df.iloc[:, 0].head(10)
    ):
        df = df.rename(columns={df.columns[0]: "Annotation"})
        wide = df.set_index("Annotation").T
    elif (df.columns.to_series().astype(str)
              .str.contains(r"_count$", case=False, regex=True).mean() > 0.5):
        wide = df.T 
    elif (df.index.to_series().astype(str)
              .str.contains(r"_count$", case=False, regex=True).mean() > 0.5):
        wide = df
    else:
        print("[WARN] Couldn't detect Bakta orientation.")
        wide = df

    norm = normalize_id_series(wide.index.to_series())
    wide.index = pd.Index(norm.values, name="Genome")

    # clean feature names
    wide.columns = wide.columns.str.strip()
    for c in wide.columns:
        wide[c] = pd.to_numeric(wide[c], errors="coerce")

    return wide

def stable_jitter(keys: pd.Series, scale=0.06) -> np.ndarray:
    """
    Deterministic jitter in [-scale, scale] per key (string).
    Same key -> same jitter across runs and independent of plot order.
    """
    out = np.empty(len(keys), dtype=float)
    for i, k in enumerate(keys.astype(str)):
        h = hashlib.blake2b(k.encode("utf-8"), digest_size=4).digest()
        u = int.from_bytes(h, "little") / 2**32
        out[i] = (2*u - 1) * scale
    return out


def bakta_annotation_dashboard(dfs: dict, output_dir: str, metrics=None):
    """
    Fuctional annotation Bakta vs completeness (CheckM2)
    x = completeness in bins, y = Bakta feature
    """

    os.makedirs(output_dir, exist_ok=True)

    # --- Normalize Bakta to wide ----
    bakta_raw = dfs["bakta"].copy()
    bakta = _coerce_bakta_to_wide(bakta_raw)
    bakta.columns = bakta.columns.str.strip().str.lower()

    # --- CheckM2 ---
    checkm2 = dfs["checkm2"].copy()
    compl_col = _first_present_col(checkm2, ["Completeness", "completeness", "CheckM2 completeness"])

    # build normalized join keys
    bakta = bakta.reset_index().rename(columns={"Genome": "__id__"})
    id_col = _first_present_col(checkm2, ["Genome", "user_genome", "Bin Id", "Name", "genome"])
    checkm2_ids = normalize_id_series(checkm2[id_col]) if id_col else normalize_id_series(checkm2.index.to_series())
    checkm2 = checkm2.assign(__id__=checkm2_ids)

    # --- Choose features ---
    # Keys = display names, values = column names in bakta
    if metrics is None:
        metrics = {
            "CDS": "cdss",
            "Hypotheticals": "hypotheticals",
            "rRNAs": "rrnas",
            "tRNAs": "trnas",
        }
        metrics = {disp: col for disp, col in metrics.items() if col in bakta.columns}

    # ---- Merge Bakta features with completeness ----
    metric_cols = list(metrics.values())
    bakta_sub = bakta[["__id__"] + metric_cols]
    checkm2_sub = checkm2[["__id__", compl_col]]

    merged = bakta_sub.merge(checkm2_sub, on="__id__", how="inner")

    # ---- Bin completeness ----
    bins = [-np.inf, 50, 70, 90, np.inf]
    labels = ["<50%", "50-70%", "70-90%", "≥90%"]
    merged["Completeness"] = pd.cut(
        pd.to_numeric(merged[compl_col], errors="coerce"),
        bins=bins, labels=labels,
        right=True, include_lowest=True
    )
    merged["Completeness"] = pd.Categorical(merged["Completeness"], categories=labels, ordered=True)
    order = labels

    # ---- long format for plotting ----
    display_map = {v: k for k, v in metrics.items()}
    long = merged.melt(
        id_vars=["__id__", "Completeness"],
        value_vars=metric_cols,
        var_name="Metric", value_name="Value"
    ).dropna(subset=["Value", "Completeness"])
    long["Metric"] = long["Metric"].map(lambda col: display_map.get(col, col))

    metrics_display = [dsp for dsp in metrics.keys() if dsp in long["Metric"].unique()]

    # ---- Plotting ----
    n_metrics = len(metrics_display)
    fig, axes = plt.subplots(1, n_metrics, figsize=(5*n_metrics, 5), sharey=False)
    if n_metrics == 1:
        axes = [axes]

    palette = sns.color_palette("Set2", n_colors=len(order))
    color_map = {lab: palette[i] for i, lab in enumerate(order)}

    for ax, metric in zip(axes, metrics_display):
        data = long[long["Metric"] == metric].copy()
        if data.empty:
            ax.set_axis_off()
            continue

        sns.boxplot(
            data=data, x="Completeness", y="Value",
            hue="Completeness", dodge=False, order=order,
            palette=palette, ax=ax, showfliers=False, legend=False,
            zorder=1
        )

        for patch in ax.artists:
            fc = patch.get_facecolor()
            patch.set_facecolor((fc[0], fc[1], fc[2], 0.6))
            patch.set_zorder(1)
        for line in ax.lines:
            line.set_zorder(1)

        # ---- Scatter Layer ----
        cats = pd.Categorical(data["Completeness"], categories=order, ordered=True)
        xpos = cats.codes.astype(float)
        keys  = data["__id__"].astype(str) + "|" + data["Completeness"].astype(str)
        x_jit = xpos + stable_jitter(keys, scale=0.06)
        point_colors = [color_map[str(cat)] for cat in cats.astype(str)]

        ax.scatter(
            x_jit, data["Value"].to_numpy(),
            s=22, facecolors=point_colors, edgecolors="black",
            linewidths=0.5, alpha=0.9, zorder=3
        )

        # ---- Styling ----
        ax.set_xlabel("Completeness", fontsize=11)
        ax.set_ylabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', rotation=0)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)

        # Log-Scale
        if metric in ["CDS", "Hypotheticals"]:
            ax.set_yscale("log")

    plt.tight_layout()
    out_path = os.path.join(output_dir, "bakta_annotation_dashboard.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')

    debug_csv = os.path.join(output_dir, "bakta_annotation_long.csv")
    long.to_csv(debug_csv, index=False)
