import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec, cm
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from id_normalizer import normalize_genome_id


def extract_rank(tax, rank: str):
    """Extract a GTDB rank from a taxonomy string."""
    if pd.isna(tax):
        return "Unclassified"
    
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

def get_cluster_majority_quality(cluster_df):
    """
    Determine majority quality class of cluster.
    """
    # Nur MAGs mit bekannter Quality
    known = cluster_df[cluster_df['Quality'].isin(['High', 'Medium', 'Low'])]
    if known.empty:
        return 'Unknown'
    
    counts = known['Quality'].value_counts()
    return counts.idxmax()

# ---- Main plot function ----
def drep_cluster_plot(drep_df: pd.DataFrame, gtdb_df: pd.DataFrame, output_path: str, 
                     tax_level: str = "genus", top_n: int = 30, fmt: str = "png",
                     fig_size=None, checkm2_df: pd.DataFrame = None,
                     quast_df: pd.DataFrame = None, bakta_df: pd.DataFrame = None,):
    """
    Create a horizontal bar plot showing the top N clusters by member count,
    annotated with:
    - representative tax (majority within cluster)
    - majority checkm2 quality per cluster
    - average genome size and GC % from QUAST
    - average CDS, rRNA counts and hypotheticals/CDS ratio
    """
    
    os.makedirs(output_path, exist_ok=True)
    
    # ---- dRep & gtdb ----
    drep = drep_df.copy()
    gtdb = gtdb_df.copy()
    
    print(f"[INFO] Creating dRep cluster plot (top {top_n}, tax_level={tax_level})")
    
    if drep.index.name is not None:
        drep = drep.reset_index()
    
    # Normalize genome IDs
    if 'genome' in drep.columns:
        drep['genome_normalized'] = drep['genome'].apply(normalize_genome_id)
    else:
        raise ValueError("dRep DataFrame must contain 'genome' column")
    
    # Count members per cluster
    cluster_counts = drep.groupby('primary_cluster').size().reset_index(name='n_members')
    cluster_counts = cluster_counts.sort_values('n_members', ascending=False).head(top_n)
    
    print(f"[INFO] Found {len(drep['primary_cluster'].unique())} total clusters, showing top {top_n}")
    
    if gtdb.index.name is not None:
        gtdb = gtdb.reset_index()
    
    if 'user_genome' not in gtdb.columns or 'classification' not in gtdb.columns:
        raise ValueError("GTDB file must contain 'user_genome' and 'classification' columns")
    
    # Normalize genome IDs in gtdb
    gtdb['user_genome_normalized'] = gtdb['user_genome'].apply(normalize_genome_id)
    
    # Extract taxonomic information
    gtdb[tax_level] = gtdb['classification'].apply(lambda x: extract_rank(x, tax_level))
    
    # Merge drep with GTDB
    drep_with_tax = drep.merge(
        gtdb[['user_genome_normalized', tax_level]], 
        left_on='genome_normalized', 
        right_on='user_genome_normalized', 
        how='left'
    )

    # For each cluster, find the best taxonomy (prefer classified over "Unclassified")
    def get_cluster_rep_and_tax(cluster_df):
        """
        Pick representative MAG and majority tax for cluster.
        - Prefer MAGs with classified taxonomy != unclassified
        - Take majority tax among classified MAGs
        - If every MAG is unclassified in cluster: take first MAG
        """
        tmp = cluster_df.copy()
        tmp[tax_level] = tmp[tax_level].fillna('Unclassified')
        
        classified = tmp[tmp[tax_level] != 'Unclassified']
        
        if not classified.empty:
            majority_tax = classified[tax_level].mode().iloc[0]
            rep_row = classified[classified[tax_level] == majority_tax].iloc[0]
        else:
            rep_row = tmp.iloc[0]
            majority_tax = rep_row[tax_level] if pd.notna(rep_row[tax_level]) else 'Unclassified'
        
        return pd.Series({
            'representative': rep_row['genome_normalized'],
            tax_level: majority_tax,
        })

    cluster_rep = drep_with_tax.groupby('primary_cluster').apply(
        get_cluster_rep_and_tax
        ).reset_index()

    # Combine cluster size & representative taxonomy
    cluster_data = cluster_counts.merge(cluster_rep, on='primary_cluster', how='left')
    cluster_data[tax_level] = cluster_data[tax_level].fillna('Unclassified')
    cluster_data = cluster_data.sort_values('n_members', ascending=True)
    
    print(f"[INFO] Taxonomic breakdown:")
    print(cluster_data[tax_level].value_counts())

    # ---- Checkm2 - majority cluster quality ----
    rep_quality_available = False
    if checkm2_df is not None:
        print("[INFO] Processing CheckM2 data for representative quality")
        checkm2 = checkm2_df.copy()
        
        # Reset index if needed
        if checkm2.index.name is not None:
            checkm2 = checkm2.reset_index()
            checkm2.rename(columns={checkm2.columns[0]: 'Name'}, inplace=True)

        if 'Name' in checkm2.columns:
            checkm2['genome_normalized'] = checkm2['Name'].apply(normalize_genome_id)
        else:
            print("[WARN] CheckM2 DataFrame has no 'Name' column, skipping quality")
            checkm2 = None
        
        rep_quality_available = False
        
        if checkm2 is not None:
            drep_with_quality = drep_with_tax.merge(
                checkm2[['genome_normalized', 'Completeness', 'Contamination']],
                on='genome_normalized',
                how='left'
            )

            def classify_quality(row):
                if pd.isna(row['Completeness']) or pd.isna(row['Contamination']):
                    return 'Unknown'
                comp = row['Completeness']
                cont = row['Contamination']
                
                if comp >= 90 and cont <= 5:
                    return 'High'
                elif comp >= 70:
                    return 'Medium'
                else:
                    return 'Low'
            
            drep_with_quality['Quality'] = drep_with_quality.apply(classify_quality, axis=1)

            cluster_quality_majority = drep_with_quality.groupby('primary_cluster').apply(
                get_cluster_majority_quality
            ).reset_index()

            cluster_quality_majority.columns = ['primary_cluster', 'ClusterQuality']
            cluster_data = cluster_data.merge(cluster_quality_majority, on='primary_cluster', how='left')
            cluster_data['ClusterQuality'] = cluster_data['ClusterQuality'].fillna('Unknown')
            rep_quality_available = True
    
    # ---- QUAST: avg_genome_size, avg_gc ----
    quast_available = False
    size_norm = None
    gc_norm = None
    size_cmap = None
    gc_cmap = None

    if quast_df is not None:
        quast = quast_df.copy()

        gc_row_name = None
        size_row_name = None

        for name in quast.index:
            lname = str(name).lower()
            if gc_row_name is None and 'gc' in lname and '%' in lname:
                gc_row_name = name
            if size_row_name is None and 'total length' in lname:
                size_row_name = name

        if gc_row_name is not None and size_row_name is not None:
            gc_series = quast.loc[gc_row_name]
            size_series = quast.loc[size_row_name]

            quast_long = pd.DataFrame({
                'genome_raw': gc_series.index,
                'gc': pd.to_numeric(gc_series.values, errors='coerce'),
                'genome_size': pd.to_numeric(size_series.values, errors='coerce'),
            })

            quast_long['genome_normalized'] = quast_long['genome_raw'].apply(normalize_genome_id)

            drep_with_quast = drep.merge(
                quast_long[['genome_normalized', 'gc', 'genome_size']],
                on='genome_normalized',
                how='left'
            )

            cluster_quast = (
                drep_with_quast
                .groupby('primary_cluster')[['gc', 'genome_size']]
                .mean()
                .reset_index()
                .rename(columns={'gc': 'avg_gc', 'genome_size': 'avg_genome_size'})
            )

            cluster_data = cluster_data.merge(cluster_quast, on='primary_cluster', how='left')
            quast_available = True
        else:
            print("[WARN] Could not find GC (%) and/or Total length rows in QUAST file; skipping avg GC/genome size.")
    
    # ---- Bakta: avg_cds, avg_rrna, hypotheticals/CDS ----
    bakta_available = False
    if bakta_df is not None:
        bakta = bakta_df.copy()

        # help function
        def _find_row(name_substr):
            name_substr = name_substr.lower()
            for idx in bakta.index:
                if name_substr in str(idx).lower():
                    return idx
            return None

        cds_row_name  = _find_row("cds")
        rrna_row_name = _find_row("rrna")
        hypo_row_name = _find_row("hypothetical")

        if cds_row_name is not None and rrna_row_name is not None:
            cds_series  = bakta.loc[cds_row_name]
            rrna_series = bakta.loc[rrna_row_name]
            hypo_series = bakta.loc[hypo_row_name] if hypo_row_name is not None else None

            bakta_long = pd.DataFrame({
                "genome_raw": cds_series.index,
                "cds":  pd.to_numeric(cds_series.values,  errors="coerce"),
                "rrna": pd.to_numeric(rrna_series.values, errors="coerce"),
            })

            # Add hypotheticals if available
            has_hypotheticals = False
            if hypo_series is not None:
                bakta_long["hypotheticals"] = pd.to_numeric(hypo_series.values, errors="coerce")
                has_hypotheticals = True

            # remove "_Count" suffix and normalize
            bakta_long["genome_raw_clean"] = bakta_long["genome_raw"].str.replace("_Count", "", regex=False)
            bakta_long["genome_normalized"] = bakta_long["genome_raw_clean"].apply(normalize_genome_id)

            merge_cols = ["genome_normalized", "cds", "rrna"]
            if has_hypotheticals:
                merge_cols.append("hypotheticals")

            drep_with_bakta = drep.merge(
                bakta_long[merge_cols],
                on="genome_normalized",
                how="left",
            )

            agg_cols = ["cds", "rrna"]
            if has_hypotheticals and "hypotheticals" in drep_with_bakta.columns:
                agg_cols.append("hypotheticals")
            
            cluster_bakta = (
                drep_with_bakta
                .groupby("primary_cluster")[agg_cols]
                .mean()
                .reset_index()
            )

            # Calculate ratio AFTER aggregation
            if has_hypotheticals and "hypotheticals" in cluster_bakta.columns:
                cluster_bakta["hypo_cds_ratio"] = cluster_bakta["hypotheticals"] / cluster_bakta["cds"]
                cluster_bakta["hypo_cds_ratio"] = cluster_bakta["hypo_cds_ratio"].replace([np.inf, -np.inf], np.nan).fillna(0)
            
            # Rename columns
            rename_dict = {"cds": "avg_cds", "rrna": "avg_rrna"}
            if has_hypotheticals and "hypotheticals" in cluster_bakta.columns:
                rename_dict["hypotheticals"] = "avg_hypotheticals"
                if "hypo_cds_ratio" in cluster_bakta.columns:
                    rename_dict["hypo_cds_ratio"] = "avg_hypo_cds_ratio"

            cluster_bakta = cluster_bakta.rename(columns=rename_dict)

            # in cluster_data aufnehmen
            cluster_data = cluster_data.merge(cluster_bakta, on="primary_cluster", how="left")
            bakta_available = True
        else:
            print("[WARN] Could not find CDS and/or rRNA rows in Bakta file; skipping Bakta columns.")

    # ---- Figure & Layout ----
    # Create figure
    if fig_size is None:
        if rep_quality_available is not None:
            fig_size = (12, max(7, top_n * 0.28))
        else:
            fig_size = (12, max(8, top_n * 0.3))

    fig = plt.figure(figsize=fig_size)

    # Create gridspec
    if rep_quality_available:
        # [taxonomy labels | bars | spacer | quality]
        gs = gridspec.GridSpec(1, 4, figure=fig, width_ratios=[0.1, 2, 0.3, 0.4], wspace=0.02)
        ax_tax = fig.add_subplot(gs[0, 0])
        ax_bars = fig.add_subplot(gs[0, 1])
        ax_heatmap = fig.add_subplot(gs[0, 3])
    else:
        gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[0.1, 2], wspace=0.02)
        ax_tax = fig.add_subplot(gs[0, 0])
        ax_bars = fig.add_subplot(gs[0, 1])
        ax_heatmap = None

    # ---- Bar Plot ----
    # Y positions
    y_pos = np.arange(len(cluster_data))
    
    # Plot bars
    bars = ax_bars.barh(
        y_pos,
        cluster_data['n_members'].values,
        height=0.7,
        color='#1f2937',
        edgecolor='#111827',
        linewidth=0.5
    )
    
    # Configure bars axis
    max_members = int(cluster_data['n_members'].max())
    nice_max = int(np.ceil(max_members * 1.05))
    step = 5 if nice_max <= 50 else 10

    ax_bars.set_xlim(0, nice_max)
    ax_bars.set_xticks(np.arange(0, nice_max + 1, step))
    ax_bars.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_bars.set_xlabel('# of MAGs in species-level cluster', fontsize=11)

    # Add grid
    ax_bars.grid(axis='x', color='#e5e7eb', linewidth=0.5, zorder=0)
    ax_bars.set_axisbelow(True)
    
    # Remove y-axis labels and ticks from bars
    ax_bars.set_yticks([])
    ax_bars.spines['left'].set_visible(False)
    ax_bars.spines['top'].set_visible(False)
    ax_bars.spines['right'].set_visible(False)
    
    # Add member count labels on bars
    for i, (idx, row) in enumerate(cluster_data.iterrows()):
        count = int(row['n_members'])
        ax_bars.text(
            count + nice_max * 0.01,
            i,
            str(count),
            va='center',
            ha='left',
            fontsize=9,
            fontweight='bold'
        )
    
    # ---- Configure taxonomy axis ----
    ax_tax.set_xlim(0, 1)
    ax_tax.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_tax.set_yticks(y_pos)
    ax_tax.set_yticklabels(cluster_data[tax_level].values, fontsize=9)
    ax_tax.set_xticks([])
    ax_tax.invert_xaxis()
    
    # Remove spines
    for spine in ax_tax.spines.values():
        spine.set_visible(False)
    
    # ---- Heatmap - quality, quast, bakta ---- #
    # Plot quality heatmap if available
    if rep_quality_available and ax_heatmap is not None:

        qual_colors = {
            'High':   "#b64a4a",
            'Medium': "#7f7f7f",
            'Low':    "#86cbd5",
            'Unknown':"#e5e7eb",
        }

        ax_heatmap.set_ylim(-0.5, len(cluster_data) - 0.5)
        ax_heatmap.set_yticks([])
        ax_heatmap.spines['left'].set_visible(False)
        ax_heatmap.spines['top'].set_visible(False)
        ax_heatmap.spines['right'].set_visible(False)

        # --- define horizontal blocks (quality, size, gc,....) ---
        block_width = 1.0
        block_pos = {}
        x = 0

        block_pos['quality'] = x
        x += 1

        if quast_available:
            block_pos['size'] = x
            block_pos['gc']   = x + 1
            x += 2

        if bakta_available:
            block_pos['cds']  = x
            block_pos['rrna'] = x + 1
            x += 2
            if 'avg_hypo_cds_ratio' in cluster_data.columns:
                block_pos['hypo_ratio'] = x
                x += 1

        n_blocks = x
        ax_heatmap.set_xlim(0, n_blocks)
        ax_heatmap.set_xticks([])

        # ---- color norms and cmaps ----
        size_norm = gc_norm = cds_norm = rrna_norm = hypo_ratio_norm = None
        size_cmap = gc_cmap = cds_cmap = rrna_cmap = hypo_ratio_cmap = None

        if quast_available:
            size_vals = cluster_data['avg_genome_size'].values.astype(float)
            gc_vals   = cluster_data['avg_gc'].values.astype(float)
            size_vals_valid = size_vals[~np.isnan(size_vals)]
            gc_vals_valid   = gc_vals[~np.isnan(gc_vals)]
            if len(size_vals_valid) > 0:
                size_norm = Normalize(vmin=size_vals_valid.min(), vmax=size_vals_valid.max())
                size_cmap = cm.get_cmap('Blues')
            if len(gc_vals_valid) > 0:
                gc_norm = Normalize(vmin=gc_vals_valid.min(), vmax=gc_vals_valid.max())
                gc_cmap = cm.get_cmap('Greens')

        if bakta_available:
            cds_vals  = cluster_data['avg_cds'].values.astype(float)
            rrna_vals = cluster_data['avg_rrna'].values.astype(float)
            cds_vals_valid  = cds_vals[~np.isnan(cds_vals)]
            rrna_vals_valid = rrna_vals[~np.isnan(rrna_vals)]
            if len(cds_vals_valid) > 0:
                cds_norm = Normalize(vmin=cds_vals_valid.min(), vmax=cds_vals_valid.max())
                cds_cmap = cm.get_cmap('Purples')
            if len(rrna_vals_valid) > 0:
                rrna_norm = Normalize(vmin=rrna_vals_valid.min(), vmax=rrna_vals_valid.max())
                rrna_cmap = cm.get_cmap('Oranges')
            # Hypotheticals/CDS ratio
            if 'avg_hypo_cds_ratio' in cluster_data.columns:
                hypo_ratio_vals = cluster_data['avg_hypo_cds_ratio'].values.astype(float)
                hypo_ratio_vals_valid = hypo_ratio_vals[~np.isnan(hypo_ratio_vals)]
                if len(hypo_ratio_vals_valid) > 0:
                    hypo_ratio_norm = Normalize(vmin=hypo_ratio_vals_valid.min(), vmax=hypo_ratio_vals_valid.max())
                    hypo_ratio_cmap = cm.get_cmap('Reds')

        # ---- draw rows ----
        for i, row in cluster_data.iterrows():
            # quality
            q = row['ClusterQuality']
            ax_heatmap.barh(
                i,
                block_width,
                left=0,
                height=0.7,
                color=qual_colors.get(q, "#e5e7eb"),
                edgecolor='white',
                linewidth=0.5
            )
            
            # genome size
            if quast_available and size_norm is not None:
                size = row['avg_genome_size']
                col_size = size_cmap(size_norm(size)) if pd.notna(size) else 'white'
                ax_heatmap.barh(
                    i,
                    block_width,
                    left=block_pos['size'],
                    height=0.7,
                    color=col_size,
                    edgecolor='white',
                    linewidth=0.2
                )

            # GC
            if quast_available and gc_norm is not None:
                gc_val = row['avg_gc']
                col_gc = gc_cmap(gc_norm(gc_val)) if pd.notna(gc_val) else 'white'
                ax_heatmap.barh(
                    i,
                    block_width,
                    left=block_pos['gc'],
                    height=0.7,
                    color=col_gc,
                    edgecolor='white',
                    linewidth=0.2
                )

            # CDS
            if bakta_available and cds_norm is not None:
                cds = row['avg_cds']
                col_cds = cds_cmap(cds_norm(cds)) if pd.notna(cds) else 'white'
                ax_heatmap.barh(
                    i,
                    block_width,
                    left=block_pos['cds'],
                    height=0.7,
                    color=col_cds,
                    edgecolor='white',
                    linewidth=0.2
                )
            
            # rRNAs
            if bakta_available and rrna_norm is not None:
                rr = row['avg_rrna']
                col_rr = rrna_cmap(rrna_norm(rr)) if pd.notna(rr) else 'white'
                ax_heatmap.barh(
                    i,
                    block_width,
                    left=block_pos['rrna'],
                    height=0.7,
                    color=col_rr,
                    edgecolor='white',
                    linewidth=0.2
                )
            
            # Hypotheticals/CDS ratio
            if bakta_available and 'hypo_ratio' in block_pos and hypo_ratio_norm is not None:
                ratio = row['avg_hypo_cds_ratio']
                col_ratio = hypo_ratio_cmap(hypo_ratio_norm(ratio)) if pd.notna(ratio) else 'white'
                ax_heatmap.barh(
                    i,
                    block_width,
                    left=block_pos['hypo_ratio'],
                    height=0.7,
                    color=col_ratio,
                    edgecolor='white',
                    linewidth=0.2
                )

        # Labels over heatmap
        ax_heatmap.set_xlabel('')
        label_y = len(cluster_data) - 0.5 + 0.5  # Etwas über der obersten Zeile
        
        # Quality
        ax_heatmap.text(
            block_pos['quality'] + 0.5,
            label_y,
            'Quality',
            ha='center',
            va='bottom',
            fontsize=8,
            rotation=90
        )
        
        # QUAST labels
        if quast_available:
            ax_heatmap.text(
                block_pos['size'] + 0.5,
                label_y,
                'Genome\nsize',
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=90
            )
            ax_heatmap.text(
                block_pos['gc'] + 0.5,
                label_y,
                'GC %',
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=90
            )
        
        # Bakta labels
        if bakta_available:
            ax_heatmap.text(
                block_pos['cds'] + 0.5,
                label_y,
                'CDS',
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=90
            )
            ax_heatmap.text(
                block_pos['rrna'] + 0.5,
                label_y,
                'rRNA',
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=90
            )
            
            if 'hypo_ratio' in block_pos:
                ax_heatmap.text(
                    block_pos['hypo_ratio'] + 0.5,
                    label_y,
                    'Hypo/\nCDS',
                    ha='center',
                    va='bottom',
                    fontsize=8,
                    rotation=90
                )
                
        # ---- Legends and colorbars under heatmap ----
        pos = ax_heatmap.get_position()
        legend_y = pos.y0 - 0.15

        # hight, width of colorbars
        cbar_h = 0.10
        cbar_w = 0.015
        spacing = 0.05

        # Quality legend
        legend_handles = [
            Patch(color=qual_colors['High'],   label='High'),
            Patch(color=qual_colors['Medium'], label='Medium'),
            Patch(color=qual_colors['Low'],    label='Low'),
            Patch(color=qual_colors['Unknown'],label='Unknown'),
        ]

        ax_leg_quality = fig.add_axes([
            pos.x0,         # links
            legend_y,       # y
            0.12,           # Breite
            cbar_h          # Höhe
        ])

        ax_leg_quality.axis("off")
        ax_leg_quality.legend(
            handles=legend_handles,
            loc="center left",
            frameon=False,
            fontsize=7,
            ncol=1,
        )

        # collect available scalar metrics for colorbars
        cbar_specs = []

        if quast_available and size_norm is not None:
            cbar_specs.append(("size", size_norm, size_cmap, "Avg genome size (bp)"))
        if quast_available and gc_norm is not None:
            cbar_specs.append(("gc", gc_norm, gc_cmap, "Avg GC (%)"))
        if bakta_available and cds_norm is not None:
            cbar_specs.append(("cds", cds_norm, cds_cmap, "Avg CDS"))
        if bakta_available and rrna_norm is not None:
            cbar_specs.append(("rrna", rrna_norm, rrna_cmap, "Avg rRNAs"))
        if bakta_available and hypo_ratio_norm is not None:
            cbar_specs.append(
                ("hypo_ratio", hypo_ratio_norm, hypo_ratio_cmap, "Hypo/CDS ratio")
            )

        for idx, (_, norm, cmap, label) in enumerate(cbar_specs):
            ax_cbar = fig.add_axes(
                [
                    pos.x0 + 0.08 + idx * (cbar_w + spacing),
                    legend_y,
                    cbar_w,
                    cbar_h,
                ]
            )
            sm = cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=ax_cbar, orientation="vertical")
            cb.ax.tick_params(labelsize=6)
            cb.set_label(label, fontsize=7, labelpad=1)
    
    # ---- Info box ----
    # Calculate statistics
    total_mags = len(drep)
    total_clusters = drep['primary_cluster'].nunique()
    mags_in_top_clusters = drep[drep['primary_cluster'].isin(cluster_data['primary_cluster'])].shape[0]
    
    # Create info text
    info_text = (
        f"Total MAGs: {total_mags}\n"
        f"Dereplicated (representatives): {total_clusters}\n"
        f"Total clusters: {total_clusters}\n"
        f"MAGs in top {top_n} clusters: {mags_in_top_clusters} "
        f"({mags_in_top_clusters/total_mags*100:.1f}%)"
    )
    
    # Add text box in figure coordinates
    fig.text(
        0.02, - 0.05,
        info_text,
        fontsize=9,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, pad=0.5),
        family='monospace'
    )

    plt.suptitle(
        f'Top {top_n} species-level clusters by MAG count (by {tax_level})',
        fontsize=13,
        y=0.98,
        fontweight='bold'
    )
    
    out_file = os.path.join(output_path, f"drep_cluster_top{top_n}_{tax_level}.{fmt}")
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] dRep cluster plot saved to {out_file}")
    return out_file