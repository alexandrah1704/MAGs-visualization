from json import load
from networkx import number_of_nodes
import pandas as pd
import numpy as np
import argparse
import time
import os
from version import __version__
from sanky_taxa import generate_taxa_sanky,taxa_sanky_rank
from comp_conta_plot import completeness_contamination_plot, rank_completeness_contamination_plot, metadata_completeness_contamination_plot
from species_level_plot import species_level_plot
from mag_heatmap import mag_detection_heatmap
from heatmap import mag_heatmap
from histogram_plots import create_n50_histogram, number_of_contigs, create_assambly_info_histo
from rank_dist_plot import rank_distribution_pie
from amber_plots import binner_plot
from assembly_quality import assembly_quality_plot
from bakta_checkm2_plot import bakta_annotation_plot
# from drep_taxa_plot import drep_taxa_plot
from drep_cluster_plot import drep_cluster_plot

# ---- Argument parsing ---- #
def positive_int(value):
    ivalue = int(value)
    if ivalue < 5:
        raise argparse.ArgumentTypeError(f"{value} must be >= 5")
    return ivalue

def parse_arguments():

    parser = argparse.ArgumentParser(
        prog="MAGs-visualization",
        description="TODO",
        usage="TODO",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True,
    )

    parser.add_argument(
        '-v',
        '-V',
        '--version',
        action='version',
        version=__version__
    )

    # Input
    parser.add_argument(
        '--coverm',
        help="Path to the CoverM files",
        dest="coverm_path"
    )

    parser.add_argument(
        '--checkm',
        help="Path to the CheckM file",
        dest="checkm_file",
        default=None
    )

    parser.add_argument(
        '--checkm2',
        help="Path to the CheckM2 file",
        dest="checkm2_file",
        default=None
    )

    parser.add_argument(
        '--gtdb',
        help="Path to the GTDB file",
        dest="gtdb_file",
        default=None
    )

    parser.add_argument(
        '--drep',
        help='Path to the dRep file',
        dest='drep_file',
        default=None
    )

    parser.add_argument(
        '--gtdb_bac',
        help="Input GTDB bacteria result file",
        dest='gtdb_bac_file',
        default=None
    )

    parser.add_argument(
        '--gtdb_ar',
        help="Input GTDB archaea result file",
        dest='gtdb_ar_file',
        default=None
    )

    parser.add_argument(
        '--amber',
        help="Input CAMI amber result file for different plots",
        dest="amber_file",
        default=None
    )

    parser.add_argument(
        '--quast',
        help="'Path to QUAST table",
        dest="quast_file",
        default=None,
    )

    parser.add_argument(
        '--bakta',
        help='Path to Bakta annotation table',
        dest='bakta_file',
        default=None
    )

    parser.add_argument(
        '--metadata',
        help='Path to METADATA',
        dest="metadata_file",
        default=None,
    )

    parser.add_argument(
        '--metadata_heatmap',
        help='Path to METADATA_HEATMAP',
        dest='metadata_heatmap_file',
        default=None
    )

    parser.add_argument(
        '--metadata_heatmap_new',
        help="Path to METADATA_HEATMAP_NEW",
        dest="metadata_heatmap_new_file",
        default=None,
    )

    # Output
    parser.add_argument(
        '-o',
        '--output',
        help="Path to the output folder to save the plots",
        dest="output",
        default=None
    )

    # Taxonomic rank
    parser.add_argument(
        '-r',
        '--rank',
        help="Select which rank should be used for the rank sanky",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        default="phylum",
    )

    parser.add_argument(
        '-n',
        '--top_n_counts',
        help='For the rank distribution, only plot the top n. Min: 5 DEFAULT:5',
        type=positive_int,
        dest="n",
        default=5
    )

    parser.add_argument(
        '--test'
    )

    # Plot arguments and styles
    parser.add_argument(
        '--quality',
        action="store_true",
        help="Create plots colored by quality catergories"
    )

    parser.add_argument(
        '--tax',
        action="store_true",
        help="Create plots colored by taxonomy"
    )

    parser.add_argument(
        '--fig_size',
        nargs=2,
        type=float,
        metavar=('WIDTH', 'HEIGHT'),
        help="Figure size in inches for plots",
        dest='fig_size'
    )

    parser.add_argument(
        '--tax_level',
        help="GTDB rank used for taxonomy-colored plots",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        dest="tax_level"
    )

    parser.add_argument(
        '--format',
        help="Output image format for plots",
        choices=["png", "pdf", "svg"],
        default="png",
        dest="format"
    )

    # Metadata options
    parser.add_argument(
        '--meta_col',
        help="Metadata-Column to color by (e.g. 'weather, 'temperatur')",
        dest="meta_col",
        default=None
    )

    parser.add_argument(
        '--meta_cols', 
        nargs='+', 
        help='List of metadata columns'
        )

    parser.add_argument(
        '--meta_bin_width',
        help='Bin-Width for numeric metadata',
        dest='meta_bin_width',
        type=float,
        default=5.0
    )

    # Quast columns
    parser.add_argument(
        '--column_choice',
        nargs='+',
        metavar='METRIC',
        help=("QUAST metrics to plot in the assembly quality dashboard. "
            "Use the exact QUAST names, e.g. "
            "\"N50\" \"GC\" \"# contigs\" \"Total length (>= 0 bp)\""
        )
    )

    parser.add_argument(
        '--color_by',
        choices=['quality', 'tax', 'meta'],
        help=(
            "How to color points in the assembly quality plots: "
            "'quality' = High/Medium/Low quality, "
            "'tax' = taxonomy (uses GTDB) "
            "'meta' = metadata. "
            "If omitted, points are colored by completeness bins."
        )
    )

    # Bakta metric
    parser.add_argument(
        '--bakta_metrics',
        nargs='+',
        metavar='METRIC',
        default=None,
        help="Use the exact BAKTA names ('CDSs', 'hypotheticals'...)"
    )

    parser.add_argument(
        "--ratio",
        action="store_true",
        help="Plot with feature/cds ratios in bakta plot"
    )

    parser.print_usage = parser.print_help

    args = parser.parse_args()

    return args

# ---- Data loading ---- #
def load_dfs(coverm, checkm, checkm2, gtdb, drep, bakta=None, quast=None, 
             metadata=None, metadata_heatmap=None, metadata_heatmap_new=None):
    dfs = {}
    coverm_dfs = {}

    if checkm.endswith(".csv"):
        dfs['checkm'] = pd.read_csv(checkm, index_col=0)
    elif checkm.endswith(".tsv") or checkm.endswith(".tabular"):
        dfs['checkm'] = pd.read_csv(checkm, sep="\t", index_col=0)
    else:
        dfs['checkm'] = pd.read_csv(checkm, index_col=0)
    print(f"[INFO] checkm loaded: {dfs['checkm'].shape} rows x columns")
    
    if checkm2.endswith(".csv"):
        dfs['checkm2'] = pd.read_csv(checkm2, index_col=0)
    elif checkm2.endswith(".tsv") or checkm2.endswith(".tabular"):
        dfs['checkm2'] = pd.read_csv(checkm2, sep="\t", index_col=0)
    else:
        dfs['checkm2'] = pd.read_csv(checkm2, index_col=0)
    print(f"[INFO] checkm2 loaded: {dfs['checkm2'].shape} rows x columns")
    
    if drep.endswith(".csv"):
        dfs['drep'] = pd.read_csv(drep, index_col=0)
    elif drep.endswith(".tsv") or drep.endswith(".tabular"):
        dfs['drep'] = pd.read_csv(drep, sep="\t", index_col=0)
    else:
        dfs['drep'] = pd.read_csv(drep, index_col=0)
    print(f"[INFO] drep loaded: {dfs['drep'].shape} rows x columns")
    
    if gtdb.endswith(".csv"):
        dfs['gtdb'] = pd.read_csv(gtdb, index_col=0)
    elif gtdb.endswith(".tsv") or gtdb.endswith(".tabular"):
        dfs['gtdb'] = pd.read_csv(gtdb, sep="\t", index_col=0)
    else:
        dfs['gtdb'] = pd.read_csv(gtdb, index_col=0)
    print(f"[INFO] gtdb loaded: {dfs['gtdb'].shape} rows x columns")
    coverm_dfs = {}

    if os.path.isdir(coverm):
        files = [f for f in os.listdir(coverm) if not f.startswith('.')]
        if not files:
            print(f"[WARN] No files found in coverm directory: {coverm}")
        for i, file in enumerate(sorted(files)):
            path = os.path.join(coverm, file)
            if path.endswith(".csv"):
                coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, index_col=0)
            elif path.endswith(".tsv") or path.endswith(".tabular"):
                coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, sep="\t", index_col=0)
            else:
                coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, index_col=0)
            print(f"[INFO] coverm_{i} loaded: {coverm_dfs[f'coverm_{i}'].shape} rows x columns")
    else:
        path = coverm
        if path.endswith(".csv"):
            df_cov = pd.read_csv(path, index_col=0)
        elif path.endswith(".tsv") or path.endswith(".tabular"):
            df_cov = pd.read_csv(path, sep="\t", index_col=0)
        else:
            df_cov = pd.read_csv(path, index_col=0)
        coverm_dfs['coverm_0'] = df_cov
        print(f"[INFO] coverm_0 loaded (single file): {df_cov.shape} rows x columns")

    dfs['coverm'] = coverm_dfs

    if bakta is not None:
        if bakta.endswith(".csv"):
            dfs["bakta"] = pd.read_csv(bakta, index_col=0)
        else:
            dfs["bakta"] = pd.read_csv(bakta, sep="\t", index_col=0)
        print(f"[INFO] bakta loaded: {dfs['bakta'].shape} rows x columns")

    if quast is not None:
        if quast.endswith(".csv"):
            dfs["quast"] = pd.read_csv(quast, index_col=0)
        else:
            dfs["quast"] = pd.read_csv(quast, sep="\t", index_col=0)
        print(f"[INFO] quast loaded: {dfs['quast'].shape} rows x columns")
    
    if metadata is not None:
        if metadata.endswith(".csv"):
            dfs['metadata'] = pd.read_csv(metadata, index_col=0)
        elif metadata.endswith(".tsv") or metadata.endswith(".tabular"):
            dfs['metadata'] = pd.read_csv(metadata, sep="\t", index_col=0)
        else:
            dfs['metadata'] = pd.read_csv(metadata, index_col=0)
        print(f"[INFO] metadata loaded: {dfs['metadata'].shape} rows x columns")
    else:
        print("[INFO] No metadata file provided.")

    if metadata_heatmap is not None:
        if metadata_heatmap.endswith(".csv"):
            dfs['metadata_heatmap'] = pd.read_csv(metadata_heatmap, index_col=0)
        elif metadata_heatmap.endswith(".tsv") or metadata_heatmap.endswith(".tabular"):
            dfs['metadata_heatmap'] = pd.read_csv(metadata_heatmap, sep="\t", index_col=0)
        else:
            dfs['metadata_heatmap'] = pd.read_csv(metadata_heatmap, index_col=0)
        print(f"[INFO] metadata_heatmap loaded: {dfs['metadata_heatmap'].shape} rows x columns")
        print(f"[INFO] example samples in metadata_heatmap: {dfs['metadata_heatmap'].index[:5].tolist()}")
    else:
        print("[INFO] No heatmap metadata file provided.")
    
    if metadata_heatmap_new is not None:
        if metadata_heatmap_new.endswith(".csv"):
            dfs['metadata_heatmap_new'] = pd.read_csv(metadata_heatmap_new, index_col=0)
        elif metadata_heatmap_new.endswith(".tsv") or metadata_heatmap_new.endswith(".tabular"):
            dfs['metadata_heatmap_new'] = pd.read_csv(metadata_heatmap_new, sep="\t", index_col=0)
        else:
            dfs['metadata_heatmap_new'] = pd.read_csv(metadata_heatmap_new, index_col=0)
        print(f"[INFO] metadata_heatmap_new loaded: {dfs['metadata_heatmap_new'].shape} rows x columns")
        print(f"[INFO] example samples in metadata_heatmap_new: {dfs['metadata_heatmap_new'].index[:5].tolist()}")
    else:
        print("[INFO] No heatmap metadata file provided.")

    return dfs


def load_single_df(file_path):
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path, index_col=0)
        print(f"[INFO] {file_path} loaded: {df.shape} rows x columns")
        return df
    elif file_path.endswith('.tsv') or file_path.endswith('.tabular'):
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        print(f"[INFO] {file_path} loaded: {df.shape} rows x columns")
        return df
    else:
        df = pd.read_csv(file_path, index_col=0)
        print(f"[INFO] {file_path} loaded: {df.shape} rows x columns")
        return df

    
def merged_coverm(coverm_dfs):
    """
    Merge multiple coverm tables column wise.
    """
    if not coverm_dfs:
        print("[WARN] No CoverM DataFrames to merge.")
        return pd.DataFrame()

    if len(coverm_dfs) == 1:
        df = coverm_dfs.get('coverm_0', next(iter(coverm_dfs.values())))
        return df

    dfs = [coverm_dfs[k] for k in sorted(coverm_dfs.keys())]
    coverm_merged = pd.concat(dfs, axis=1)
    print(f"[INFO] coverm merged: {coverm_merged.shape} rows x columns")
    return coverm_merged


def check_path(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"[INFO] Created folder: {output_path}")
    else:
        print(f"[INFO] Folder already exists: {output_path}")

# ---- Main ---- #
if __name__ == '__main__':
    start_time = time.time()
    args = parse_arguments()

    dfs = load_dfs(
        args.coverm_path,
        args.checkm_file,
        args.checkm2_file,
        args.gtdb_file,
        args.drep_file,
        bakta=args.bakta_file,
        quast=args.quast_file,
        metadata=args.metadata_file,
        metadata_heatmap=args.metadata_heatmap_file,
        metadata_heatmap_new=args.metadata_heatmap_new_file,
    )

    dfs['coverm'] = merged_coverm(dfs['coverm'])
    check_path(args.output)

    # ---- Figure size ----
    if args.fig_size is not None:
        comp_fig_size = tuple(args.fig_size)
    else:
        comp_fig_size = None

    # ---- Taxonomic rank for various plots (heatmap, assembly_quality...)
    if args.tax_level is not None:
        tax_rank = args.tax_level
    elif args.rank is not None:
        tax_rank = args.rank
    else:
        tax_rank = "phylum"

    # same
    if args.tax_level is not None:
        aq_tax_rank = args.tax_level
    elif args.rank is not None:
        aq_tax_rank = args.rank
    else:
        aq_tax_rank = "phylum"


    # ---- Batka metric selection ----
    if args.bakta_metrics:
        # Eingaben wie "CDS", "cds", "CDs" -> alles nach lowercase
        user_metrics = [m.strip().lower() for m in args.bakta_metrics]

        # canonical_name -> (Display-Name, Spaltenname in bakta)
        bakta_metric_map = {
            "cds": ("CDS", "cdss"),
            "hypotheticals": ("Hypotheticals", "hypotheticals"),
            "rrnas": ("rRNAs", "rrnas"),
            "trnas": ("tRNAs", "trnas"),
            "crispr": ("CRISPR", "crispr arrays"),
            "pseudogenes": ("Pseudogenes", "pseudogenes"),
        }

        metrics_dict = {}
        for m in user_metrics:
            if m in bakta_metric_map:
                disp, col = bakta_metric_map[m]
                metrics_dict[disp] = col

        if not metrics_dict:
            print(
                "[WARN] None of the requested Bakta metrics matched. "
                "Falling back to default metrics (CDS, Hypotheticals, rRNAs, tRNAs)."
            )
            metrics_dict = None
    else:
        metrics_dict = None

    bakta_ratio_flag = args.ratio

    # ---- Sankey plots ----
    generate_taxa_sanky(dfs['gtdb'], args.output, args.rank)
    taxa_sanky_rank(dfs['gtdb'], args.output, args.rank)

    # Decide which comp_conta_plot to run
    run_quality = args.quality or (not args.quality and not args.tax)
    run_tax     = args.tax or (not args.quality and not args.tax)

    # ---- Completeness-Contamination Plots ----
    # ---- Quality completeness/contamination plots ----
    if run_quality:
        completeness_contamination_plot(
            dfs['checkm'],
            args.output,
            tag="checkm",
            title="CheckM: Completeness vs Contamination",
            fig_size=comp_fig_size if comp_fig_size is not None else (9, 8),
            fmt=args.format
        )
        completeness_contamination_plot(
            dfs['checkm2'],
            args.output,
            tag="checkm2",
            title="CheckM2: Completeness vs Contamination",
            fig_size=comp_fig_size if comp_fig_size is not None else (9, 8),
            fmt=args.format
        )
    
    # ---- Taxonomic completeness/contamination plots ----
    if run_tax and args.gtdb_file is not None:
        print("[RUN] Rank plot (using single GTDB table for both ar/bac) ...")
        rank_completeness_contamination_plot(
            dfs["checkm"],
            dfs["checkm2"],
            dfs["gtdb"],
            tax_rank,
            args.output,
            args.n,
            fig_size=comp_fig_size if comp_fig_size is not None else (10, 8),
            fmt=args.format
        )
    
    # --- Metadata completeness/contamination plots ---
    if args.metadata_file is not None and "metadata" in dfs and args.meta_col is not None:
        if args.meta_col not in dfs["metadata"].columns:
            print(f"[WARN] Metadata column '{args.meta_col}' not found in metadata table. "
                  f"Available columns: {list(dfs['metadata'].columns)}")
        else:
            print(f"[RUN] Metadata plot colored by '{args.meta_col}' ...")
            metadata_completeness_contamination_plot(
                dfs["checkm"],
                dfs["checkm2"],
                dfs["metadata"],
                meta_col=args.meta_col,
                output_path=args.output,
                fig_size=comp_fig_size if comp_fig_size is not None else (10, 8),
                fmt=args.format,
                bin_width=args.meta_bin_width,
                top_n=args.n
            )
    else:
        if args.meta_col is not None and args.metadata_file is None:
            print("[WARN] --meta_col set but no --metadata file was given."
                  "→ skipping metadata-based plots.")

    # drep and heatmap
    species_level_plot(dfs['drep'], args.output)
    

    # ---- dRep cluster plot ----
    drep_cluster_plot(
        dfs["drep"],
        dfs["gtdb"],
        args.output,
        tax_level=tax_rank,
        top_n=args.n,
        fig_size=comp_fig_size,
        fmt=args.format,
        checkm2_df=dfs.get("checkm2"),
        quast_df=dfs.get("quast"),
        bakta_df=dfs.get("bakta"),
    )

    mag_detection_heatmap(dfs["coverm"], args.output)

    # ---- Heatmap with bars ----
    mag_heatmap(
        dfs["coverm"],
        dfs["gtdb"],
        args.output,
        rank=tax_rank,
        metadata_df=dfs.get("metadata"),
        meta_cols=args.meta_cols,
        meta_bin_width=args.meta_bin_width,
        fmt=args.format,
    )

    # ---- Assembly histograms from checkm2 ----
    create_n50_histogram(dfs['checkm2'], args.output)
    number_of_contigs(dfs["checkm2"], args.output)
    create_assambly_info_histo(dfs["checkm2"], args.output)

    # ---- Quast-based plots ----
    if "quast" in dfs:
        print("[INFO] QUAST found → creating assembly quality plots.")
        assembly_quality_plot(
            dfs,
            args.output,
            metrics=args.column_choice, 
            color_by=args.color_by,
            tax_rank=aq_tax_rank,
            meta_col=args.meta_col,
            fig_size=comp_fig_size if comp_fig_size is not None else (5, 5),
            fmt=args.format,
        )
    else:
        print("[WARN] No Quast file loaded → skipping assembly quality plots.")

    # ---- Rank distribution pie plot ----
    rank_distribution_pie(dfs["gtdb"], args.output, args.rank, args.n)

    # ---- Bakta-based plots ----
    if "bakta" in dfs:
        print("[INFO] Bakta found → creating Bakta plots.")
        bakta_annotation_plot(dfs,
                                   args.output,
                                   metrics=metrics_dict,
                                   color_by=args.color_by,
                                   tax_rank=aq_tax_rank,
                                   meta_col=args.meta_col,
                                   plot_ratio=args.ratio,
                                   )
    else:
        print("[WARN] No Bakta file loaded → skipping Bakta plots.")

   # ---- Amber plots ----
    if args.amber_file is not None:
        binner_plot(load_single_df(args.amber_file), args.output)


    end_time = time.time()
    print(f'[INFO] Run time: {time.strftime("%H:%M:%S", time.gmtime(end_time - start_time))}')