from json import load
import pandas as pd
import numpy as np
import argparse
import time
import os
import csv
from .version import __version__
from .sanky_taxa import generate_taxa_sanky,taxa_sanky_rank
from .comp_conta_plot import completeness_contamination_plot, rank_completeness_contamination_plot, metadata_completeness_contamination_plot
from .species_level_plot import species_level_plot
from .mag_heatmap import mag_detection_heatmap
from .heatmap import mag_heatmap
from .histogram_plots import create_n50_histogram, number_of_contigs, create_assambly_info_histo
from .rank_dist_plot import rank_distribution_pie
# from .scripts.amber_plots import binner_plot
from .assembly_quality import assembly_quality_plot
from .bakta_checkm2_plot import bakta_annotation_plot
# from drep_taxa_plot import drep_taxa_plot
from .drep_cluster_plot import drep_cluster_plot


def read_table(path, index_col=None, prefer_tsv=None):
    if path is None:
        return None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()

    # If caller explicitly prefers TSV, respect it
    if prefer_tsv is True:
        df = pd.read_csv(path, sep="\t", engine="python")
    elif prefer_tsv is False:
        # Prefer CSV, fallback to sniffing
        try:
            df = pd.read_csv(path, sep=",", engine="python")
        except Exception:
            df = pd.read_csv(path, sep=None, engine="python")
    else:
        # Auto mode: tabs -> TSV else CSV/sniff
        if "\t" in first:
            df = pd.read_csv(path, sep="\t", engine="python")
        else:
            try:
                df = pd.read_csv(path, sep=",", engine="python")
            except Exception:
                df = pd.read_csv(path, sep=None, engine="python")

    # Only set index if requested
    if index_col is not None and df.shape[1] > index_col:
        df = df.set_index(df.columns[index_col])

    return df



# ---- Argument parsing ---- #
def positive_int(value):
    ivalue = int(value)
    if ivalue < 5:
        raise argparse.ArgumentTypeError(f"{value} must be >= 5")
    return ivalue

def parse_arguments(argv=None):

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
        '--top_n',
        help='For the rank distribution, only plot the top n. Min: 5 DEFAULT:5',
        type=positive_int,
        dest="n",
        default=5
    )

    parser.add_argument(
        '--test'
    )

    # Heatmap Layout Options
    heatmap_layout = parser.add_argument_group('Heatmap Layout Options',
                                                'Fine-tune the heatmap visualization layout')
    
    heatmap_layout.add_argument(
        '--fig_height_per_row',
        type=float,
        default=0.50,
        dest='fig_height_per_row',
        help='Height per sample row in heatmap (default: 0.50). Increase for many samples (e.g., 0.55-0.60)'
    )
    
    heatmap_layout.add_argument(
        '--fig_base_height',
        type=float,
        default=2.5,
        dest='fig_base_height',
        help='Base height for top/bottom margins (default: 2.5). Decrease for less space above (e.g., 2.0)'
    )
    
    heatmap_layout.add_argument(
        '--top_bar_height',
        type=float,
        default=0.7,
        dest='top_bar_height',
        help='Relative height of top histogram (default: 0.7). Smaller = more compact (e.g., 0.5-0.6)'
    )
    
    heatmap_layout.add_argument(
        '--hspace',
        type=float,
        default=0.25,
        dest='hspace',
        help='Vertical spacing between top and heatmap (default: 0.20). Reduce if too much space (0.15)'
    )
    
    heatmap_layout.add_argument(
        '--heatmap_width',
        type=float,
        default=11.0,
        dest='heatmap_width',
        help='Width ratio of the heatmap (default: 11.0). Increase for larger heatmap (e.g., 13.0)'
    )

    heatmap_layout.add_argument(
        '--max_col',
        type=int,
        default=10,
        dest='max_col',
        help="Max shown taxonomy names."
    )
    
    heatmap_layout.add_argument(
        '--meta_bar_thickness',
        type=float,
        default=0.8,
        dest='meta_bar_thickness',
        help='Thickness of metadata bars (default: 0.7). Smaller = more space between bars (e.g., 0.5-0.6)'
    )
    
    heatmap_layout.add_argument(
        '--title_y',
        type=float,
        default=0.95,
        dest='title_y',
        help='Vertical position of title, 0.9-1.0 (default: 0.95). Lower = more space from top'
    )
    
    heatmap_layout.add_argument(
        '--heatmap_preset',
        choices=['default', 'compact', 'spacious', 'many_samples', 'few_samples', 'focus_heatmap'],
        default='default',
        dest='heatmap_preset',
        help='Use a predefined layout preset for common scenarios'
    )

    heatmap_layout.add_argument(
        '--spacer_legend',
        type=float,
        default=0.3,
        dest='spacer_legend',
        help='Spacer between Legend and meta_bars.'
    )

    heatmap_layout.add_argument(
        '--spacer_meta',
        type=float,
        default=2.0,
        dest='spacer_meta',
        help='Spacer between Meta_bars and Heatmap.'
    )

    heatmap_layout.add_argument(
        '--spacer_heatmap',
        type=float,
        default=0.10,
        dest='spacer_heatmap',
        help='Spacer between Heatmap and histogram.'
    )

    heatmap_layout.add_argument(
        '--legend',
        type=float,
        default=2.5,
        dest='legend',
        help='Legend width'
    )

    heatmap_layout.add_argument(
        '--meta_bar_add',
        type=float,
        default=1.5,
        dest='meta_bar_add',
        help='Adding width to meta_bars'
    )

    heatmap_layout.add_argument(
        '--top_bar_spacer',
        type=float,
        default=0.0,
        dest='top_bar_spacer',
        help='Space between top histogram title.'
    )

    heatmap_layout.add_argument(
        '--no_log',
        action="store_true",
        help="Disable log10 for MAGs/rank bar plot."
    )

    # drep
    parser.add_argument(
        '--tax_levels_space',
        type=float,
        default=0.3,
        dest='tax_levels_space',
        help='Space for tax_levels columns in drep_cluster_plot.'
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
        '--tax_levels',
        nargs='+',
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        default=["phylum", "genus"],
        help="GTDB taxonomy levels to show in dRep cluster plot"
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

    args = parser.parse_args(argv)

    return args

# ---- Data loading ---- #
def load_dfs(coverm, checkm, checkm2, gtdb, drep, bakta=None, quast=None, metadata=None):
    dfs = {}

    # Required
    dfs["checkm"]  = read_table(checkm,  index_col=0)
    print(f"[INFO] checkm loaded:  {dfs['checkm'].shape}")

    dfs["checkm2"] = read_table(checkm2, index_col=0)
    print(f"[INFO] checkm2 loaded: {dfs['checkm2'].shape}")

    dfs["drep"] = read_table(drep, index_col=None, prefer_tsv=False)
    print(f"[INFO] drep loaded:   {dfs['drep'].shape}")


    dfs["gtdb"] = read_table(gtdb, index_col=0, prefer_tsv=True)
    print(f"[INFO] gtdb loaded:   {dfs['gtdb'].shape}")

    if dfs["drep"] is not None:
        cols = [c.strip() for c in dfs["drep"].columns.astype(str)]
        dfs["drep"].columns = cols

    # --- dRep: ensure genome column exists ---
    if dfs["drep"] is not None:
        dfs["drep"].columns = [str(c).strip() for c in dfs["drep"].columns]

        if "genome" not in dfs["drep"].columns and "Genome" in dfs["drep"].columns:
            dfs["drep"] = dfs["drep"].rename(columns={"Genome": "genome"})

        if "genome" not in dfs["drep"].columns:
            dfs["drep"] = dfs["drep"].reset_index().rename(columns={"index": "genome"})


    # Optional: coverm 
    coverm_dfs = {}
    if coverm is None:
        print("[INFO] No coverm provided.")
    else:
        if os.path.isdir(coverm):
            files = [f for f in os.listdir(coverm) if not f.startswith(".")]
            if not files:
                print(f"[WARN] No files found in coverm directory: {coverm}")
            for i, file in enumerate(sorted(files)):
                path = os.path.join(coverm, file)
                coverm_dfs[f"coverm_{i}"] = read_table(path, index_col=0)
                print(f"[INFO] coverm_{i} loaded: {coverm_dfs[f'coverm_{i}'].shape}")
        else:
            coverm_dfs["coverm_0"] = read_table(coverm, index_col=0)
            print(f"[INFO] coverm_0 loaded: {coverm_dfs['coverm_0'].shape}")

    dfs["coverm"] = coverm_dfs

    # Optional tables (only load + print if provided)
    if bakta is not None:
        dfs["bakta"] = read_table(bakta, index_col=0)
        print(f"[INFO] bakta loaded:  {dfs['bakta'].shape}")

    if quast is not None:
        dfs["quast"] = read_table(quast, index_col=0)
        print(f"[INFO] quast loaded:  {dfs['quast'].shape}")

    if metadata is not None:
        dfs["metadata"] = read_table(metadata, index_col=0)
        print(f"[INFO] metadata loaded:{dfs['metadata'].shape}")
    else:
        print("[INFO] No metadata file provided.")

    return dfs


def load_single_df(file_path):
    return read_table(file_path, index_col=0)

    
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
def main(argv=None):
    start_time = time.time()
    args = parse_arguments(argv)

    dfs = load_dfs(
        args.coverm_path,
        args.checkm_file,
        args.checkm2_file,
        args.gtdb_file,
        args.drep_file,
        bakta=args.bakta_file,
        quast=args.quast_file,
        metadata=args.metadata_file,
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
        user_metrics = [m.strip().lower() for m in args.bakta_metrics]

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
    # species_level_plot(dfs['drep'], args.output)
    

    # ---- dRep cluster plot ----
    drep_cluster_plot(
        dfs["drep"],
        dfs["gtdb"],
        args.output,
        tax_levels=args.tax_levels,
        top_n=args.n,
        fig_size=comp_fig_size,
        fmt=args.format,
        tax_levels_space=args.tax_levels_space,
        checkm2_df=dfs.get("checkm2"),
        quast_df=dfs.get("quast"),
        bakta_df=dfs.get("bakta"),
    )

    # mag_detection_heatmap(dfs["coverm"], args.output)
    
    mag_heatmap(
        dfs["coverm"],
        dfs["gtdb"],
        args.output,
        rank=tax_rank,
        metadata_df=dfs.get("metadata"),
        meta_cols=args.meta_cols,
        meta_bin_width=args.meta_bin_width,
        fmt=args.format,
        # Layout parameters
        top_bar_height=args.top_bar_height,
        hspace=args.hspace,
        heatmap_width=args.heatmap_width,
        spacer_legend=args.spacer_legend,
        spacer_meta=args.spacer_meta,
        spacer_heatmap=args.spacer_heatmap,
        legend=args.legend,
        meta_bar_add=args.meta_bar_add,
        top_bar_spacer=args.top_bar_spacer,
        max_col=args.max_col,
        log_top=not args.no_log,
    )

    # ---- Assembly histograms from checkm2 ----
    # create_n50_histogram(dfs['checkm2'], args.output)
    # number_of_contigs(dfs["checkm2"], args.output)
    # create_assambly_info_histo(dfs["checkm2"], args.output)

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
    # rank_distribution_pie(dfs["gtdb"], args.output, args.rank, args.n)

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
    # if args.amber_file is not None:
    #     binner_plot(load_single_df(args.amber_file), args.output)


    end_time = time.time()
    print(f'[INFO] Run time: {time.strftime("%H:%M:%S", time.gmtime(end_time - start_time))}')

    return 0

if __name__ == '__main__':
    raise SystemExit(main())
