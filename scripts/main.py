from json import load
from networkx import number_of_nodes
import pandas as pd
import numpy as np
import argparse
import time
import os
from version import __version__
from sanky_taxa import generate_taxa_sanky,taxa_sanky_rank
from comp_conta_plot import completeness_contamination_plot, rank_completeness_contamination_plot
from species_level_plot import species_level_plot
from mag_heatmap import mag_detection_heatmap
from heatmap import mag_heatmap
from histogram_plots import create_n50_histogram, number_of_contigs, create_assambly_info_histo
from rank_dist_plot import rank_distribution_pie
from amber_plots import binner_plot

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
        '-o',
        '--output',
        help="Path to the output folder to save the plots",
        dest="output",
        default=None
    )

    parser.add_argument(
        '-r',
        '--rank',
        help="Select which rank should be used for the rank sanky",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"]
    )

    parser.add_argument(
        '-n',
        '--top_n_counts',
        help='For the rank distrubionen, only plot the top n. Min: 5 DEFAULT:5',
        type=positive_int,
        dest="n",
        default=5
    )

    parser.add_argument(
        '--amber',
        help="Input CAMI amber result file for different plots",
        dest="amber_file",
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
        '--test'
    )

    parser.print_usage = parser.print_help

    args = parser.parse_args()

    return args

def load_dfs(coverm, checkm, checkm2, gtdb, drep):
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

    for i, file in enumerate(os.listdir(coverm)):
        path = os.path.join(coverm, file)
        if file.endswith(".csv"):
            coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, index_col=0)
        elif checkm.endswith(".tsv") or path.endswith(".tabular"):
            coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, sep="\t", index_col=0)
        else:
            coverm_dfs[f'coverm_{i}'] = pd.read_csv(path, index_col=0)
        print(f"[INFO] coverm_{i} loaded: {coverm_dfs[f'coverm_{i}'].shape} rows x columns")
    
    dfs['coverm'] = coverm_dfs

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
    clean_dfs = []
    if len(coverm_dfs.keys()) == 1:
        return coverm_dfs['coverm_0']
    else:
        dfs = [coverm_dfs[k] for k in coverm_dfs.keys()]

        for df in dfs:
            clean_dfs.append(df)

        coverm_merged = pd.concat(clean_dfs, axis=1)

        print(f"[INFO] gtdb loaded: {coverm_merged.shape} rows x columns")

        return coverm_merged

def check_path(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"[INFO] Created folder: {output_path}")
    else:
        print(f"[INFO] Folder already exists: {output_path}")

if __name__ == '__main__':
    start_time = time.time()
    args = parse_arguments()

    dfs = load_dfs(args.coverm_path, args.checkm_file, args.checkm2_file, args.gtdb_file, args.drep_file)

    dfs['coverm'] = merged_coverm(dfs['coverm'])

    check_path(args.output)

    generate_taxa_sanky(dfs['gtdb'], args.output)
    taxa_sanky_rank(dfs['gtdb'], args.output, args.rank)

    completeness_contamination_plot(dfs['checkm'], args.output)

    species_level_plot(dfs['drep'], args.output)

    mag_detection_heatmap(dfs["coverm"], args.output)

    mag_heatmap(dfs["coverm"], dfs["gtdb"], args.output)

    create_n50_histogram(dfs['checkm2'], args.output)
    number_of_contigs(dfs["checkm2"], args.output)
    create_assambly_info_histo(dfs["checkm2"], args.output)

    rank_distribution_pie(dfs["gtdb"], args.output, args.rank, args.n)

    if args.amber_file is not None:
        binner_plot(load_single_df(args.amber_file), args.output)

    if args.gtdb_ar_file is not None and args.gtdb_bac_file is not None:
        rank_completeness_contamination_plot(load_single_df(args.test), load_single_df(args.gtdb_bac_file), load_single_df(args.gtdb_ar_file), args.rank, args.output, args.n)

    end_time = time.time()
    print(f'[INFO] Run time: {time.strftime("%H:%M:%S", time.gmtime(end_time - start_time))}')