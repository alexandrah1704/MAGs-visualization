import os
import time
import pandas as pd
from main import (
    load_dfs,
    merged_coverm,
    check_path,
    completeness_contamination_plot,
    rank_completeness_contamination_plot,
    species_level_plot,
    mag_detection_heatmap,
    mag_heatmap,
    create_n50_histogram,
    number_of_contigs,
    create_assambly_info_histo,
    rank_distribution_pie,
    generate_taxa_sanky,
    taxa_sanky_rank,
    # box_swarm_plot,
)

INPUT_DIR = "test-data"
COVERM_DIR = "test-data/coverm.tsv"
CHECKM = "test-data/checkm.tsv"
CHECKM2 = "test-data/checkm2.tsv"
GTDB = "test-data/gtdb.tsv"
DREP = "test-data/drep.csv"
BAKTA = "test-data/bakta.tsv"
CAMI_AMBER = "test-data/cami_amber.tsv"
QUAST = "test-data/quast.tsv"
OUTPUT_DIR = "test-plots"
RANK_INPUT = input("Select taxonomic rank (d/p/c/o/f/g/s): ").strip().lower()
TOP_N = 10

coverm = pd.read_csv(os.path.join(INPUT_DIR, "coverm.tsv"), sep="\t")
gtdb   = pd.read_csv(os.path.join(INPUT_DIR, "gtdb.tsv"),   sep="\t")

rank_map = {
    "d": "domain", "p": "phylum", "c": "class", "o": "order",
    "f": "family", "g": "genus", "s": "species"
}
RANK = rank_map.get(RANK_INPUT, "phylum")
print(f"-> Using rank: {RANK}")
selected_rank = rank_map.get(RANK, "phylum")

if __name__ == "__main__":
    start_time = time.time()
    print("[RUN] Test script started...")

    # ---- Load Data ----
    dfs = load_dfs(COVERM_DIR, CHECKM, CHECKM2, GTDB, DREP, BAKTA, QUAST)
    dfs["coverm"] = merged_coverm(dfs["coverm"])


    check_path(OUTPUT_DIR)

    # ---- Execute all plots ----
    completeness_contamination_plot(dfs["checkm"], OUTPUT_DIR, tag="checkm", title="CheckM: Completeness vs Contamination")
    completeness_contamination_plot(dfs["checkm2"], OUTPUT_DIR, tag="checkm2", title="CheckM2: Completeness vs Contamination")

    rank_completeness_contamination_plot(
      dfs["checkm"],
      dfs["checkm2"],
        # dfs["gtdb"],
        dfs["gtdb"],
        RANK,
        OUTPUT_DIR,
        TOP_N  
    )

    species_level_plot(dfs["drep"], OUTPUT_DIR)
    mag_detection_heatmap(dfs["coverm"], OUTPUT_DIR)
    mag_heatmap(coverm, gtdb, OUTPUT_DIR, rank=RANK)
    create_n50_histogram(dfs["checkm2"], OUTPUT_DIR)
    number_of_contigs(dfs["checkm2"], OUTPUT_DIR)
    create_assambly_info_histo(dfs["checkm2"], OUTPUT_DIR)
    # box_swarm_plot(dfs, OUTPUT_DIR)
    rank_distribution_pie(dfs["gtdb"], OUTPUT_DIR, RANK, TOP_N)
    generate_taxa_sanky(dfs['gtdb'], OUTPUT_DIR, RANK)
    taxa_sanky_rank(dfs['gtdb'], OUTPUT_DIR, RANK)

    print(f"[INFO] All test plots saved in '{OUTPUT_DIR}'")
    print(f"[INFO] Total runtime: {time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time))}")
