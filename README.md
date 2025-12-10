<p align="left">
  <h1 align="left">Visualizations of MAGs</h1>
  <p align="left">A toolkit for visualizing MAG quality, taxonomy, clustering, assembly metrics, detection patterns and genome annotations.</p>
</p>

---

## Install

With conda/mamba:

```
conda env create -f environment.yml
conda activate mags
```

---

## What is this tool ?

This tool generates a variety of visualizations for MAGs, including:

- Taxonomic sankey diagrams
- Completeness-/Contamination-Plots
- Heatmaps
- dRep-cluster visualization
- assembly quality plot
- bakta plot
- Rank distribution diagram...

All plots are saved in an user-defined output directory.

---

## Input files

Below are the inputs for a complete visualization run:

| Argument           | Description                           |
|--------------------|---------------------------------------|
| --coverm           | CoverM table                          |
| --checkm           | CheckM result file                    |
| --checkm2          | CheckM2 result file                   |
| --gtdb             | GTDB annotation table                 |
| --drep             | dRep cluster table                    |
| -o                 | Output folder for all generated plots |

Optional:

| Argument               | Description                           |
|------------------------|---------------------------------------|
| --quast                | QUAST assembly statistics             |
| --bakta                | Bakta annotation table                |
| --metadata             | Metadata table for coloring plots     |
| --metadata_heatmap_new | Metadata for heatmap visualization    |
| --amber                | CAMI Amber binning evaluation         |

## Command-Line usage

### Show help
```bash
python main.py --help
```

### Basic example
```bash
python scripts/main.py \
--coverm test-data/coverm.tsv \
--checkm test-data/checkm.tsv \
--checkm2 test-data/checkm2.tsv \
--gtdb test-data/gtdb.tsv \
--drep test-data/drep.csv \
-o new-test-plots
```

## Plot Configurations

### Taxonomix rank
```bash
--rank phylum
```

Available ranks:
```pgsql
domain, phylum, class, order, family, genus, species
```

### Top N taxa for plots
```bash
--top_n_counts 10
```
Minimum and Default = 5

### Plot size
```bash
--fig_size WIDTH HEIGHT
```

### Output format
```bash
--format png    # png, pdf or svg
```

### Coloring mode
```bash
--quality   # color points by quality categories hq, mq, lq
or
--color_by quality

--tax       # color by taxonomy
--color_by tax
--tax_level genus

--color_by meta  # color by metadata
--meta_col temperature  # weather or others
--meta_bin_width 5  # for numeric columns
```
To show in the heatmap more than one metadata column:
```bash
--meta_cols weather temp ground #example columns
```

## Bakta Options
### Choose annotation features
```bash
--bakta_metrics cds hypotheticals rrnas trnas crispr
```

### Plot feature
```bash
--ratio  # e.g. hypotheticals/CDs
```
Example:
```bash
--bakta_metrics hypotheticals
--ratio
```

## Assembly Options
### Choose annotation features
```bash
--column_choice "N50" "GC (%)" ...
--color_by quality  # or tax, meta
```

## dRep Options
```bash
--top_n 30  # show top 30 clusters with most cluster members
```

## Examples
### Example without configurations
```bash
python scripts/main.py \
--coverm test-data/coverm.tsv \
--checkm test-data/checkm.tsv \
--checkm2 test-data/checkm2.tsv \
--gtdb test-data/gtdb.tsv \
--drep test-data/drep.csv \
--quast test-data/quast.tsv \
--bakta test-data/bakta.tsv \
--metadata test-data/metadata.tsv \
--metadata_heatmap_new test-data/metadata_heatmap_new.tsv \
-o out \
```
The plots you will get:
Completeness-Contamination plots colored by quality and taxonomy 'pyhlum' by default.
Sankey_plots, heatmap with abundance, rank distribution pie, some assemlby_quality plots and bakta plots.

### Example with color configuration
```bash
python scripts/main.py \
--coverm ... --checkm ... checkm2 ... etc. \
-o out \
--quality   # or color_by quality

--tax             # → Checkm/Checkm2 plots colored by taxonomy \
--color_by tax    # → Assembly- and -Bakta plots colored by taxonomy \
--tax_level genus \

--color_by meta \
--meta_col weather \
--bakta_metrics hypotheticals rrnas \
--ratio \
```
The plots you will get in addition:
Completeness-Contamination plots colored by meta weather, Bakta hypotheticals and rrnas with ratios and also assembly_quality colored by meta column weather.

If you want a specific assembly_quality plot from quast, you need to specify the column:
```bash
--column_choice "N50" "GC (%)" ... \
--color_by tax \
--tax_level genus \
```

### Example to color heatmap with more than one meta column
```bash
python scripts/main.py \
--coverm ... --checkm ... checkm2 ... etc. \
--metadata_heatmap_new test-data/metadata_heatmap_new.tsv \
-o out \
--meta_cols weather temp ground
```

### Example dRep to show most representative MAGs
```bash
python scripts/main.py \
--coverm test-data/coverm.tsv \
--checkm test-data/checkm.tsv \
--checkm2 test-data/checkm2.tsv \
--gtdb test-data/gtdb.tsv \
--drep test-data/drep.csv \
--quast test-data/quast.tsv \
--bakta test-data/bakta.tsv \
-o new-test-plots \
--top_n 30
```



