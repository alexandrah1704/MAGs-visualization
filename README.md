<p align="left">
  <h1 align="left">Visualizations of MAGs</h1>
  <p align="left">A toolkit for visualizing MAG quality, taxonomy, clustering, assembly metrics, detection patterns and genome annotations.</p>
</p>

---

## Installation

### Prerequisites

- Python ≥ 3.11
- Conda (Miniconda, Miniforge, Mambaforge)
- Git

### Option 1: Conda / Mamba

Clone the repository and change into the project directory:

```bash
git clone https://github.com/alexandrah1704/MAGs-visualization.git
cd MAGs-visualization
```

Create conda environment and activate it:

```bash
conda env create -f environment.yml
conda activate mags
```

### Option 2: Python virtual environment (pip)

```powershell
# Change into project directory
cd MAGs-visualization

# Create virtual environment
python -m venv .venv

# Allow script execution for this session
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass

# Activate virtual environment
.\.venv\Scripts\Activate.ps1

# Install dependencies
pip install -r requirements.txt
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

All plots are saved in a user-defined output directory.

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
use_case_folder="../test-data/marine"

python ../scripts/main.py \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  -o new-test-plots
```

## Plot Configurations

### Taxonomic rank
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
--meta_cols weather temp ground # example columns
```

## Heatmap Options
### Plot features
```bash
--top_bar_height 0.8  # Height of top bar

--hspace 0.25 # Gap between top bar and heatmap

--heatmap_width 11.0

--spacer_legend 0.3 # Gap between legend and meta_bar

--spacer_meta 2.0 # Gap between meta_bar and heatmap

--spacer_heatmap # Gap between heatmap and histogram

--legend 2.5  # Size of legend

--meta_bar_add 1.5  # Additional width for meta_bar

--top_bar_spacer 0.0  # Gap between header and top bar

--max_col 10  # How many taxonomy names are shown (top 10)
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
Full examples can be found in ['use-cases/README.md'](use-cases/README.md)

### Example with color configuration
```bash
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


