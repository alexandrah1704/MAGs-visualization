<p align="left">
  <h1 align="left">Visualizations of MAGs</h1>
  <p align="left">A toolkit for visualizing MAG quality, taxonomy, clustering, assembly metrics, detection patterns and genome annotations.</p>
</p>

---

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

| Argument           | Description                           |
|--------------------|---------------------------------------|
| --quast            | QUAST assembly statistics             |
| --bakta            | Bakta annotation table                |
| --metadata         | Metadata table for coloring plots     |
| --metadata_heatmap | Metadata for heatmap visualization    |
| --amber            | CAMI Amber binning evaluation         |

## Command-Line usage

### Show help
```bash
python main.py --help
```

### Basic example
```bash
python main.py \
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

## Full Example
```bash
python main.py \
--coverm test-data/coverm.tsv \
--checkm test-data/checkm.tsv \
--checkm2 test-data/checkm2.tsv \
--gtdb test-data/gtdb.tsv \
--drep test-data/drep.csv \
--quast test-data/quast.tsv \
--bakta test-data/bakta.tsv \
--metadata test-data/metadata.tsv \
--metadata_heatmap test-data/metadata_heatmap.tsv \
-o new-test-plots \
--column_choice "N50" \
--color_by meta \
--meta_col weather \
--bakta_metrics hypotheticals \
--ratio
```









  
