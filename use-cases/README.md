# Example use case
All commands below assume that the command-line tool `mags-visualization`is installed and available in your environment.

Each plot is generated using a dedicated **subcommand**.

The working directory should be **MAGs-visualization/use-cases** 
so that the example paths resolve correctly.

If the command is not available, install the package first (['README.md'](../README.md))

# Bash
## Bee Gut Microbiome

```bash
use_case_folder="bee-use-case"
```

### Sample heatmap

* Genus

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --metadata "$use_case_folder/data/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --tax_level genus \
  --top_bar_spacer -1.0 \
  --spacer_meta 2.5 \
  --output "$use_case_folder/plots" \
  --no_log \
  --top_bar_height 2.0
```

* Specis

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --metadata "$use_case_folder/data/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --tax_level species \
  --top_bar_spacer -1.0 \
  --spacer_meta 2.5 \
  --output "$use_case_folder/plots" \
  --no_log \
  --top_bar_height 2.0
```

### dRep cluster plot

```bash
mags-visualization drep-cluster \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Assembly quality (QUAST)

```bash
mags-visualization assembly-quality \
  --quast "$use_case_folder/quast.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --color_by tax \
  --tax_level phylum \
  --output out
```

### Bakta annotation

```bash
mags-visualization bakta-annotation \
  --bakta "$use_case_folder/bakta.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```

## Macroalgal microbiome

```bash
use_case_folder="marine-use-case"
```

### Sample heatmap

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --tax_level family \
  --no_log \
  --top_bar_height 2.0 \
  --top_bar_spacer -1 \
  --output "$use_case_folder/plots"
```

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --tax_level species \
  --output "$use_case_folder/plots"
```

### dRep cluster plot

```bash
mags-visualization drep-cluster \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Assembly quality (QUAST)

```bash
mags-visualization assembly-quality \
  --quast "$use_case_folder/quast.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --color_by tax \
  --tax_level phylum \
  --output out
```

### Bakta annotation

```bash
mags-visualization bakta-annotation \
  --bakta "$use_case_folder/bakta.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```


## Cloud

```bash
use_case_folder="cloud-use-case/data"
```

### Sample heatmap

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" \
  --tax_level genus \
  --top_bar_spacer -0.5 \
  --no_log \
  --top_bar_height 3.5 \
  --output out
```

If without --no_log:
```bash
--top_bar_height 1.4
```

### dRep cluster plot

```bash
mags-visualization drep-cluster \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Assembly quality (QUAST)

```bash
mags-visualization assembly-quality \
  --quast "$use_case_folder/quast.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --color_by tax \
  --tax_level phylum \
  --output out
```

### Bakta annotation

```bash
mags-visualization bakta-annotation \
  --bakta "$use_case_folder/bakta.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```

## Termite Head Microbiome

```bash
use_case_folder="termite-use-case"
```

### Sample heatmap

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --metadata "$use_case_folder/data/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level genus \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --metadata "$use_case_folder/data/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level species \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/data/coverm.tsv" \
  --gtdb "$use_case_folder/data/gtdb.tsv" \
  --metadata "$use_case_folder/data/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level family \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

### dRep cluster plot

```bash
mags-visualization drep-cluster \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Assembly quality (QUAST)

```bash
mags-visualization assembly-quality \
  --quast "$use_case_folder/quast.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --color_by tax \
  --tax_level genus \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```


# Powershell
If you want to use powershell instead of bash, it's just a different syntax.
In the following examples the **subcommand** "all" is used, which produces all plots.

## Macroalgal microbiome

```powershell
$use_case_folder = "marine-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --color_by tax `
  --tax_level genus `
  --max_col 10 `
  --top_n 30 `
  -o test
```

## Bee

```powershell
$use_case_folder = "bee-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" `
  --color_by tax `
  --tax_level phylum `
  --top_n 30 `
  --top_bar_spacer -0.5 `
  --spacer_meta 2.5 `
  -o test
```

```powershell
--no_log `
--top_bar_height 2.0
```

## Cloud

```powershell
$use_case_folder = "cloud-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" `
  --color_by tax `
  --tax_level phylum `
  --top_n 30 `
  --top_bar_spacer -0.5 `
  --no_log `
  --top_bar_height 3.5 `
  -o test
```

If without --no_log:
```powershell
--top_bar_height 1.4
```

## Termite

```powershell
$use_case_folder = "termite-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Species" "Casts" "Colony" `
  --color_by tax `
  --tax_level genus `
  --max_col 10 `
  --top_n 30 `
  --spacer_meta 4.0 `
  --no_log `
  -o test
```
