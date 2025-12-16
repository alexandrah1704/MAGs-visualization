# Example use case
All commands below assume that you are working from the directory **MAGs-visualization/use-cases**.

# Bash
## Marine

```bash
use_case_folder="../test-data/marine"

python ../scripts/main.py \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --color_by tax \
  --tax_level genus \
  --max_col 10 \
  --top_n 30 \
  -o test
```

## Bee

```bash
use_case_folder="../test-data/bee"

python ../scripts/main.py \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --color_by tax \
  --tax_level phylum \
  --top_n 30 \
  --top_bar_spacer -0.5 \
  --spacer_meta 2.5 \
  -o test
```

## Cloud

```bash
use_case_folder="../test-data/cloud"

python ../scripts/main.py \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" \
  --color_by tax \
  --tax_level phylum \
  --top_n 30 \
  --top_bar_spacer -0.5 \
  --no_log \
  --top_bar_height 3.5 \
  -o test
```

## Termite

```bash
use_case_folder="../test-data/termite"

python ../scripts/main.py \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --color_by tax \
  --tax_level phylum \
  --top_n 30 \
  --spacer_meta 4.0 \
  --no_log \
  -o test
```

# Powershell
If you want to use powershell instead of bash, here are the examples.

## Marine

```powershell
$use_case_folder = "..\test-data\marine"

python ..\scripts\main.py `
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
$use_case_folder = "..\test-data\bee"

python ..\scripts\main.py `
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
$use_case_folder = "..\test-data\cloud"

python ..\scripts\main.py `
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
$use_case_folder = "..\test-data\termite"

python ..\scripts\main.py `
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