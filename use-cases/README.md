# Example use case

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
  -o test2
```