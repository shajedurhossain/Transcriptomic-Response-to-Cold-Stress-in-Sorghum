# Data

## `metadata/sample_sheet.csv`

One row per RNA-seq library (36 rows).

| column | meaning |
|---|---|
| `sample_id` | column name in the count matrix, `<genotype>_<LT/HT>_<C/S>_<replicate>` |
| `library_id` | original library ID, `SbX-C/S<n>`; numbers 4–6 are libraries from regrown plants |
| `genotype` | line code, Sb1–Sb6 |
| `line_name` | breeding line / accession name |
| `tolerance`, `tolerance_code` | cold tolerance group: Low (LT) or High (HT) |
| `condition`, `condition_code` | Control (C) or Cold stress (S) |
| `replicate` | replicate 1–3 within line × condition |
| `growth_batch` | `original` or `regrown` (Sb4 and Sb6 all libraries, Sb3 all except Sb3-S1) |
| `plants_pooled` | plants pooled per library (3) |

## `metadata/gene_descriptions.csv`

Gene symbol, Entrez ID and description for 5,020 genes, collected from the gene-info tables that iDEP exported (only genes iDEP reported as differentially expressed). The scripts use it to annotate the result tables.

## `raw/sorghum_gene_count_matrix.csv` (not included)

The raw count matrix is **not distributed** with this repository. The scripts expect a CSV file:

- first column `gene_id` (v5.1 gene IDs renamed to the Ensembl v3 pattern for iDEP, e.g. `SORBI_3001G000200`);
- one integer column per sample, named exactly as `sample_id` in the sample sheet;
- 30,945 genes × 36 samples.

To re-run the analysis, place the file at `data/raw/sorghum_gene_count_matrix.csv`. The file is listed in `.gitignore`, so it will not be committed by accident. For access to the data, please contact the author.
