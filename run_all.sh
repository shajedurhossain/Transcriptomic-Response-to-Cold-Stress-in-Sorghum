#!/usr/bin/env bash
# Run the full analysis from the repository root.
set -euo pipefail
cd "$(dirname "$0")"
Rscript scripts/01_reproduce_idep.R
Rscript scripts/02_deseq2_genotype_aware.R
Rscript scripts/03_figures.R
Rscript scripts/04_batch_sensitivity.R
echo "Done. Tables in results/tables, figures in results/figures."
