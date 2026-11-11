#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 01_reproduce_idep.R
#
# Re-runs the DESeq2 models that iDEP v2.4.3 fitted for the thesis, so the
# web-app results can be checked from the command line. Two fixes relative to
# the iDEP-generated scripts in idep/generated_code/:
#   * the DESeqDataSet is actually assigned to `dds` (the iDEP code omits it);
#   * contrasts are written as Cold vs Control and High vs Low tolerance, so a
#     positive log2FC means "higher under cold" / "higher in tolerant lines".
#     (iDEP's "C-S" contrast was Control/Cold, i.e. the opposite sign.)
#
# Models (as in iDEP; genotype is NOT in the model):
#   A  ~ tolerance + condition               -> tolerance and cold main effects
#   B  ~ tolerance + condition + tolerance:condition -> interaction
#
# Usage (from repo root):  Rscript scripts/01_reproduce_idep.R
# -----------------------------------------------------------------------------

source(file.path("scripts", "utils.R"))

dat <- load_counts()
counts <- filter_genes(dat$counts)
samples <- dat$samples
annot <- load_annotation()

# ---- Model A: additive ------------------------------------------------------
dds_a <- DESeqDataSetFromMatrix(counts, samples, design = ~ tolerance + condition)
dds_a <- DESeq(dds_a, quiet = TRUE)

res_tol <- tidy_results(results(dds_a, contrast = c("tolerance", "High", "Low"), alpha = FDR), annot)
res_cold <- tidy_results(results(dds_a, contrast = c("condition", "Cold", "Control"), alpha = FDR), annot)

# ---- Model B: with interaction -----------------------------------------------
dds_b <- DESeqDataSetFromMatrix(counts, samples, design = ~ tolerance + condition + tolerance:condition)
dds_b <- DESeq(dds_b, quiet = TRUE)
res_int <- tidy_results(results(dds_b, name = "toleranceHigh.conditionCold", alpha = FDR), annot)

write_results(res_tol, "idep_tolerance_High_vs_Low")
write_results(res_cold, "idep_cold_vs_control")
write_results(res_int, "idep_interaction")

summary_tab <- rbind(
  count_degs(res_tol, "idep_tolerance_High_vs_Low", "High vs Low tolerance lines (all samples pooled)", "~ tolerance + condition"),
  count_degs(res_cold, "idep_cold_vs_control", "Cold vs Control (all lines pooled)", "~ tolerance + condition"),
  count_degs(res_int, "idep_interaction", "Tolerance x cold interaction (no genotype term)", "~ tolerance * condition")
)

# Numbers reported by the iDEP web app (sig_gene_stats.png in idep/results/).
# iDEP signs were LT/HT and Control/Cold, so up/down are swapped here.
summary_tab$idep_up <- c(936, 2137, 148)
summary_tab$idep_down <- c(1131, 1930, 136)
write.csv(summary_tab, file.path(TABLE_DIR, "summary_idep_reproduction.csv"), row.names = FALSE)
print(summary_tab)

saveRDS(list(dds_a = dds_a, dds_b = dds_b), file.path(OBJ_DIR, "idep_models.rds"))
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo_01.txt"))
