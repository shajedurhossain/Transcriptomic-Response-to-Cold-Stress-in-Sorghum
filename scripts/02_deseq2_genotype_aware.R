#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 02_deseq2_genotype_aware.R
#
# Corrected differential-expression analysis that accounts for the structure of
# the experiment: 6 genotypes (3 low- and 3 high-cold-tolerance lines), each
# sampled 3x under Control and 3x under Cold.
#
# Why: the iDEP models (script 01) ignore genotype. Genotype is the dominant
# source of variance (PCA: PC1-PC3 correlate with genotype), so pooling lines
# (i) inflates tolerance-group differences with genotype-specific genes
#     (many |log2FC| > 20 "presence/absence" genes), and
# (ii) makes the "LT (C vs S)" and "HT (C vs S)" runs identical, because an
#     additive model has one cold effect for everybody.
#
# Models fitted here
#   M1  nested design (DESeq2 vignette, "individuals nested within groups"):
#       ~ tolerance + tolerance:genotype_n + tolerance:condition
#       -> cold response within Low lines, within High lines, their difference
#          (interaction), and the average cold response.
#   M2  per-genotype cold response: ~ genotype + genotype:condition
#   M3  constitutive tolerance-group difference at the genotype level:
#       replicates summed per genotype x condition (n = 3 vs 3 lines),
#       ~ tolerance, fitted separately for Control and Cold samples.
#
# Usage (from repo root):  Rscript scripts/02_deseq2_genotype_aware.R
# -----------------------------------------------------------------------------

source(file.path("scripts", "utils.R"))
set.seed(SEED)

dat <- load_counts()
counts <- filter_genes(dat$counts)
samples <- dat$samples
annot <- load_annotation()

summary_rows <- list()
add_summary <- function(df, name, description, model) {
  summary_rows[[name]] <<- count_degs(df, name, description, model)
  write_results(df, name)
}

# ---- M1: nested design ------------------------------------------------------
dds1 <- DESeqDataSetFromMatrix(counts, samples,
  design = ~ tolerance + tolerance:genotype_n + tolerance:condition
)
dds1 <- DESeq(dds1, quiet = TRUE)
rn <- resultsNames(dds1)
message("M1 coefficients: ", paste(rn, collapse = ", "))

cold_low <- tidy_results(results(dds1, name = "toleranceLow.conditionCold", alpha = FDR), annot)
cold_high <- tidy_results(results(dds1, name = "toleranceHigh.conditionCold", alpha = FDR), annot)
cold_diff <- tidy_results(results(dds1,
  contrast = list("toleranceHigh.conditionCold", "toleranceLow.conditionCold"), alpha = FDR
), annot)
w <- setNames(numeric(length(rn)), rn)
w[c("toleranceLow.conditionCold", "toleranceHigh.conditionCold")] <- 0.5
cold_avg <- tidy_results(results(dds1, contrast = w, alpha = FDR), annot)

add_summary(cold_avg, "cold_vs_control_all_lines", "Cold vs Control, average over all 6 lines", "M1 nested")
add_summary(cold_low, "cold_vs_control_low_tolerance", "Cold vs Control within low-tolerance lines (Sb1, Sb4, Sb5)", "M1 nested")
add_summary(cold_high, "cold_vs_control_high_tolerance", "Cold vs Control within high-tolerance lines (Sb2, Sb3, Sb6)", "M1 nested")
add_summary(cold_diff, "cold_response_High_vs_Low", "Interaction: cold response in High minus Low tolerance lines", "M1 nested")

# ---- M2: cold response per genotype -----------------------------------------
dds2 <- DESeqDataSetFromMatrix(counts, samples, design = ~ genotype + genotype:condition)
dds2 <- DESeq(dds2, quiet = TRUE)
per_geno_lfc <- list()
for (g in levels(samples$genotype)) {
  coef <- paste0("genotype", g, ".conditionCold")
  df <- tidy_results(results(dds2, name = coef, alpha = FDR), annot)
  tol <- as.character(unique(samples$tolerance[samples$genotype == g]))
  add_summary(df, paste0("cold_vs_control_", g),
              sprintf("Cold vs Control in %s (%s, %s tolerance)", g,
                      unique(samples$line_name[samples$genotype == g]), tolower(tol)),
              "M2 per genotype")
  per_geno_lfc[[g]] <- setNames(df$log2FoldChange, df$gene_id)
}

# ---- M3: tolerance group at the genotype level (pseudo-bulk) -----------------
pseudo <- function(cond) {
  keep <- samples$condition == cond
  sub <- samples[keep, ]
  cs <- sapply(split(sub$sample_id, droplevels(sub$genotype)), function(ids) rowSums(counts[, ids]))
  cd <- unique(sub[, c("genotype", "line_name", "tolerance")])
  rownames(cd) <- cd$genotype
  cd <- cd[colnames(cs), ]
  dds <- DESeqDataSetFromMatrix(cs, cd, design = ~ tolerance)
  DESeq(dds, quiet = TRUE)
}
for (cond in c("Control", "Cold")) {
  dds3 <- pseudo(cond)
  df <- tidy_results(results(dds3, contrast = c("tolerance", "High", "Low"), alpha = FDR), annot)
  add_summary(df, paste0("tolerance_High_vs_Low_", tolower(cond)),
              sprintf("High vs Low tolerance lines under %s (genotype-level, n = 3 vs 3)", cond),
              "M3 pseudo-bulk")
}

# ---- Summary ----------------------------------------------------------------
summary_tab <- do.call(rbind, summary_rows)
rownames(summary_tab) <- NULL
write.csv(summary_tab, file.path(TABLE_DIR, "summary_DEG_counts.csv"), row.names = FALSE)
print(summary_tab[, c("contrast", "up", "down", "total")])

# Overlap of the cold response between tolerance groups
sig_low <- cold_low$gene_id[cold_low$significant]
sig_high <- cold_high$gene_id[cold_high$significant]
overlap <- data.frame(
  set = c("low_only", "high_only", "shared", "shared_same_direction"),
  n = c(length(setdiff(sig_low, sig_high)), length(setdiff(sig_high, sig_low)),
        length(intersect(sig_low, sig_high)),
        sum(sapply(intersect(sig_low, sig_high), function(g)
          cold_low$direction[cold_low$gene_id == g] == cold_high$direction[cold_high$gene_id == g])))
)
write.csv(overlap, file.path(TABLE_DIR, "overlap_cold_low_vs_high.csv"), row.names = FALSE)
print(overlap)

# Per-genotype log2FC matrix (for figures) -- effect sizes only, no counts
lfc_mat <- do.call(cbind, lapply(per_geno_lfc, function(x) x[rownames(counts)]))
saveRDS(list(dds1 = dds1, lfc_mat = lfc_mat), file.path(OBJ_DIR, "genotype_aware_models.rds"))
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo_02.txt"))
