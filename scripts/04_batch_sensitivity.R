#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 04_batch_sensitivity.R
#
# Three genotypes were regrown, and their libraries come from a second growth
# batch (library numbers 4-6): Sb4 and Sb6 completely, and Sb3 except Sb3-S1.
# The batch is nested in genotype, so the genotype terms in the main models
# already absorb it. This script checks what is left:
#   (a) Sb3 is the only genotype with mixed batches (Sb3-S1 original, all other
#       Sb3 libraries regrown). Re-fit M1 without Sb3-S1 and compare.
#   (b) The tolerance groups are unbalanced for batch (Low: 1 of 3 lines
#       regrown; High: 2 of 3). Re-fit M3 (pseudo-bulk) with batch as a
#       covariate. Sb3 is counted as regrown, the batch of most of its libraries.
#
# Usage (from repo root, after 02):  Rscript scripts/04_batch_sensitivity.R
# -----------------------------------------------------------------------------

source(file.path("scripts", "utils.R"))

dat <- load_counts()
counts <- filter_genes(dat$counts)
samples <- dat$samples
samples$growth_batch <- factor(samples$growth_batch, levels = c("original", "regrown"))

sig_ids <- function(df) df$gene_id[df$significant]
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
read_main <- function(name) read.csv(gzfile(file.path(TABLE_DIR, paste0("DE_", name, ".csv.gz"))))
rows <- list()

# ---- (a) M1 without the cross-batch library Sb3-S1 -------------------------
keep <- samples$library_id != "Sb3-S1"
dds <- DESeqDataSetFromMatrix(counts[, keep], samples[keep, ],
  design = ~ tolerance + tolerance:genotype_n + tolerance:condition)
dds <- DESeq(dds, quiet = TRUE)
alt <- list(
  cold_vs_control_low_tolerance = results(dds, name = "toleranceLow.conditionCold", alpha = FDR),
  cold_vs_control_high_tolerance = results(dds, name = "toleranceHigh.conditionCold", alpha = FDR),
  cold_response_High_vs_Low = results(dds, contrast = list("toleranceHigh.conditionCold",
                                                           "toleranceLow.conditionCold"), alpha = FDR)
)
for (nm in names(alt)) {
  a <- tidy_results(alt[[nm]]); m <- read_main(nm)
  rows[[nm]] <- data.frame(check = "M1 without Sb3-S1", contrast = nm,
    main_DEGs = sum(m$significant), sensitivity_DEGs = sum(a$significant),
    shared = length(intersect(sig_ids(a), sig_ids(m))), jaccard = round(jaccard(sig_ids(a), sig_ids(m)), 3))
}

# ---- (b) M3 pseudo-bulk with batch covariate --------------------------------
geno_batch <- c(Sb1 = "original", Sb2 = "original", Sb3 = "regrown",
                Sb4 = "regrown", Sb5 = "original", Sb6 = "regrown")
for (cond in c("Control", "Cold")) {
  sub <- samples[samples$condition == cond, ]
  cs <- sapply(split(sub$sample_id, droplevels(sub$genotype)), function(ids) rowSums(counts[, ids]))
  cd <- unique(sub[, c("genotype", "tolerance")]); rownames(cd) <- cd$genotype; cd <- cd[colnames(cs), ]
  cd$batch <- factor(geno_batch[as.character(cd$genotype)], levels = c("original", "regrown"))
  d3 <- DESeq(DESeqDataSetFromMatrix(cs, cd, design = ~ batch + tolerance), quiet = TRUE)
  a <- tidy_results(results(d3, contrast = c("tolerance", "High", "Low"), alpha = FDR))
  nm <- paste0("tolerance_High_vs_Low_", tolower(cond)); m <- read_main(nm)
  rows[[nm]] <- data.frame(check = "M3 + batch covariate", contrast = nm,
    main_DEGs = sum(m$significant), sensitivity_DEGs = sum(a$significant),
    shared = length(intersect(sig_ids(a), sig_ids(m))), jaccard = round(jaccard(sig_ids(a), sig_ids(m)), 3))
}

out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(TABLE_DIR, "batch_sensitivity.csv"), row.names = FALSE)
print(out)
