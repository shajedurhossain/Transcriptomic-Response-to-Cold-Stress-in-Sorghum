# -----------------------------------------------------------------------------
# utils.R -- shared settings and helpers for the sorghum cold-stress RNA-seq
# analysis. Sourced by the numbered scripts; run those from the repo root.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- Paths ------------------------------------------------------------------
COUNTS_FILE <- file.path("data", "raw", "sorghum_gene_count_matrix.csv")
SAMPLES_FILE <- file.path("data", "metadata", "sample_sheet.csv")
TABLE_DIR <- file.path("results", "tables")
FIG_DIR <- file.path("results", "figures")
OBJ_DIR <- file.path("results", "rds") # intermediate objects (git-ignored)

for (d in c(TABLE_DIR, FIG_DIR, OBJ_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Analysis parameters (kept identical to the iDEP analysis in the thesis) --
MIN_CPM <- 0.5      # keep genes with >= 0.5 CPM ...
MIN_SAMPLES <- 1    # ... in at least this many samples
FDR <- 0.1          # adjusted p-value cutoff
FC <- 2             # fold-change cutoff (|log2FC| > 1)
SEED <- 42

# ---- Colours (validated categorical palette; colour-blind safe pairs) -------
COL_UP <- "#eb6834"    # higher in numerator (e.g. cold-induced)
COL_DOWN <- "#2a78d6"  # lower in numerator (e.g. cold-repressed)
COL_NS <- "#b9b8b2"
COL_TOL <- c(Low = "#eb6834", High = "#2a78d6")

# ---- Data loading -----------------------------------------------------------
load_counts <- function() {
  if (!file.exists(COUNTS_FILE)) {
    stop(
      "Count matrix not found at '", COUNTS_FILE, "'.\n",
      "The raw counts are not distributed with this repository; see data/README.md."
    )
  }
  counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
  samples <- read.csv(SAMPLES_FILE, stringsAsFactors = FALSE)
  stopifnot(setequal(colnames(counts), samples$sample_id))
  counts <- as.matrix(counts[, samples$sample_id])
  storage.mode(counts) <- "integer"

  samples$genotype <- factor(samples$genotype, levels = paste0("Sb", 1:6))
  samples$tolerance <- factor(samples$tolerance, levels = c("Low", "High"))
  samples$condition <- factor(samples$condition, levels = c("Control", "Cold"))
  samples$growth_batch <- factor(samples$growth_batch, levels = c("original", "regrown"))
  # genotype index *within* its tolerance group (1..3), for the nested design
  samples$genotype_n <- factor(ave(as.integer(samples$genotype), samples$tolerance,
    FUN = function(x) as.integer(factor(x))
  ))
  rownames(samples) <- samples$sample_id
  list(counts = counts, samples = samples)
}

# Same pre-filter as iDEP: >= MIN_CPM counts-per-million in >= MIN_SAMPLES samples
filter_genes <- function(counts, min_cpm = MIN_CPM, min_samples = MIN_SAMPLES) {
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  keep <- rowSums(cpm >= min_cpm) >= min_samples
  message(sprintf("Filter: %d of %d genes kept (>= %.1f CPM in >= %d sample[s])",
                  sum(keep), nrow(counts), min_cpm, min_samples))
  counts[keep, , drop = FALSE]
}

# ---- Annotation -------------------------------------------------------------
# Chromosome is encoded in the Sorghum bicolor v3 gene ID: SORBI_3<chr>G<num>
# (IDs of the form SORBI_3000xxxxx are unplaced / scaffold genes).
chromosome_from_id <- function(ids) {
  chr <- suppressWarnings(as.integer(sub("^SORBI_3(\\d{3})G.*$", "\\1", ids)))
  ifelse(is.na(chr) | chr == 0, "unplaced", as.character(chr))
}

# Optional gene descriptions (collected from the iDEP gene-info exports)
load_annotation <- function(file = file.path("data", "metadata", "gene_descriptions.csv")) {
  if (!file.exists(file)) return(NULL)
  read.csv(file, stringsAsFactors = FALSE)
}

# ---- Results helpers --------------------------------------------------------
tidy_results <- function(res, annot = NULL) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$chromosome <- chromosome_from_id(df$gene_id)
  df$significant <- !is.na(df$padj) & df$padj < FDR & abs(df$log2FoldChange) > log2(FC)
  df$direction <- ifelse(!df$significant, "ns", ifelse(df$log2FoldChange > 0, "up", "down"))
  if (!is.null(annot)) df <- merge(df, annot, by = "gene_id", all.x = TRUE, sort = FALSE)
  first <- c("gene_id", "chromosome", "baseMean", "log2FoldChange", "lfcSE", "stat",
             "pvalue", "padj", "significant", "direction")
  df <- df[, c(first, setdiff(colnames(df), first))]
  df[order(df$padj, -abs(df$log2FoldChange), na.last = TRUE), ]
}

count_degs <- function(df, contrast, description, model) {
  data.frame(
    contrast = contrast, description = description, model = model,
    up = sum(df$direction == "up"), down = sum(df$direction == "down"),
    total = sum(df$significant), tested = sum(!is.na(df$padj))
  )
}

write_results <- function(df, name) {
  out <- file.path(TABLE_DIR, paste0("DE_", name, ".csv.gz"))
  con <- gzfile(out, "w")
  write.csv(df, con, row.names = FALSE)
  close(con)
  sig <- df[df$significant, ]
  write.csv(sig, file.path(TABLE_DIR, paste0("DEG_", name, ".csv")), row.names = FALSE)
  invisible(out)
}
