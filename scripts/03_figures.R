#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 03_figures.R -- publication-style figures from scripts 01/02.
# Usage (from repo root, after 01 and 02):  Rscript scripts/03_figures.R
# -----------------------------------------------------------------------------

source(file.path("scripts", "utils.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

theme_set(theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ecebe7", linewidth = 0.4),
        axis.text = element_text(colour = "#52514e"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "#52514e"),
        legend.position = "top", strip.text = element_text(face = "bold")))

save_fig <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(FIG_DIR, name), p, width = w, height = h, dpi = 200, bg = "white")
}
read_de <- function(name) read.csv(gzfile(file.path(TABLE_DIR, paste0("DE_", name, ".csv.gz"))))

m <- readRDS(file.path(OBJ_DIR, "genotype_aware_models.rds"))
dds1 <- m$dds1
samples <- as.data.frame(colData(dds1))

# ---- 1. PCA on VST data ------------------------------------------------------
vsd <- vst(dds1, blind = TRUE)
pc <- prcomp(t(assay(vsd)[head(order(-apply(assay(vsd), 1, var)), 500), ]))
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pcd <- cbind(samples, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
cent <- aggregate(cbind(PC1, PC2) ~ genotype + tolerance + condition, pcd, mean)
regrown <- c("Sb3", "Sb4", "Sb6") # second growth batch (library numbers 4-6)
cent$lab <- ifelse(cent$genotype %in% regrown, paste0(cent$genotype, "*"), as.character(cent$genotype))
p <- ggplot(pcd, aes(PC1, PC2)) +
  geom_point(aes(colour = tolerance, shape = condition), size = 3, stroke = 1.1) +
  geom_text_repel(data = cent, aes(label = lab), size = 3.2, colour = "#2b2b2a",
                  box.padding = 0.6, min.segment.length = 0.2, seed = SEED) +
  scale_colour_manual(values = COL_TOL, name = "Frost tolerance") +
  scale_shape_manual(values = c(Control = 16, Cold = 2), name = "Condition") +
  labs(title = "Frost treatment is the main axis of variation (PC1)",
       subtitle = "PCA of VST-normalised expression, 500 most variable genes",
       caption = "* regrown genotypes (second growth batch); PC2 separates them from the others",
       x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"))
save_fig(p, "fig1_pca.png", 8, 6)

# ---- 2. DEG counts per contrast ---------------------------------------------
s <- read.csv(file.path(TABLE_DIR, "summary_DEG_counts.csv"))
lab <- c(
  cold_vs_control_all_lines = "Cold vs Control, all lines",
  cold_vs_control_low_tolerance = "Cold vs Control, low-tolerance lines",
  cold_vs_control_high_tolerance = "Cold vs Control, high-tolerance lines",
  cold_response_High_vs_Low = "Interaction: cold response High minus Low",
  tolerance_High_vs_Low_control = "High vs Low lines, under Control",
  tolerance_High_vs_Low_cold = "High vs Low lines, under Cold"
)
sm <- s[s$contrast %in% names(lab), ]
long <- rbind(data.frame(contrast = sm$contrast, dir = "Up", n = sm$up),
              data.frame(contrast = sm$contrast, dir = "Down", n = -sm$down))
long$label <- factor(lab[long$contrast], levels = rev(lab))
p <- ggplot(long, aes(n, label, fill = dir)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = abs(n), hjust = ifelse(n > 0, -0.15, 1.15)), size = 3.4, colour = "#2b2b2a") +
  geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.4) +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = NULL,
                    labels = c(Up = "Up (higher in Cold / High-tolerance)", Down = "Down")) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.15)) +
  labs(title = "Differentially expressed genes (genotype-aware models)",
       subtitle = sprintf("DESeq2, FDR < %.2g and |fold change| > %g", FDR, FC), x = "Number of genes", y = NULL)
save_fig(p, "fig2_deg_counts.png", 9, 4.6)

# ---- 3. Cold-responsive genes per genotype ----------------------------------
pg <- s[grepl("^cold_vs_control_Sb", s$contrast), ]
pg$genotype <- sub("cold_vs_control_", "", pg$contrast)
pg <- merge(pg, unique(samples[, c("genotype", "line_name", "tolerance")]), by = "genotype")
pg$label <- paste0(pg$genotype, "\n", pg$line_name)
pg$label <- factor(pg$label, levels = pg$label[order(pg$tolerance, pg$genotype)])
pgl <- rbind(data.frame(pg[, c("label", "tolerance")], dir = "Up", n = pg$up),
             data.frame(pg[, c("label", "tolerance")], dir = "Down", n = pg$down))
p <- ggplot(pgl, aes(label, n, fill = dir)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.4, size = 3, colour = "#2b2b2a") +
  facet_grid(~ tolerance, scales = "free_x", space = "free_x",
             labeller = labeller(tolerance = c(Low = "Low frost tolerance", High = "High frost tolerance"))) +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = NULL,
                    labels = c(Up = "Cold-induced", Down = "Cold-repressed")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Cold-responsive genes in each genotype", x = NULL, y = "Number of DEGs",
       subtitle = sprintf("Cold vs Control within each line, FDR < %.2g, |FC| > %g", FDR, FC))
save_fig(p, "fig3_cold_degs_per_genotype.png", 9, 5)

# ---- 4. Volcano plots -------------------------------------------------------
vol <- function(name, title) {
  d <- read_de(name)
  d <- d[!is.na(d$padj), ]
  d$y <- pmin(-log10(d$padj), 100)
  d$dir <- factor(d$direction, levels = c("up", "down", "ns"))
  top <- head(d[d$significant & !is.na(d$description) & d$description != "Hypothetical protein", ], 8)
  ggplot(d, aes(log2FoldChange, y, colour = dir)) +
    geom_point(data = d[d$dir == "ns", ], size = 0.6, alpha = 0.5) +
    geom_point(data = d[d$dir != "ns", ], size = 0.9, alpha = 0.7) +
    geom_text_repel(data = top, aes(label = substr(description, 1, 32)), size = 2.6,
                    colour = "#2b2b2a", max.overlaps = 20, seed = SEED) +
    geom_vline(xintercept = c(-1, 1) * log2(FC), linetype = 2, colour = "#8a8984", linewidth = 0.3) +
    geom_hline(yintercept = -log10(FDR), linetype = 2, colour = "#8a8984", linewidth = 0.3) +
    scale_colour_manual(values = c(up = COL_UP, down = COL_DOWN, ns = COL_NS), guide = "none") +
    coord_cartesian(xlim = c(-12, 12)) +
    labs(title = title, x = expression(log[2] ~ "fold change"), y = expression(-log[10] ~ "FDR (capped at 100)"),
         subtitle = sprintf("%d up, %d down", sum(d$dir == "up"), sum(d$dir == "down")))
}
for (v in list(c("cold_vs_control_low_tolerance", "Cold vs Control, low-tolerance lines"),
               c("cold_vs_control_high_tolerance", "Cold vs Control, high-tolerance lines"),
               c("cold_response_High_vs_Low", "Interaction: cold response, High minus Low tolerance"))) {
  save_fig(vol(v[1], v[2]), paste0("fig4_volcano_", v[1], ".png"), 7, 5.5)
}

# ---- 5. Cold response: low vs high tolerance lines ---------------------------
lo <- read_de("cold_vs_control_low_tolerance")
hi <- read_de("cold_vs_control_high_tolerance")
it <- read_de("cold_response_High_vs_Low")
d <- merge(lo[, c("gene_id", "log2FoldChange", "significant")],
           hi[, c("gene_id", "log2FoldChange", "significant")], by = "gene_id", suffixes = c("_low", "_high"))
d <- merge(d, it[, c("gene_id", "significant")], by = "gene_id")
d$class <- with(d, ifelse(significant, "Different response (interaction FDR < 0.1)",
                   ifelse(significant_low & significant_high, "Cold-responsive in both",
                   ifelse(significant_low | significant_high, "Cold-responsive in one group", "Not cold-responsive"))))
d$class <- factor(d$class, levels = c("Not cold-responsive", "Cold-responsive in one group",
                                      "Cold-responsive in both", "Different response (interaction FDR < 0.1)"))
resp <- d$significant_low | d$significant_high
r <- cor(d$log2FoldChange_low[resp], d$log2FoldChange_high[resp], use = "complete.obs", method = "spearman")
both <- d$significant_low & d$significant_high
same <- mean(sign(d$log2FoldChange_low[both]) == sign(d$log2FoldChange_high[both]))
p <- ggplot(d[order(d$class), ], aes(log2FoldChange_low, log2FoldChange_high, colour = class)) +
  geom_abline(slope = 1, intercept = 0, colour = "#8a8984", linetype = 2, linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c("#d9d8d3", "#8fb8ea", COL_DOWN, COL_UP), name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  labs(title = "The cold response is largely shared between tolerance groups",
       subtitle = sprintf("Spearman rho = %.2f over cold-responsive genes\n%.1f%% of genes responsive in both groups change in the same direction",
                         r, 100 * same),
       x = expression("Cold " * log[2] * "FC, low-tolerance lines"),
       y = expression("Cold " * log[2] * "FC, high-tolerance lines"))
save_fig(p, "fig5_cold_lfc_low_vs_high.png", 7, 7)

# ---- 6. Heatmap: top interaction genes, per-genotype cold log2FC --------------
top <- head(it[it$significant & it$baseMean > 50, ], 40)
mat <- m$lfc_mat[top$gene_id, ]
ann <- unique(samples[, c("genotype", "tolerance")])
ann <- data.frame(Tolerance = ann$tolerance, row.names = ann$genotype)
ord <- rownames(ann)[order(ann$Tolerance)]
desc <- ifelse(is.na(top$description) | top$description == "Hypothetical protein", "",
               paste0("  ", substr(top$description, 1, 35)))
rownames(mat) <- paste0(top$gene_id, desc)
lim <- 6
png(file.path(FIG_DIR, "fig6_heatmap_interaction_genes.png"), width = 8, height = 10, units = "in", res = 200)
pheatmap(pmax(pmin(mat[, ord], lim), -lim),
         color = colorRampPalette(c(COL_DOWN, "#f2f1ee", COL_UP))(101),
         breaks = seq(-lim, lim, length.out = 102), cluster_cols = FALSE,
         annotation_col = ann, annotation_colors = list(Tolerance = COL_TOL),
         gaps_col = sum(ann$Tolerance == "Low"), border_color = NA, fontsize_row = 7,
         main = "Top 40 genes whose cold response differs between tolerance groups\n(cold vs control log2FC per genotype)")
invisible(dev.off())

message("Figures written to ", FIG_DIR)
