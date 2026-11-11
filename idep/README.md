# iDEP analysis (thesis version)

These are the outputs of the analysis done in the **iDEP v2.4.3** web application ([Ge et al. 2018](https://doi.org/10.1186/s12859-018-2486-6)), cleaned up for sharing.
Before using these files, read the notes on [direction](#direction-of-the-contrasts-read-this-first) and [duplicated runs](#five-runs-three-distinct-results).

## iDEP settings

| Step | Setting |
|---|---|
| Input | gene count matrix, 30,945 genes × 36 samples (v5.1 gene IDs renamed to `SORBI_3…` so iDEP can recognise them) |
| Filter | ≥ 0.5 CPM in ≥ 1 sample → **22,810 genes** (21,051 mapped to iDEP's Ensembl/STRING database) |
| Transformation (QC, PCA, clustering) | EdgeR log2-CPM, pseudocount 4; missing values imputed by gene median |
| Clustering | top 1,000 variable genes, average linkage, Pearson distance, genes centred, z-score cutoff 3 |
| Differential expression | DESeq2, Wald test, independent filtering, **FDR < 0.1, fold change > 2** |
| Enrichment | iDEP GO Biological Process, separately for up- and down-regulated genes |

## Direction of the contrasts (read this first)

iDEP wrote the contrasts as `LT-HT` and `C-S`, meaning **LT / HT** and **Control / Stress**.
The genes iDEP called "Upregulated" in the C vs S runs are therefore **higher in Control**, i.e. *repressed* by cold. This was verified against the count matrix: all 1,793 "Upregulated" genes have lower mean expression under cold.
The cleaned files are named by biological meaning, and each table keeps the original label in an `idep_label` column:

| Contrast | iDEP "Upregulated" means | iDEP "Downregulated" means |
|---|---|---|
| `tolerance_LT_vs_HT` | `higher_in_low_tolerance_lines` | `higher_in_high_tolerance_lines` |
| `cold_stress` (C vs S) | `cold_repressed` | `cold_induced` |
| `interaction` (`BHT.CS`) | `interaction_positive`: cold response stronger/more positive in high-tolerance lines | `interaction_negative` |

## Five runs, three distinct results

The thesis folder contained five DESeq2 runs. With the additive model `~ tolerance + condition`, changing the reference level does not change the estimated effect, so:

- **C (LT vs HT)** and **S (LT vs HT)** are the same pooled tolerance contrast (684/684 and 938/946 genes shared). The cleaned results use the *S (LT vs HT)* run.
- **LT (C vs S)** and **HT (C vs S)** are the same pooled cold contrast (2,042/2,042 and 1,786/1,793 genes shared). The cleaned results use the *LT (C vs S)* run.
- **Interaction** (`~ tolerance * condition`) is distinct.

The genotype-aware analysis in `scripts/02_deseq2_genotype_aware.R` estimates true per-group cold responses.

DEG counts reported by iDEP:

| Contrast | iDEP up | iDEP down | Meaning |
|---|---:|---:|---|
| LT-HT | 1,131 | 936 | 1,131 higher in low-tolerance lines, 936 higher in high-tolerance lines |
| C-S | 1,930 | 2,137 | 1,930 cold-repressed, 2,137 cold-induced |
| Interaction | 148 | 136 | |

(The gene tables contain slightly fewer rows, because iDEP exports gene info only for genes in its database.)

## Folder contents

```
exploratory/
  preprocessing/   library sizes, read distribution by chromosome, density/box plots, mean-SD plot
  pca/             PCA, 3D PCA, biplot, scree, t-SNE, MDS, top-loading genes for PC1/PC2
  clustering/      hierarchical and k-means heatmaps (top 1,000 genes), sample tree
results/<contrast>/
  genes_<meaning>.csv          DEG lists with iDEP gene info
  GO_enrichment.csv            enriched GO terms (FDR, fold enrichment, genes)
  GO_enrichment_<meaning>.png  enrichment bar plots
  enrichment_tree.png, enrichment_network.jpeg
  volcano.png, ma_plot.png, heatmap_degs.png, deg_counts.png, venn.png, upset.png
  chromosome_positions.csv     (tolerance only) DEG positions exported from iDEP
generated_code/    R scripts exported by iDEP, unmodified (see note below)
reports/           iDEP HTML reports (pre-processing, clustering, PCA, DEG)
```

**Note on `generated_code/`:** these scripts are kept exactly as iDEP exported them. They do not run as-is: `DESeqDataSetFromMatrix()` is never assigned to `dds`, and they expect iDEP's `converted_counts_data.csv`. Use `scripts/01_reproduce_idep.R` instead. It fits the same models and reproduces iDEP's DEG counts (within 7 genes).

PCA note: iDEP's PCA uses all filtered genes. PC1–PC3 correlate with genotype (p < 1e-7) and PC4 with condition (p = 7e-11). The PCA in `results/figures/fig1_pca.png` uses the 500 most variable VST genes, where the cold effect dominates PC1.
