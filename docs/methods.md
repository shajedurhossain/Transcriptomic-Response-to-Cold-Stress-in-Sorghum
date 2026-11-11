# Methods

## Plant material and experimental design

Six *Sorghum bicolor* lines were selected for contrasting frost tolerance:

- **Low frost tolerance:** SB14011 (Sb1), GG12 338 (Sb4) and B18107 (Sb5).
- **High frost tolerance:** GG12 5051 (Sb2), ES-17-02 (Sb3) and GG12 477 (Sb6).

Two-week-old seedlings were exposed to **−5 °C for 2 h in the dark** (frost stress). Control plants were grown under the same conditions without the frost treatment.
Each line × condition combination has three biological replicates. **Each library is a pool of three plants.**

**Growth batch.** Three genotypes had to be regrown. Their libraries carry the numbers 4–6 (`SbX-C4/5/6`, `SbX-S4/5/6`), and pools 1–3 of those plants were not sequenced:

- Sb4 and Sb6: all libraries regrown.
- Sb3: all libraries regrown except Sb3-S1.

The sample sheet records the original library ID and the batch of each sample (`library_id`, `growth_batch`). The analysis names (`Sb3_HT_C_1` …) number the replicates 1–3 within each line and condition.

## Data

The data are gene-level read counts for the 36 libraries (108–172 million assigned reads each; 30,945 genes).
Counts were generated against the *S. bicolor* v5.1 annotation. For iDEP, the gene IDs were renamed to the Ensembl v3 pattern (`Sobic.001G000300.v5.1` → `SORBI_3001G000300`; `Sbicolor…` and `Sobic.K…` scaffold genes accordingly). The counts themselves are unchanged.

## Pre-processing

Genes were kept if they had ≥ 0.5 counts per million in at least one sample (22,810 genes), the same filter as iDEP.
DESeq2 median-of-ratios normalisation was used for differential expression.
The variance-stabilising transformation (`vst`, blind) was used for PCA.

## Differential expression

All models were fitted with DESeq2 1.46 (R 4.4) using the Wald test and independent filtering (α = 0.1).
Genes were called differentially expressed at **Benjamini–Hochberg FDR < 0.1 and |log2 fold change| > 1**, the thresholds used in the thesis.
All contrasts are oriented so that a positive log2FC means *higher under Cold* or *higher in High-tolerance lines*.

### Reproduction of the iDEP models (`scripts/01_reproduce_idep.R`)

| Model | Design | Contrasts |
|---|---|---|
| A | `~ tolerance + condition` | High vs Low tolerance; Cold vs Control |
| B | `~ tolerance + condition + tolerance:condition` | interaction `toleranceHigh.conditionCold` |

### Genotype-aware models (`scripts/02_deseq2_genotype_aware.R`)

**M1, lines nested in tolerance groups.** This follows the DESeq2 vignette section *"Group-specific condition effects, individuals nested within groups"*.
`genotype_n` numbers the lines 1–3 within each tolerance group:

```r
design = ~ tolerance + tolerance:genotype_n + tolerance:condition
```

- `toleranceLow.conditionCold`: cold response in low-tolerance lines
- `toleranceHigh.conditionCold`: cold response in high-tolerance lines
- `list("toleranceHigh.conditionCold", "toleranceLow.conditionCold")`: difference in cold response (interaction)
- 0.5 × each of the above two coefficients: average cold response over all lines

**M2, per-genotype cold response.** `~ genotype + genotype:condition`, with one Cold vs Control coefficient per line.

**M3, tolerance groups at the line level.** For each condition, the three replicates of each line were summed (pseudo-bulk; 6 samples), and `~ tolerance` was fitted.
This treats the line, not the library, as the unit of replication. That is the appropriate test of a tolerance-group difference when there are only three lines per group.

**Batch sensitivity (`scripts/04_batch_sensitivity.R`).** (a) M1 re-fitted without Sb3-S1, the only Sb3 library from the original batch. (b) M3 re-fitted as `~ batch + tolerance`, with Sb3 counted as regrown.

## Annotation and enrichment

Gene descriptions come from the gene-info tables exported by iDEP. Chromosome is parsed from the gene ID (`SORBI_3<chr>G…`).
GO Biological Process enrichment was done in iDEP for the pooled contrasts (`idep/results/*/GO_enrichment.csv`).

## Software

R 4.4.3, DESeq2 1.46.0, ggplot2 4.0.3, ggrepel 0.9.8, pheatmap 1.0.13 (full session info in `results/sessionInfo_*.txt`); iDEP v2.4.3 (DESeq2 1.48.2, R 4.5.1).
