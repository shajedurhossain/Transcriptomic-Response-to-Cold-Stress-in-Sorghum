# Research plan

## Project brief

> Investigate transcriptomic responses to cold (frost) stress in sorghum. Based on the information given and the analyses possible with the iDEP web server, come up with a research plan.

| Genotype | Library | Frost tolerance |
|---|---|---|
| SB14011 | Sb1 | Low |
| GG12 5051 | Sb2 | High |
| ES-17-02 | Sb3 | High |
| GG12 338 | Sb4 | Low |
| B18107 | Sb5 | Low |
| GG12 477 | Sb6 | High |

- **Stress:** 2 h at −5 °C in the dark, seedlings at two weeks of age.
- **Control:** grown under the same conditions, without the stress treatment.
- **Libraries:** `SbX-C1/2/3` are control replicates and `SbX-S1/2/3` are stressed replicates. `SbX-C/S-4/5/6` are also replicates, from genotypes that had to be regrown (their pools 1–3 were not sequenced).
- **Pools:** each library represents 3 pooled plants.

## Objectives

1. Identify genes whose expression changes after a short frost exposure in sorghum seedlings (the general frost response).
2. Compare that response between low- and high-tolerance lines. Which genes respond only, or more strongly, in one group?
3. Look for baseline (constitutive) expression differences between tolerance groups that could contribute to tolerance.
4. Interpret the gene sets functionally (GO enrichment, pathways) and propose candidate genes for follow-up.

## Hypotheses

- **H1:** Frost triggers a large, shared transcriptional response (signalling, membrane and cell-wall remodelling, protective proteins) in all lines.
- **H2:** Tolerant lines differ from sensitive lines in part of that response. They may induce protective genes more strongly or limit damage-related responses.
- **H3:** Tolerant lines differ constitutively in genes that prime them for frost.

## Analysis plan (iDEP, with checks in R)

| Step | iDEP module | Purpose | Status / output |
|---|---|---|---|
| 1. Data preparation | Load data | Count matrix plus design file (factors: tolerance, condition) | ✅ `data/metadata/sample_sheet.csv` |
| 2. Pre-processing and QC | Pre-process | Filter low counts (≥ 0.5 CPM in ≥ 1 sample), check library sizes, transform (log-CPM) | ✅ `idep/exploratory/preprocessing/` |
| 3. Exploratory analysis | Clustering, PCA | Do samples group by condition, genotype or tolerance? Spot outliers | ✅ `idep/exploratory/`, `results/figures/fig1_pca.png` |
| 4. Differential expression | DEG1 / DEG2 (DESeq2) | Frost vs control; low vs high tolerance; tolerance × frost interaction | ✅ `idep/results/`; corrected in `results/tables/` |
| 5. Functional enrichment | Pathway / enrichment | GO Biological Process for up- and down-regulated genes in each contrast | ✅ `idep/results/*/GO_enrichment.csv` (pooled contrasts) |
| 6. Candidate genes | Gene lists, heatmaps | Rank frost-responsive and interaction genes; check known cold genes | 🔶 started, limited by annotation (see README) |
| 7. Reproducibility and robustness | (R, outside iDEP) | Re-fit the models in DESeq2 with genotype and batch accounted for | ✅ `scripts/01`–`04` |

### Design decisions

- **Genotype must be in the model.** Each tolerance group has only three lines, so samples are not independent within a group. Per-group frost responses are estimated with lines nested in tolerance groups. Group-level comparisons use the line as the unit (see `docs/methods.md`).
- **Regrown genotypes (libraries 4–6)** form a second growth batch. The batch is nested in genotype for Sb4 and Sb6, so the genotype terms absorb it. Sb3 is mixed (Sb3-S1 original), and the tolerance groups are unbalanced for batch. Both are handled by sensitivity checks.
- **Thresholds:** FDR < 0.1 and |fold change| > 2, as in the iDEP analysis.

### Expected outputs

- Lists of frost-induced and frost-repressed genes, overall and per tolerance group.
- Genes with a tolerance-dependent frost response (the interaction) as candidates for qPCR validation.
- Enriched biological processes for each gene set.
- A reproducible, documented analysis (this repository).
