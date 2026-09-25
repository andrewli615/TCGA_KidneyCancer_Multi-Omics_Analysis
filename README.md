# TCGA Kidney Cancer Multi-Omics Analysis

A 2025 R project by Andrew Li using clinical, mutation, and RNA-seq data from TCGA kidney renal clear cell carcinoma. The matched cohort has **354 patients**; survival analyses use **352** with positive follow-up.

## Run

The compressed inputs are included. In RStudio, run:

```r
install.packages(c("rmarkdown", "dplyr", "ggplot2", "survival", "BiocManager"))
BiocManager::install("DESeq2")
rmarkdown::render("project.Rmd")
```

[project.Rmd](project.Rmd) runs Kaplan–Meier and Cox survival analysis, PCA, hierarchical clustering, VHL mutant versus wild-type expression tests, DESeq2 comparison of the two expression clusters, and GO biological-process enrichment. [loader.qmd](loader.qmd) checks the cohort intersection.

## Outputs

- [Mutation figure](top20_mutated_genes.png), [PCA](pca_top_variable_genes.png), [cluster heatmap](expression_cluster_heatmap.png), [DESeq2 volcano plot](volcano_DESeq2_expression_clusters.png), and [GO enrichment plot](go_bp_enrichment.png)
- [Differential-expression and enrichment tables](results/)
- [Input provenance](data/README.md) and [GO reference provenance](reference/README.md)

The GO results include immune-response, ion-transport, and ion-homeostasis terms. These are **exploratory findings**: the clusters were learned from the same expression data used for differential expression, and their biological interpretation needs independent validation.
