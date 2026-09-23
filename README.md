# Single-cell RNA-seq analysis of PBMC 3k

## Overview

This project presents a complete single-cell RNA sequencing (scRNA-seq) analysis workflow applied to the **PBMC 3k** dataset using **R** and **Seurat**.

The objective was to build a clear, reproducible and biologically interpretable workflow covering the main stages of a standard scRNA-seq analysis, from quality control to cell-type annotation and downstream biological exploration.

The analysis includes:

- quality control of single-cell transcriptomic data;
- detection and removal of predicted doublets;
- evaluation of cell-cycle effects;
- normalization using **SCTransform**;
- dimensionality reduction using **PCA** and **UMAP**;
- graph-based clustering;
- identification of cluster-specific marker genes;
- manual cell-type annotation based on canonical markers;
- sub-clustering of the T/NK compartment;
- pseudotime analysis of the monocyte compartment using **Slingshot**;
- exploratory analysis of candidate ligand–receptor expression;
- gene-signature scoring using **AUCell**.

## Dataset

The project uses the **PBMC 3k** dataset distributed through the SeuratData ecosystem. The starting dataset contains approximately **2,700 cells** and **13,714 detected genes**.

PBMCs (Peripheral Blood Mononuclear Cells) include several major immune-cell populations such as T lymphocytes, B lymphocytes, NK cells, monocytes and dendritic cells.

Useful resources:

- Seurat PBMC 3k tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial
- SeuratData: https://github.com/satijalab/seurat-data

## Analysis workflow

### 1. Quality control

Several commonly used scRNA-seq quality metrics were examined: the number of genes detected per cell (`nFeature_RNA`), total RNA counts (`nCount_RNA`), mitochondrial transcript percentage and ribosomal transcript percentage.

Cells with very low transcriptomic complexity, unusually high gene counts or excessive mitochondrial RNA were filtered before downstream analysis. Relationships between RNA counts, detected genes and mitochondrial content were also inspected to identify atypical cells.

### 2. Doublet detection

Potential doublets were detected using **scDblFinder**. Doublets correspond to droplets containing RNA from more than one cell and can create artificial hybrid transcriptomic profiles that distort dimensionality reduction, clustering and marker-gene identification.

Predicted doublets were removed before the main normalization and clustering workflow.

Reference: https://bioconductor.org/packages/scDblFinder/

### 3. Cell-cycle assessment

Cell-cycle scores were calculated using canonical S-phase and G2/M-phase gene sets. A PCA restricted to cell-cycle-associated genes was used as a diagnostic analysis to assess whether cell-cycle variation strongly structured the dataset.

This step was used for interpretation rather than as an automatic reason to regress out biological variation.

### 4. SCTransform normalization

The filtered singlet cells were normalized using **SCTransform v2**. SCTransform models UMI counts using a regularized negative-binomial framework and reduces the dependence between sequencing depth and normalized gene expression.

Mitochondrial percentage was included as a variable to regress during this step.

Reference: https://satijalab.org/seurat/articles/sctransform_vignette

### 5. Principal component analysis

PCA was performed on the SCTransform-normalized data. An elbow plot was examined to evaluate the amount of information captured by successive principal components.

The first **20 principal components** were retained for construction of the neighborhood graph, clustering and UMAP representation.

### 6. Neighborhood graph and clustering

A nearest-neighbor graph was constructed from the selected principal components. Graph-based clustering was then performed using the Louvain algorithm through Seurat.

The resulting clusters represent groups of cells with similar transcriptomic profiles. Cluster numbers themselves have no intrinsic biological meaning and were interpreted only after marker-gene analysis.

### 7. UMAP visualization

UMAP was used to project the high-dimensional transcriptomic structure into two dimensions.

The UMAP representation provides an intuitive visualization of transcriptionally similar cells but should not be interpreted as a quantitative measure of biological distance between clusters.

### 8. Marker-gene identification

Cluster-specific marker genes were identified using Seurat's differential-expression framework. Only positively enriched markers were retained, with minimum expression-frequency and log-fold-change thresholds.

Marker genes were compared with known immune-cell markers to support biological interpretation of each cluster.

Examples of markers examined include:

- T cells: `CD3D`, `CD3E`, `CCR7`, `IL7R`, `CD8A`;
- B cells: `MS4A1`, `CD79A`, `CD79B`;
- NK cells: `NKG7`, `GNLY`, `PRF1`;
- CD14+ monocytes: `CD14`, `LYZ`, `S100A8`, `S100A9`;
- FCGR3A+ monocytes: `FCGR3A`, `MS4A7`, `LST1`;
- dendritic cells: `FCER1A`, `CST3`;
- platelets: `PPBP`, `PF4`.

### 9. Cell-type annotation

Cell types were assigned manually after examination of cluster-specific differential-expression results, canonical immune-cell markers and global expression patterns shown in DotPlots and UMAP representations.

Manual annotation was deliberately preferred over assigning identities solely from cluster numbers.

### 10. T/NK sub-clustering

The T/NK compartment was isolated and reprocessed independently. A new SCTransform normalization, PCA, neighborhood graph, clustering and UMAP were performed on this subset.

The goal was to investigate finer transcriptional heterogeneity that may be obscured during clustering of all PBMC populations together.

Markers including `CCR7`, `IL7R`, `CD8A`, `GZMK`, `CCL5`, `NKG7`, `GNLY`, `PRF1`, `IL2RA`, `CTLA4` and `FOXP3` were examined.

This analysis is exploratory and fine cell-state annotation requires validation using multiple consistent markers.

### 11. Monocyte pseudotime analysis

A targeted trajectory analysis was performed on the monocyte compartment using **Slingshot**.

Rather than constructing a trajectory across unrelated PBMC lineages, the analysis was restricted to transcriptionally related monocyte populations. A pseudotime value was inferred for each cell to represent its relative position along the inferred transcriptional continuum.

Pseudotime should not be interpreted as real chronological time or as direct experimental proof of lineage progression.

Reference: https://bioconductor.org/packages/slingshot/

### 12. Candidate ligand–receptor exploration

Expression of selected ligand and receptor genes was compared across annotated cell populations.

This analysis provides a descriptive view of potentially compatible ligand–receptor expression patterns. Importantly, co-expression of a ligand and its receptor does **not** demonstrate actual cell-cell communication. A dedicated interaction-inference framework would be required for that conclusion.

### 13. AUCell gene-signature analysis

**AUCell** was used to estimate the enrichment of predefined gene signatures at the individual-cell level.

The analysis included illustrative myeloid and B-cell signatures and provided an independent way to examine whether expected transcriptional programs were enriched in their corresponding cell populations.

Reference: https://bioconductor.org/packages/AUCell/

## Main output files

The analysis generates several figures and tables, including:

- QC violin plots;
- RNA-count / gene-count scatterplots;
- predicted-doublet visualizations;
- cell-cycle PCA;
- PCA elbow plot;
- global UMAP clustering;
- canonical-marker DotPlot;
- annotated UMAP;
- T/NK sub-clustering UMAP;
- T/NK marker DotPlot;
- monocyte pseudotime representation;
- ligand–receptor candidate-expression DotPlot;
- AUCell signature maps;
- complete marker-gene tables;
- top marker genes per cluster;
- cell-type abundance summary;
- serialized Seurat object.

## Project structure

```text
project/
│
├── 01_projet.R
├── README.md
│
├── 03_results/
│   ├── figures
│   ├── marker tables
│   ├── cell-type summary
│   └── Seurat object
│
└── report/
    └── scRNA-seq analysis report
```

## Software and main R packages

The workflow was developed in R using mainly:

- Seurat
- SeuratData
- scDblFinder
- SingleCellExperiment
- Slingshot
- AUCell
- ggplot2
- dplyr
- patchwork

A random seed was set during the analysis to improve reproducibility of stochastic steps.

## Interpretation and limitations

This project is designed as a complete analytical workflow and educational study of PBMC single-cell transcriptomics.

Several limitations should be considered:

- the PBMC 3k dataset represents a small reference dataset rather than a biological cohort;
- cell-type annotation is based primarily on transcriptomic marker expression;
- UMAP is a visualization method and does not provide a direct quantitative biological distance;
- pseudotime represents an inferred transcriptional ordering and not real chronological time;
- candidate ligand–receptor expression does not establish functional cell-cell communication;
- fine immune-cell subtypes should ideally be validated using additional markers, reference mapping or complementary experimental information.

## Reproducibility

The project was organized so that the main analysis can be reproduced from the original PBMC 3k dataset.

For full reproducibility, the following should be retained with the analysis:

- R version;
- package versions;
- random seed;
- QC thresholds;
- number of principal components;
- clustering resolution;
- R session information.

## Use of artificial intelligence

Artificial intelligence was used as a **support tool** during the development of this project.

More specifically, **ChatGPT (OpenAI)** was used to assist with:

- reviewing and debugging parts of the R code;
- improving code clarity and organization;
- correcting and optimizing sections of the analysis workflow;
- identifying potential inconsistencies or fragile implementation choices;
- improving the structure, wording and scientific clarity of the written report;
- helping ensure that biological interpretations were appropriately qualified and did not exceed what the analyses could support.

The analyses were executed and inspected by the author, and the final methodological choices, biological interpretation and validation of the project remained under the responsibility of the author.

AI assistance was therefore used for **code correction, optimization and editorial support**, and not as a substitute for scientific validation, critical interpretation or authorship responsibility.

## Author

**Nahi El Akoum**  
Master 2 — Artificial Intelligence and Data Analysis in Biology (AIDA)  
Sorbonne Université

## Selected references

1. Stuart T. et al. *Comprehensive Integration of Single-Cell Data*. Cell, 2019. https://doi.org/10.1016/j.cell.2019.05.031
2. Hao Y. et al. *Integrated analysis of multimodal single-cell data*. Cell, 2021. https://doi.org/10.1016/j.cell.2021.04.048
3. Hafemeister C., Satija R. *Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression*. Genome Biology, 2019. https://doi.org/10.1186/s13059-019-1874-1
4. Germain P.-L. et al. *Doublet identification in single-cell sequencing data using scDblFinder*. F1000Research, 2021. https://doi.org/10.12688/f1000research.73600.1
5. Street K. et al. *Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics*. BMC Genomics, 2018. https://doi.org/10.1186/s12864-018-4772-0
6. Aibar S. et al. *SCENIC: single-cell regulatory network inference and clustering*. Nature Methods, 2017. https://doi.org/10.1038/nmeth.4463

## License

This repository is intended for educational and academic use.

If the project or part of the workflow is reused, please cite the original software packages and publications on which the analysis is based.
