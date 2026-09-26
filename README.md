# Single-cell RNA-seq analysis of PBMC 3k

[![R](https://img.shields.io/badge/R-single--cell%20analysis-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-scRNA--seq-6A5ACD)](https://satijalab.org/seurat/)
[![Dataset](https://img.shields.io/badge/Dataset-PBMC%203k-2E8B57)](https://satijalab.org/seurat/articles/pbmc3k_tutorial)

A portfolio project demonstrating a complete **single-cell RNA-seq analysis workflow** on the PBMC 3k reference dataset using **R, Seurat and Bioconductor**.

The objective is not only to recover the major PBMC populations, but to show the reasoning behind a standard scRNA-seq workflow: quality control, doublet detection, normalization, dimensionality reduction, clustering, marker analysis, biological annotation and targeted downstream exploration.

## Project at a glance

**Biological system:** human peripheral blood mononuclear cells  
**Data type:** single-cell RNA sequencing  
**Main language:** R  
**Core tools:** Seurat, SCTransform, scDblFinder, SingleCellExperiment, Slingshot, AUCell  
**Main tasks:** QC, clustering, annotation, sub-clustering, trajectory analysis and gene-signature scoring

## Biological question

Can major immune-cell populations and finer transcriptional structure be recovered from a reference PBMC scRNA-seq dataset using a rigorous and biologically interpretable analysis workflow?

The analysis focuses on both the **global immune-cell landscape** and selected downstream questions within transcriptionally related compartments.

## Dataset

The project uses the **PBMC 3k** dataset distributed through the Seurat ecosystem. It contains approximately **2,700 cells** and **13,714 detected genes**.

PBMCs include several major immune populations, including:

- CD4+ and CD8+ T cells
- B cells
- NK cells
- CD14+ monocytes
- FCGR3A+ monocytes
- dendritic cells
- platelets

Reference dataset and tutorial:  
https://satijalab.org/seurat/articles/pbmc3k_tutorial

## Analysis workflow

### 1. Quality control

Cell-level quality metrics were examined before downstream analysis:

- number of detected genes per cell
- total RNA counts
- mitochondrial transcript percentage
- ribosomal transcript percentage

Low-complexity cells, cells with unusually high feature counts and cells with excessive mitochondrial RNA were excluded using predefined QC criteria.

### 2. Doublet detection

Potential doublets were identified with **scDblFinder** and removed before the final normalization and clustering workflow.

Doublet removal is important because multiplets can generate artificial hybrid transcriptional profiles and distort both clustering and marker-gene identification.

### 3. Cell-cycle assessment

S-phase and G2/M scores were calculated using canonical cell-cycle gene sets.

A PCA restricted to cell-cycle-associated genes was used as a diagnostic step to assess whether cell-cycle variation strongly structured the dataset. Cell-cycle signal was interpreted rather than automatically removed.

### 4. SCTransform normalization

Filtered singlets were normalized with **SCTransform v2**, with mitochondrial percentage included as a technical covariate.

SCTransform was used to stabilize variance and reduce the dependence between sequencing depth and normalized expression.

### 5. PCA and dimensionality selection

Principal component analysis was performed on the normalized expression matrix.

An elbow plot was used to inspect the information carried by successive components, and the first **20 principal components** were retained for the main neighborhood graph, clustering and UMAP representation.

### 6. Graph-based clustering and UMAP

A nearest-neighbor graph was built from the selected PCs and clustered using Seurat's graph-based framework.

UMAP was used for two-dimensional visualization of the transcriptional structure.

Cluster labels were treated as computational groups only; biological meaning was assigned after marker analysis.

### 7. Marker-gene analysis and cell-type annotation

Cluster-specific markers were identified and compared with canonical immune-cell genes.

Representative markers included:

| Cell population | Example markers |
|---|---|
| T cells | CD3D, CD3E, CCR7, IL7R, CD8A |
| B cells | MS4A1, CD79A, CD79B |
| NK cells | NKG7, GNLY, PRF1 |
| CD14+ monocytes | CD14, LYZ, S100A8, S100A9 |
| FCGR3A+ monocytes | FCGR3A, MS4A7, LST1 |
| Dendritic cells | FCER1A, CST3 |
| Platelets | PPBP, PF4 |

Cell identities were assigned manually using multiple consistent markers rather than relying on cluster numbers alone.

### 8. T/NK sub-clustering

The T/NK compartment was isolated and reprocessed independently to explore finer immune heterogeneity that can be masked in the global PBMC analysis.

Markers examined included CCR7, IL7R, CD8A, GZMK, CCL5, NKG7, GNLY, PRF1, IL2RA, CTLA4 and FOXP3.

This analysis is exploratory; fine cell-state annotation requires support from multiple markers and, ideally, external validation.

### 9. Monocyte trajectory analysis

A targeted pseudotime analysis was performed with **Slingshot** on transcriptionally related myeloid populations.

The trajectory was restricted to the monocyte compartment rather than being forced across unrelated PBMC lineages.

Pseudotime is interpreted as an inferred transcriptional ordering, **not** as real chronological time or direct proof of lineage progression.

### 10. Candidate ligand-receptor exploration

Selected ligand and receptor genes were compared across annotated cell populations.

This analysis is descriptive: compatible ligand-receptor expression does **not** by itself demonstrate functional cell-cell communication.

### 11. AUCell gene-signature scoring

**AUCell** was used to score predefined transcriptional signatures at the single-cell level.

This provides an additional view of whether expected biological programs are enriched in the corresponding cell populations.

## Main skills demonstrated

- single-cell RNA-seq quality control
- Seurat-based preprocessing and clustering
- SCTransform normalization
- PCA and UMAP
- graph-based clustering
- differential marker analysis
- manual immune-cell annotation
- doublet detection
- sub-clustering
- pseudotime analysis
- gene-signature scoring
- biologically cautious interpretation of exploratory analyses

## Interpretation and limitations

This project uses a small public reference dataset and is intended as a methodological single-cell analysis project rather than a discovery cohort study.

Important limitations include:

- PBMC 3k is a reference dataset rather than an independent biological cohort;
- cell-type annotation is primarily transcriptome-based;
- UMAP is a visualization and should not be interpreted as a quantitative biological distance;
- pseudotime is an inferred ordering, not real time;
- ligand-receptor co-expression does not prove cellular communication;
- fine immune-cell states should be validated with additional markers, reference mapping or complementary experimental information.

## Reproducibility principles

The workflow was designed around explicit analytical choices and a fixed random seed. For a fully reproducible scRNA-seq project, the following should always be recorded with the analysis:

- R version
- package versions
- random seed
- QC thresholds
- number of principal components
- clustering resolution
- session information

## Selected references

1. Stuart T. et al. **Comprehensive Integration of Single-Cell Data.** Cell, 2019.  
   https://doi.org/10.1016/j.cell.2019.05.031

2. Hao Y. et al. **Integrated analysis of multimodal single-cell data.** Cell, 2021.  
   https://doi.org/10.1016/j.cell.2021.04.048

3. Hafemeister C., Satija R. **Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.** Genome Biology, 2019.  
   https://doi.org/10.1186/s13059-019-1874-1

4. Germain P.-L. et al. **Doublet identification in single-cell sequencing data using scDblFinder.** F1000Research, 2021.  
   https://doi.org/10.12688/f1000research.73600.1

5. Street K. et al. **Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics.** BMC Genomics, 2018.  
   https://doi.org/10.1186/s12864-018-4772-0

## Author

**Nahi El Akoum**  
Master 2 — Artificial Intelligence and Data Analysis in Biology (AIDA)  
Sorbonne Université
