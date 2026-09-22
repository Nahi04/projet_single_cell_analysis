# ==========================
# single cell analysis
# ==========================

# 1. Installation des packages requis 
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SeuratData")

install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
library(SeuratData)


# charger les packages utiles
suppressPackageStartupMessages({ library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(dplyr)
  library(patchwork) })


options(future.globals.maxSize = 8000 * 1024^2)

# un seed pour la reproductibilité
set.seed(42)



# ==================================
# telecharger le jeu de données
# ==================================

# Vérification et téléchargement du dataset PBMC 3k/10k
AvailableData() 

# Installation et chargement du dataset pbmc3k 
# InstallData() + LoadData() créent directement un objet Seurat prêt à l’emploi. c'est pour cela qu'il faut que seuratdata soit installé
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k")

cat("Nombre de gènes initiaux :", nrow(pbmc), "\n")
cat("Nombre de cellules initiales :", ncol(pbmc), "\n")

class(pbmc)  # donne seurat objext , c'est une étape de verification avant de continuer
dim(pbmc)



# ============================
# quality control 
#=============================

# mito et ribosomes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# visualiser quelques metriques
VlnPlot( pbmc, 
         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
         ncol = 3, 
         pt.size = 0.1
) & theme_minimal()


# Filtrage selon la méthodologie standard Seurat sur PBMC
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# nombre de ceullules apres filtrage
cat("Nombre de cellules restantes :", ncol(pbmc), "\n")





# =========
# eliminer les doublet
# =========
install.packages("remotes")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")
}
library(scDblFinder)
library(SingleCellExperiment)

# 2. Exécution de scDblFinder (Conversion temporaire SCE)
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, verbose = FALSE)

# 3. Injection des métriques dans l'objet Seurat
pbmc$Doublet_Status <- sce$scDblFinder.class
pbmc$Doublet_Score  <- sce$scDblFinder.score

# 4. Normalisation et UMAP temporaires pour la représentation graphique
pbmc_visu <- NormalizeData(pbmc, verbose = FALSE)
pbmc_visu <- FindVariableFeatures(pbmc_visu, verbose = FALSE)
pbmc_visu <- ScaleData(pbmc_visu, verbose = FALSE)
pbmc_visu <- RunPCA(pbmc_visu, verbose = FALSE)
pbmc_visu <- RunUMAP(pbmc_visu, dims = 1:10, verbose = FALSE)

# 5. Visualisation des Résultats (2 Graphiques)
# 1. Visualisation sans viridis (utilisation de scale_color_gradient)
p1 <- DimPlot(
  pbmc_visu, 
  group.by = "Doublet_Status", 
  cols = c("singlet" = "steelblue", "doublet" = "firebrick")
) + ggtitle("Classification des Doublets") + theme_minimal()

p2 <- FeaturePlot(pbmc_visu, features = "Doublet_Score") + 
  scale_color_gradient(low = "lightgray", high = "darkred") + 
  ggtitle("Score de Probabilité de Doublet") + 
  theme_minimal()

p1 | p2

table(pbmc$Doublet_Status)

pbmc <- subset(pbmc, subset = Doublet_Status == "singlet")
cat("Cellules uniques (singlets) conservées :", ncol(pbmc), "\n")



# =============
# evaluer le cycle cellulaire 
# =============

# 1. herger les marqueurs
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# 2. Normalisation et scoring du cycle
pbmc <- NormalizeData(pbmc, verbose = FALSE)

pbmc <- CellCycleScoring( pbmc, 
                          s.features = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident = TRUE)

# Aperçu des scores attribués
head(pbmc@meta.data[, c("S.Score", "G2M.Score", "Phase")])

# 3. Intersecter les gènes du cycle présents dans la matrice pour éviter le crash
cc.genes.present <- intersect(c(s.genes, g2m.genes), rownames(pbmc))

# 4.FindVariableFeatures + ScaleData pour la PCA
pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
pbmc <- ScaleData(pbmc, features = rownames(pbmc), verbose = FALSE)

# 5. PCA ciblée sur les gènes du cycle cellulaire
pbmc <- RunPCA(pbmc, features = cc.genes.present, verbose = FALSE)

# 6. Affichage graphique propre
DimPlot(pbmc, reduction = "pca", group.by = "Phase") + 
  ggtitle("Structure du cycle cellulaire (PCA)") + 
  theme_minimal()




# ===============================
# application de sctransform
#================================
pbmc <- SCTransform( pbmc, 
                     vst.flavor = "v2", 
                     vars.to.regress = "percent.mt", 
                     verbose = FALSE )


# ==========================
# reduction de dimension
# ==========================
# 1. Calcul de la PCA sur l'assay SCT
pbmc <- RunPCA(pbmc, assay = "SCT", npcs = 50, verbose = FALSE)

# 2. Visualisation de la variance expliquée par chaque composante principale
ElbowPlot(pbmc, ndims = 50) + 
  ggtitle("Elbow Plot - Choix des dimensions PCA") + 
  theme_minimal()



# =====================================
# on va commencer le clustering
# =====================================

# 1. Construction du graphe K-NN / SNN basé sur les 20 premières PCs
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)

# 2. Clustering par algorithme de Louvain (Résolution de 0.5)
# on peut tester plusieurs resolution jusqu'a trouver un resultat qui semble etre optimal
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# 3. Projection non-linéaire UMAP 2D
pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)

# 4. Visualisation de l'UMAP avec les clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP - PBMC 10k (Clustering Louvain res = 0.5)") + 
  theme_minimal()



# 1. Préparation de l'assay pour l'analyse différentielle
PrepSCTFindMarkers(pbmc)

# 2. Recherche des marqueurs spécifiques à chaque cluster (uniquement les sur-exprimés)
pbmc.markers <- FindAllMarkers( pbmc, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25
)

# 3. Extraction du top 2 gènes marqueurs par cluster
library(dplyr)

top2_list <- lapply(split(pbmc.markers, pbmc.markers$cluster), function(x) {
  head(x[order(x$avg_log2FC, decreasing = TRUE), ], 2)
})

top2 <- do.call(rbind, top2_list)

# Affichage des résultats
print(top2[, c("cluster", "gene", "avg_log2FC", "p_val_adj")])


# 1. Définition du vecteur d'annotation des clusters
new.cluster.ids <- c(
  "0" = "CD4+ T Memory",
  "1" = "CD4+ T Naive",
  "2" = "CD14+ Monocytes",
  "3" = "CD8+ T Cells",
  "4" = "B Cells",
  "5" = "NK Cells",
  "6" = "FCGR3A+ Monocytes",
  "7" = "Dendritic Cells",
  "8" = "Platelets"
)

# 2. Re-nommage des identités dans l'objet Seurat
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$cell_type <- Idents(pbmc)

# 3. Visualisation de l'UMAP annoté
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + 
  ggtitle("UMAP - Annotations des Cellules Immunitaires (PBMC 10k)") + 
  theme_minimal() + 
  theme(legend.position = "right")




# Sélection de gènes marqueurs clés canoniques pour l'affichage
known_markers <- c(
  "CD40LG", "AQP3",    # CD4+ T Memory
  "CCR7", "SELL",      # CD4+ T Naive
  "S100A8", "CD14",    # CD14+ Monocytes
  "GZMK", "CD8A",      # CD8+ T Cells
  "MS4A1", "VPREB3",   # B Cells
  "GNLY", "NKG7",      # NK Cells
  "FCGR3A", "CDKN1C",  # FCGR3A+ Monocytes
  "FCER1A", "CST3",    # Dendritic Cells
  "PPBP", "CLDN5"      # Platelets
)

# DotPlot de publication
DotPlot(pbmc, features = known_markers, group.by = "cell_type") + 
  RotatedAxis() + 
  scale_color_gradient(low = "lightgrey", high = "darkblue") +
  ggtitle("Expression des marqueurs canoniques par type cellulaire") +
  theme_minimal()



# Calcul du nombre et pourcentage de cellules par type
cell_counts <- table(pbmc$cell_type)
cell_props  <- prop.table(cell_counts) * 100

df_summary <- data.frame(
  Count = as.numeric(cell_counts),
  Percentage = round(as.numeric(cell_props), 2),
  row.names = names(cell_counts)
)

print(df_summary)

# Sauvegarde du Seurat Object final
saveRDS(pbmc, file = "PBMC_10k_SCTv2_Annotated.rds")


# ==============
# sub clustering des compartiments t et nk
# ==============

# Sélection des cellules T et NK
t_nk_cells <- subset( pbmc, 
                      idents = c("CD4+ T Memory", "CD4+ T Naive", "CD8+ T Cells", "NK Cells"))

# SCTransform dédiée au sous-ensemble
t_nk_cells <- SCTransform( t_nk_cells, 
                           vst.flavor = "v2", 
                           vars.to.regress = "percent.mt", 
                           verbose = FALSE)

# Réduction de dimension dédiée
t_nk_cells <- RunPCA(t_nk_cells, npcs = 30, verbose = FALSE)
t_nk_cells <- FindNeighbors(t_nk_cells, dims = 1:15, verbose = FALSE)

# Resolution plus fine (0.8) pour capter les petites sous-populations (ex: Treg)
t_nk_cells <- FindClusters(t_nk_cells, resolution = 0.8, verbose = FALSE)
t_nk_cells <- RunUMAP(t_nk_cells, dims = 1:15, verbose = FALSE)

# Visualisation du sub-clustering
DimPlot(t_nk_cells, reduction = "umap", label = TRUE) + 
  ggtitle("Sub-clustering du Compartiment T / NK") + 
  theme_minimal()






install.packages(c("pkgconfig", "igraph", "tibble"), force = TRUE)
pbmc <- readRDS("PBMC_10k_SCTv2_Annotated.rds")

# trajectoire Slingshot
BiocManager::install("slingshot")
library(slingshot)

# Extraction des embeddings UMAP et des clusters
sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(pbmc, layer = "counts")),
  colData = pbmc@meta.data
)

# Trajectoire Slingshot
sce <- slingshot(sce, clusterLabels = 'cell_type', reducedDim = Embeddings(pbmc, "umap"), start.clus = "CD14+ Monocytes")

# Visualisation des courbes de trajectoire sur l'UMAP
plot(Embeddings(pbmc, "umap"), col = as.numeric(pbmc$cell_type), pch = 16, cex = 0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')


library(ggplot2)

# 1. Extraction des courbes Slingshot
curves <- SlingshotDataSet(sce)

# 2. Conversion des courbes en data.frame pour ggplot2
curve_coords <- lapply(slingCurves(curves), function(curve) {
  as.data.frame(curve$s[curve$ord, ])
})
curve_df <- do.call(rbind, lapply(seq_along(curve_coords), function(i) {
  cbind(curve_coords[[i]], Lineage = as.character(i))
}))

# 3. Superposition de la trajectoire sur l'UMAP Seurat
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + 
  geom_path(
    data = curve_df, 
    aes(x = umap_1, y = umap_2, group = Lineage), 
    color = "black", 
    linewidth = 1.2, 
    arrow = arrow(type = "closed", length = unit(0.15, "inches"))
  ) +
  ggtitle("Inférence de Trajectoire Slingshot (PBMC 10k)") + 
  theme_minimal()



# Ajout du Pseudotime dans les métadonnées de Seurat
pbmc$pseudotime <- slingPseudotime(sce)[, 1]

# FeaturePlot du Pseudotime
FeaturePlot(pbmc, features = "pseudotime", pt.size = 0.5) + 
  scale_color_viridis_c(option = "magma") +
  ggtitle("Progression du Pseudotime (Continuum de différenciation)") +
  theme_minimal()


library(slingshot)
library(SingleCellExperiment)
library(ggplot2)

# 1. Sous-ensemble du compartiment myéloïde
myeloid_pbmc <- subset(pbmc, idents = c("CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic Cells"))

# 2. Inférence Slingshot ciblée
sce_myeloid <- slingshot(
  as.SingleCellExperiment(myeloid_pbmc), 
  clusterLabels = 'cell_type', 
  reducedDim = Embeddings(myeloid_pbmc, "umap"), 
  start.clus = "CD14+ Monocytes"
)

# 3. Extraction et structuration des coordonnées de la trajectoire
curves_myeloid <- SlingshotDataSet(sce_myeloid)
curve_coords_myeloid <- lapply(slingCurves(curves_myeloid), function(curve) {
  as.data.frame(curve$s[curve$ord, ])
})

curve_df_myeloid <- do.call(rbind, lapply(seq_along(curve_coords_myeloid), function(i) {
  cbind(curve_coords_myeloid[[i]], Lineage = as.character(i))
}))

# 4. Visualisation propre et biologiquement rigoureuse
DimPlot(myeloid_pbmc, reduction = "umap", label = TRUE, pt.size = 0.8, repel = TRUE) + 
  geom_path(
    data = curve_df_myeloid, 
    aes(x = umap_1, y = umap_2, group = Lineage), 
    color = "black", 
    linewidth = 1.2, 
    arrow = arrow(type = "closed", length = unit(0.15, "inches"))
  ) +
  ggtitle("Trajectoire de Maturation Monocytaire (CD14+ -> FCGR3A+)") + 
  theme_minimal()



# =========================================
#interaction cellulaire 
# =========================================
# 1. Définition d'une liste de couples Ligand-Récepteur immunologiques canoniques
lr_pairs <- data.frame(
  Ligand = c("CD40LG", "CCL5", "GNLY", "TNF", "IL2"),
  Receptor = c("CD40", "CCR5", "CD244", "TNFRSF1A", "IL2RA")
)

# 2. Extraire la moyenne d'expression par type cellulaire dans Seurat
cluster_means <- AverageExpression(pbmc, assays = "SCT", layer = "data")$SCT

# 3. Filtrer sur les gènes de nos couples présents dans ton jeu de données
genes_present <- intersect(c(lr_pairs$Ligand, lr_pairs$Receptor), rownames(cluster_means))
sub_means <- cluster_means[genes_present, ]

# 4. Afficher la matrice d'expression des Ligands et Récepteurs par cluster
print(round(sub_means, 2))

# 5. Visualiser l'expression croisée sous forme de DotPlot Seurat
DotPlot(pbmc, features = genes_present, group.by = "cell_type") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  labs(title = "Communication Ligand-Récepteur (PBMC 10k)", x = "Gènes (L/R)", y = "Type Cellulaire")



BiocManager::install("AUCell")
library(AUCell)

# Extraction de la matrice de comptage normalisée
exprMatrix <- GetAssayData(pbmc, layer = "data")

# Calcul du ranking des gènes par cellule
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)

# Exemple : Évaluer l'activité de facteurs de transcription clés (ex: SPI1/PU.1 pour monocytes, PAX5 pour B cells)
tf_signatures <- list(
  SPI1_targets = c("CD14", "FCGR3A", "S100A8", "CSF1R"),
  PAX5_targets = c("MS4A1", "CD19", "VPREB3", "CD79A")
)

cells_AUC <- AUCell_calcAUC(tf_signatures, cells_rankings)

# Intégration de l'activité du régulon dans Seurat
pbmc$SPI1_activity <- as.numeric(getAUC(cells_AUC)["SPI1_targets", ])

# Visualisation sur UMAP
FeaturePlot(pbmc, features = "SPI1_activity", pt.size = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Activité du régulon SPI1 (PU.1) - Spécifique des Monocytes") +
  theme_minimal()
