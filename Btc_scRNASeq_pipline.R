##############################################################################
# Btc knockout mouse scRNAseq analysis pipeline
#
# Pipeline: CellRanger multi output -> Seurat -> DecontX ambient RNA removal
# -> QC filtering -> normalisation -> PCA -> Harmony integration -> UMAP/
# clustering -> cluster annotation -> WT vs Mut differential expression
# (per-cell and pseudobulk) -> GO/KEGG functional enrichment
#
# All input/output paths below use a generic C:/Data/ structure - replace
# with your own local paths before running.
##############################################################################

# ---------------------------------------------------------------------------
# 0. Install packages (run once, then comment out)
# ---------------------------------------------------------------------------
# install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "tidyverse",
#                     "clustree", "RColorBrewer", "remotes", "R.utils",
#                     "harmony", "reshape2", "ggrepel", "devtools"))
# devtools::install_github("immunogenomics/presto")
#
# install.packages("BiocManager")
# BiocManager::install(c("SingleR", "celldex", "decontX", "BiocGenerics",
#                         "DelayedArray", "DelayedMatrixStats", "limma",
#                         "S4Vectors", "SingleCellExperiment",
#                         "SummarizedExperiment", "edgeR", "fgsea", "DESeq2",
#                         "clusterProfiler", "msigdbr", "org.Mm.eg.db"))

# ---------------------------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(SingleCellExperiment)
library(decontX)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(reshape2)
library(clustree)
library(harmony)
library(edgeR)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)

##############################################################################
# 2. Load data and build Seurat object
##############################################################################

setwd("C:/Data/")

# CellRanger multi filtered feature-barcode matrix (Gene Expression + CMO)
data <- Read10X(data.dir = "C:/Data/raw_data/Data_Aggr/count/filtered_feature_bc_matrix/")
# ls(data) -> "Gene Expression", "Multiplexing Capture"

seurat_object <- CreateSeuratObject(counts = data[["Gene Expression"]],
                                     min.cells = 3, min.features = 200,
                                     project = "Btc")

# Add CMO (cell multiplexing oligo) assay for sample demultiplexing
seurat_object[["CMO"]] <- CreateAssayObject(
  counts = data[["Multiplexing Capture"]][, colnames(x = seurat_object)]
)

##############################################################################
# 3. Sample demultiplexing (CMO tag calls -> Sample_ID, Sex, Condition)
##############################################################################

# One entry per CellRanger multi run; maps CMO tags to sample IDs, sex and
# genotype (Condition), and removes CMO-called doublets (num_features > 1)
sample_info <- list(
  list(file = "C:/Data/raw_data/1593--_CMO_303_304/multi/multiplexing_analysis/tag_calls_per_cell.csv",
       mapping = c("CMO303" = "KO1593_03", "CMO304" = "KO1593_04"),
       sex = c("M", "M"), condition = "Mut"),
  list(file = "C:/Data/raw_data/1649--_CMO_301_302/multi/multiplexing_analysis/tag_calls_per_cell.csv",
       mapping = c("CMO301" = "KO1649_01", "CMO302" = "KO1649_02"),
       sex = c("M", "M"), condition = "Mut"),
  list(file = "C:/Data/raw_data/1593--_CMO_305_306/multi/multiplexing_analysis/tag_calls_per_cell.csv",
       mapping = c("CMO305" = "KO1593_05", "CMO306" = "KO1593_06"),
       sex = c("M", "M"), condition = "Mut"),
  list(file = "C:/Data/raw_data/1579++_CMO_309_310/multi/multiplexing_analysis/tag_calls_per_cell.csv",
       mapping = c("CMO309" = "WT1579_09", "CMO310" = "WT1579_10"),
       sex = c("M", "F"), condition = "WT"),
  list(file = "C:/Data/raw_data/1579++_CMO_311_312/multi/multiplexing_analysis/tag_calls_per_cell.csv",
       mapping = c("CMO311" = "WT1579_11", "CMO312" = "WT1579_12"),
       sex = c("F", "F"), condition = "WT")
)

all_samples <- lapply(sample_info, function(info) {
  df <- read.csv(info$file)

  for (cm in names(info$mapping)) {
    df$sample[df$feature_call == cm] <- info$mapping[[cm]]
  }
  for (i in seq_along(info$mapping)) {
    df$Sex[df$sample == info$mapping[[i]]] <- info$sex[i]
  }
  df$Condition <- info$condition

  # Remove CMO-called doublets
  df <- df[df$num_features == 1, ]
  return(df)
})

combined <- do.call(rbind, all_samples)

# Sanity check: equal representation across Condition/Sex, no samples lost
table(combined$Condition)
table(combined$Sex)

seurat_object <- AddMetaData(seurat_object, combined$sample, col.name = "Sample_ID")
seurat_object <- AddMetaData(seurat_object, combined$Sex, col.name = "Sex")
seurat_object <- AddMetaData(seurat_object, combined$Condition, col.name = "Condition")

##############################################################################
# 4. Remove Y-chromosome genes (avoid sex-driven clustering bias)
##############################################################################

y_genes_df <- read.csv("Sundram_Btc_mouse_scRNAseq/Y_genes.csv", header = TRUE,
                        stringsAsFactors = FALSE)
y_genes <- y_genes_df[["Gene.name"]]

present_y <- intersect(y_genes, rownames(seurat_object[["RNA"]]))
message(length(present_y), " of ", length(y_genes),
        " Y-genes found in data and will be excluded")

seurat_object[["RNA"]] <- subset(
  x = seurat_object[["RNA"]],
  features = setdiff(rownames(seurat_object[["RNA"]]), present_y)
)

##############################################################################
# 5. QC: mitochondrial content filtering
##############################################################################

seurat_object$percent.mito <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
VlnPlot(seurat_object, features = "percent.mito", group.by = "Condition")

# Distribution is heavily weighted <10%, trailing off after ~20-25%. Given
# samples are early postnatal (PND3/4, naturally higher mitochondrial
# activity), a 5% cutoff was used (retains the large majority of cells).
seurat_object <- subset(seurat_object, subset = percent.mito < 5)

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 6. Ambient RNA correction (DecontX)
##############################################################################

seurat_object <- readRDS("C:/Data/scRNAseq/seurat_object.Rds")

# Rough clustering on raw counts to give DecontX cluster labels to model
# contamination against (coarse resolution - only major cell types needed)
seurat_tmp <- NormalizeData(seurat_object)
seurat_tmp <- FindVariableFeatures(seurat_tmp, nfeatures = 2000)
seurat_tmp <- ScaleData(seurat_tmp)
seurat_tmp <- RunPCA(seurat_tmp, verbose = FALSE)
ElbowPlot(seurat_tmp)  # plateaus ~PC12
seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:12)

for (i in seq(0.2, 2, by = 0.2)) {
  seurat_tmp <- FindClusters(seurat_tmp, resolution = i, verbose = FALSE)
}
clustree(seurat_tmp, prefix = "RNA_snn_res.")
seurat_tmp <- FindClusters(seurat_tmp, resolution = 0.2)

# Run DecontX on raw counts using the coarse cluster labels
sce <- as.SingleCellExperiment(seurat_object)
sce$cluster_for_decontx <- as.character(Idents(seurat_tmp))
sce <- decontX(sce, z = sce$cluster_for_decontx)
saveRDS(sce, file = "C:/Data/scRNAseq/sce.Rds")

contam_frac <- colData(sce)$decontX_contamination
corrected_counts <- decontXcounts(sce)

# Check haemoglobin contamination before/after DecontX correction
haem_genes <- c("Hbb-bt", "Hbb-bs", "Hbb-y", "Hba-a1", "Hbq1b", "Hba-a2", "Hbq1a")
haem_genes <- haem_genes[haem_genes %in% rownames(sce)]

haem_orig <- colSums(counts(sce)[haem_genes, , drop = FALSE]) / colSums(counts(sce))
haem_corrected <- colSums(corrected_counts[haem_genes, , drop = FALSE]) / colSums(corrected_counts)

haem_df <- data.frame(Cell = colnames(sce), Original = haem_orig, Corrected = haem_corrected)
haem_melt <- melt(haem_df, id.vars = "Cell", variable.name = "Stage", value.name = "HaemFraction")

ggplot(haem_melt, aes(x = HaemFraction, fill = Stage)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  scale_x_log10() +
  theme_classic() +
  labs(title = "Haemoglobin fraction per cell before/after DecontX",
       x = "Haemoglobin fraction (log scale)", y = "Number of cells")

# Remove highly contaminated cells (>30% estimated contamination)
sce_filtered <- sce[, contam_frac < 0.3]

orig_cells <- ncol(seurat_object)
filtered_cells <- ncol(sce_filtered)
cat("Cells retained after DecontX filtering:", filtered_cells, "/", orig_cells,
    "(", round(filtered_cells / orig_cells, 3), ")\n")

# Keep only DecontX-passing cells and replace counts with corrected counts
seurat_object <- subset(seurat_object, cells = colnames(sce_filtered))
seurat_object[["RNA"]] <- SetAssayData(
  seurat_object[["RNA"]], layer = "counts",
  new.data = decontXcounts(sce_filtered)
)

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 7. Normalisation and variable feature selection
##############################################################################

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 2000)
saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

top20 <- head(VariableFeatures(seurat_object), n = 20)
LabelPoints(VariableFeaturePlot(seurat_object), points = top20, repel = TRUE)

##############################################################################
# 8. Blood contamination removal
##############################################################################

blood_genes_mouse <- c(
  "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-b1", "Hbb-b2",
  "Lyz2", "S100a8", "S100a9", "C1qa", "C1qb", "C1qc",
  "Pf4", "Ppbp", "Ptprc", "Cd14", "Csf1r", "Mpo"
)
blood_genes_mouse <- intersect(blood_genes_mouse, rownames(seurat_object[["RNA"]]))

seurat_object <- AddModuleScore(seurat_object, features = list(blood_genes_mouse),
                                 name = "Blood")
VlnPlot(seurat_object, features = "Blood1", group.by = "Condition")

# Identify additional genes strongly correlated with the blood score
# (Spearman) so they can be excluded from downstream variable features
expr <- GetAssayData(seurat_object, slot = "data")
blood_score <- seurat_object$Blood1
cor_vals <- apply(expr, 1, function(x)
  cor(x, blood_score, method = "spearman", use = "pairwise.complete.obs"))
cor_p <- sapply(rownames(expr), function(g) {
  x <- expr[g, ]
  if (sd(x) == 0) return(NA)
  cor.test(x, blood_score, method = "spearman")$p.value
})
cor_df <- data.frame(gene = rownames(expr), rho = cor_vals, pval = cor_p,
                      padj = p.adjust(cor_p, "BH"))
suspect_genes <- cor_df$gene[abs(cor_df$rho) > 0.750 & cor_df$padj < 0.05]

VariableFeatures(seurat_object) <- setdiff(VariableFeatures(seurat_object), suspect_genes)

# Remove blood marker genes from the object entirely
genes_to_keep <- setdiff(rownames(seurat_object), blood_genes_mouse)
seurat_object <- subset(seurat_object, features = genes_to_keep)

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 9. Scaling and PCA
##############################################################################

seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object$pca, dims = 1:5, nfeatures = 5)

ElbowPlot(seurat_object)  # plateaus after ~PC12

# JackStraw to confirm which PCs are statistically meaningful vs noise
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:14)
JackStrawPlot(seurat_object, dims = 1:14, reduction = "pca")

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 10. UMAP and clustering resolution selection (pre-Harmony)
##############################################################################

seurat_object <- RunUMAP(seurat_object, dims = 1:14)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:14)

for (i in seq(0.2, 2, by = 0.2)) {
  seurat_object <- FindClusters(seurat_object, resolution = i, verbose = FALSE)
}
clustree(seurat_object, prefix = "RNA_snn_res.")

# Resolution 0.4 gave the most robust, consistent clustering across the
# clustree comparison
seurat_object <- FindClusters(seurat_object, resolution = 0.4)
seurat_object$pca_clusters <- seurat_object$seurat_clusters

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 11. Harmony integration
##############################################################################

DimPlot(seurat_object, reduction = "umap", group.by = "Condition")
DimPlot(seurat_object, reduction = "umap", group.by = "Sample_ID")
# No obvious batch-driven separation by Condition/Sex/Sample prior to Harmony

seurat_object <- RunHarmony(seurat_object, c("Condition", "Sex", "Sample_ID"),
                             reduction = "pca", reduction.save = "harmony")

seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:14,
                          reduction.name = "umap_harmony")
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:14)
seurat_object <- FindClusters(seurat_object, resolution = 0.4)
seurat_object$harmony_clusters <- seurat_object$seurat_clusters

clustree(seurat_object, prefix = "RNA_snn_res.", show_axis = TRUE) +
  theme(legend.key.size = unit(0.05, "cm"))

table(seurat_object$seurat_clusters)  # cluster sizes (all >200 cells except one)

# Confirm Harmony was necessary: adjusted R^2 of each PC vs Sample_ID was
# consistently very low (no strong batch effect), but Harmony embedding
# was used for downstream clustering/UMAP for consistency
pca_embed <- Embeddings(seurat_object, "pca")[, 1:14]
sapply(colnames(pca_embed), function(pc) {
  summary(lm(pca_embed[, pc] ~ seurat_object$Sample_ID))$adj.r.squared
})

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 12. Cluster marker genes and cell annotation.
##############################################################################

Idents(seurat_object) <- seurat_object$harmony_clusters
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE,
                                         min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers_sign <- subset(seurat_object.markers, p_val_adj < 0.05)
write.csv(seurat_object.markers_sign,
          file = "C:/Data/scRNAseq/seurat_object_ClusterMarkers.csv")

Top50genes <- seurat_object.markers %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)
write.csv(Top50genes, file = "C:/Data/scRNAseq/Top50avglog2FCgenes_per_cluster.csv")

# Cells anotations are assigned manually from the top 50 marker genes 
# using GPTCelltype https://github.com/Winnie09/GPTCelltype 

##############################################################################
# 13. Sub-clustering heterogeneous clusters. (optional)
# Clusters 0, 1 2, and 3 were large populations that cluster marker genes 
# suggested representative of mixed cell types. Subclustering was performed 
# on these clusters to resolve these populations. Replace cluster numbers for 
# sub-clustering or omit sub-clustering depending on your data.
##############################################################################

if (method == "kmeans") {
  res0 <- subcluster_kmeans(seurat_object, parent_cluster = "0", k = k_for_kmeans, pcs = pcs_to_use, seed = set.seed)
  res1 <- subcluster_kmeans(seurat_object, parent_cluster = "1", k = k_for_kmeans, pcs = pcs_to_use, seed = set.seed)
  res2 <- subcluster_kmeans(seurat_object, parent_cluster = "2", k = k_for_kmeans, pcs = pcs_to_use, seed = set.seed)
  res3 <- subcluster_kmeans(seurat_object, parent_cluster = "3", k = k_for_kmeans, pcs = pcs_to_use, seed = set.seed)
} else {
  res0 <- subcluster_seurat(seurat_object, parent_cluster = "0", resolution = seurat_resolution, pcs = pcs_to_use)
  res1 <- subcluster_seurat(seurat_object, parent_cluster = "1", resolution = seurat_resolution, pcs = pcs_to_use)
  res2 <- subcluster_seurat(seurat_object, parent_cluster = "2", resolution = seurat_resolution, pcs = pcs_to_use)
  res3 <- subcluster_seurat(seurat_object, parent_cluster = "3", resolution = seurat_resolution, pcs = pcs_to_use)
}

# quick sanity print
message("Cluster 0 split -> labels (example): ", paste(unique(res0$cell_labels)[1:min(5,length(unique(res0$cell_labels)))], collapse = ", "))
message("Cluster 1 split -> labels (example): ", paste(unique(res1$cell_labels)[1:min(5,length(unique(res1$cell_labels)))], collapse = ", "))
message("Cluster 2 split -> labels (example): ", paste(unique(res2$cell_labels)[1:min(5,length(unique(res2$cell_labels)))], collapse = ", "))
message("Cluster 3 split -> labels (example): ", paste(unique(res3$cell_labels)[1:min(5,length(unique(res3$cell_labels)))], collapse = ", "))

# Inject new labels into harmony_clusters column (but preserve original)

seurat_object$harmony_clusters_not_subcl <- seurat_object$harmony_clusters

orig_col   <- seurat_object$harmony_clusters
old_levels <- levels(seurat_object$harmony_clusters)
new_subcl_col <- orig_col
names(new_subcl_col) <- rownames(seurat_object@meta.data)

# assign new labels to the cells that were in cluster 0/1/2/3/etc
new_subcl_col[names(res0$cell_labels)] <- res0$cell_labels[names(res0$cell_labels)]
new_subcl_col[names(res1$cell_labels)] <- res1$cell_labels[names(res1$cell_labels)]
new_subcl_col[names(res2$cell_labels)] <- res2$cell_labels[names(res2$cell_labels)]
new_subcl_col[names(res3$cell_labels)] <- res3$cell_labels[names(res3$cell_labels)]

# write the combined labels back into harmony_clusters (overwrites original)
new_levels <- sort(unique(c(old_levels, unique(new_subcl_col))))
seurat_object$harmony_clusters <- factor(new_subcl_col, levels = new_levels)
table(seurat_object$harmony_clusters)

Idents(seurat_object) <- seurat_object$harmony_clusters

# Identifying sub-cluster marker genes and cell annotation
seurat_object.markers_subcl <- FindAllMarkers(seurat_object, only.pos = TRUE,
                                              min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers_subcl_sign <- subset(seurat_object.markers_subcl, p_val_adj < 0.05)
write.csv(seurat_object.markers_subcl_sign,
          file = "C:/Data/scRNAseq/seurat_object_ClusterMarkers_subcl.csv")
Top50genes_subcl <- seurat_object.markers_subcl %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)
write.csv(Top50genes_subcl, file = "C:/Data/scRNAseq/Top50avglog2FCgenes_per_cluster_subcl.csv")

# Cells anotations are assigned manually to cluster 0, 1, 2 and 3 from the top 50 marker genes
# after sub-clustering using GPTCelltype https://github.com/Winnie09/GPTCelltype

##############################################################################
# 14. Cluster annotation
##############################################################################

cluster_annotations_long <- c(
  "0a" = "Homeostatic Microglia / brain-resident macrophage (activated/inflammatory-like)",
  "0b" = "Activated microglia / inflammatory macrophages",
  "1a" = "Radial glia",
  "1b" = "Radial glia / early astroglial lineage",
  "2a" = "CNS endothelial / blood-brain-barrier-like cells",
  "2b" = "Immature astrocytes",
  "3a" = "Excitatory projection neurons (likely cortical layer-specific)",
  "3b" = "GABAergic interneurons (Dlx/Gad2+, migrating)",
  "4"  = "Cajal-Retzius / Reln-expressing marginal zone neurons",
  "5"  = "Oligodendrocyte precursor cells (OPCs)",
  "6"  = "Erythroid / RBC contamination",
  "7"  = "Border-associated macrophages / perivascular macrophages",
  "8"  = "Glycinergic inhibitory neurons (likely hindbrain/cerebellar)",
  "9"  = "Meningeal / stromal fibroblast-like cells",
  "10" = "Proliferating/mitotic progenitors",
  "11" = "Proliferating neural progenitors",
  "12" = "Glutamatergic excitatory neurons (Slc17a6/VGLUT2+)",
  "13" = "Interneuron progenitors / migrating interneurons (Dlx1/2)",
  "14" = "Myeloid / infiltrating granulocyte-like cells",
  "15" = "Erythroid / RBC contamination",
  "16" = "Meningeal / mesenchymal fibroblast-like cells",
  "17" = "B lymphocytes (Ig genes, Cd79a, Ms4a1)",
  "18" = "Choroid plexus epithelial / ependymal-like cells",
  "19" = "Brain endothelial / BBB transporter-rich cells"
  
)
cluster_annotations_short <- c(
  "0a" = "Homeostatic microglia",
  "0b" = "Activated microglia",
  "1a" = "Radial glia",
  "1b" = "Radial glia early astrocyte lineage",
  "2a" = "Brain endothelial cells",
  "2b" = "Immature astrocytes",
  "3a" = "Excitatory projection neurons",
  "3b" = "GABAergic interneurons",
  "4"  = "Cajal-Retzius neurons",
  "5"  = "Oligodendrocyte precursor cells",
  "6"  = "Erythroid cells",
  "7"  = "Border-associated macrophages",
  "8"  = "Glycinergic inhibitory neurons",
  "9"  = "Meningeal fibroblast-like cells",
  "10" = "Proliferating progenitors",
  "11" = "Proliferating neural progenitors",
  "12" = "Glutamatergic excitatory neurons",
  "13" = "Interneuron progenitors",
  "14" = "Myeloid inflammatory cells",
  "15" = "Erythroid cells",
  "16" = "Meningeal mesenchymal cells",
  "17" = "B lymphocytes",
  "18" = "Choroid plexus epithelial cells",
  "19" = "Brain endothelial transporter-rich cells"
)
# NOTE: update the mapping above if your final cluster count labels differ.

Idents(seurat_object) <- "harmony_clusters"
seurat_object$celltype_short <- plyr::mapvalues(
  x = as.character(Idents(seurat_object)),
  from = names(cluster_annotations_short), to = as.vector(cluster_annotations_short)
)
seurat_object$celltype_long <- plyr::mapvalues(
  x = as.character(Idents(seurat_object)),
  from = names(cluster_annotations_long), to = as.vector(cluster_annotations_long)
)

cluster_colors <- setNames(
  c("#F8766D","#EE8043","#E18A00","#D19300","#BE9C00","#A8A400","#EB69F0",
    "#68B100","#D575FE","#00BB49","#00BE70","#00C090","#00C1AB","#00BFC4",
    "#00BBDA","#00B5ED","#00ACFC","#42A0FF","#8B93FF","#B684FF","#24B700",
    "#8CAB00","#F962DD","#FF61C6"),
  
    levels(factor(seurat_object$celltype_short))
)

# UMAP with annotated cluster centroids labelled
umap_coords <- as.data.frame(Embeddings(seurat_object, "umap_harmony"))
colnames(umap_coords) <- c("Dim1", "Dim2")
umap_coords$short <- seurat_object$celltype_short
centroids <- umap_coords %>% group_by(short) %>%
  summarise(Dim1 = mean(Dim1), Dim2 = mean(Dim2))

DimPlot(seurat_object, reduction = "umap_harmony", group.by = "celltype_long",
        cols = cluster_colors, label = FALSE, repel = TRUE) +
  ggrepel::geom_text_repel(data = centroids, aes(x = Dim1, y = Dim2, label = short),
                            color = "black", size = 4, fontface = "bold",
                            box.padding = 1, max.overlaps = Inf) +
  ggtitle("UMAP of annotated clusters") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none")

saveRDS(seurat_object, file = "C:/Data/scRNAseq/seurat_object.Rds")

##############################################################################
# 15. Low-expression gene filtering
##############################################################################

total_per_gene <- rowSums(GetAssayData(seurat_object, assay = "RNA", slot = "counts"))
seurat_object <- seurat_object[total_per_gene >= 50, ]

##############################################################################
# 16. Sex-linked gene check (confirms CMO-based sex assignment)
##############################################################################

FeaturePlot(seurat_object, features = c("Ddx3y", "Uty", "Eif2s3x"),
            split.by = "Sample_ID")

##############################################################################
# 17. Differential expression: WT vs Mut per cluster (single-cell level)
##############################################################################

output_dir_clusters <- "C:/Data/scRNAseq/DEGs_by_cluster"

Idents(seurat_object) <- seurat_object$harmony_clusters
clusters <- levels(seurat_object)

for (cluster in clusters) {
  seurat_cluster <- seurat_object[, seurat_object@active.ident == cluster]
  celltype <- unique(seurat_cluster$celltype_short)
  if (length(celltype) > 1) celltype <- paste(celltype, collapse = "_")

  Idents(seurat_cluster) <- seurat_cluster$Condition
  de_result <- FindMarkers(seurat_cluster, ident.1 = "WT", ident.2 = "Mut",
                            logfc.threshold = 0.25, min.cells.feature = 3,
                            min.cells.group = 3, min.pct = 0, only.pos = FALSE)
  de_result$AveExpr <- rowMeans(
    GetAssayData(seurat_cluster, assay = "RNA", layer = "data")[rownames(de_result), ]
  )

  significant <- de_result[!is.na(de_result$p_val_adj) & de_result$p_val_adj < 0.05, ]
  if (nrow(significant) == 0) {
    message("No significant DEGs for cluster ", cluster, " (", celltype, "), skipping.")
    next
  }

  celltype_safe <- gsub("[^A-Za-z0-9_\\-]", "_", celltype)
  write.csv(significant,
            file.path(output_dir_clusters, paste0(cluster, "-", celltype_safe, ".csv")),
            row.names = TRUE)
}

##############################################################################
# 18. Differential expression: WT vs Mut, whole dataset
##############################################################################

Idents(seurat_object) <- seurat_object$Condition
ALL_DEGs <- FindMarkers(seurat_object, ident.1 = "Mut", ident.2 = "WT",
                         logfc.threshold = 0.25, min.cells.feature = 3,
                         min.cells.group = 3, min.pct = 0, only.pos = FALSE)
write.csv(subset(ALL_DEGs, p_val_adj < 0.05),
          file = "C:/Data/scRNAseq/ALL_DEGs.csv")

##############################################################################
# 19. Pseudobulk differential expression (limma-voom), WT vs Mut per cluster
##############################################################################

output_dir_pseudobulk <- "C:/Data/Btc_scRNAseq/Pseudobulk"
padj_cutoff <- 0.1  # CONFIRM: 0.05 was also used in earlier iterations

clusters <- levels(seurat_object)

for (cluster in clusters) {
  seurat_cluster <- seurat_object[, seurat_object$harmony_clusters == cluster]
  celltype <- unique(seurat_cluster$celltype_short)
  if (length(celltype) > 1) celltype <- paste(celltype, collapse = "_")
  
  counts <- GetAssayData(seurat_cluster, assay = "RNA", slot = "counts")
  meta <- seurat_cluster@meta.data
  
  # Sum counts per sample to create pseudobulk profiles
  pseudobulk_counts <- t(rowsum(t(counts), group = meta$Sample_ID))
  
  sample_info_pb <- meta %>% distinct(Sample_ID, .keep_all = TRUE)
  if (length(unique(sample_info_pb$Condition)) < 2) {
    message("Skipping cluster ", cluster, ": only one condition present")
    next
  }
  
  if ("WT" %in% sample_info_pb$Condition) {
    group <- relevel(factor(sample_info_pb$Condition), ref = "WT")
  } else {
    warning(paste("Cluster", cluster, "has no WT cells; using first level as reference"))
    group <- factor(sample_info_pb$Condition)
  }
  design <- model.matrix(~ group)
  
  y <- DGEList(counts = pseudobulk_counts)
  y <- calcNormFactors(y)
  v <- voom(y, design)
  fit <- eBayes(lmFit(v, design))
  
  degs <- topTable(fit, coef = "groupMut", number = Inf, sort.by = "P")
  degs$AveExpr <- rowMeans(v$E[rownames(degs), ])
  
  significant <- degs[degs$adj.P.Val < padj_cutoff, ]
  if (nrow(significant) == 0) {
    message("No significant pseudobulk DEGs for cluster ", cluster, " (", celltype, ")")
    next
  }
  write.csv(significant,
            file.path(output_dir_pseudobulk, paste0(cluster, "-", celltype, "-pseudobulk.csv")),
            row.names = TRUE)
}

# Check sample balance per cluster/condition (min samples needed for pseudobulk DE)
check_pseudobulk_balance <- function(seurat_object, cluster_col = "celltype_short",
                                     sample_col = "Sample_ID", condition_col = "Condition",
                                     min_samples_per_condition = 2) {
  meta <- seurat_object@meta.data
  summary_df <- meta %>%
    group_by(!!sym(cluster_col), !!sym(sample_col), !!sym(condition_col)) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  conds <- unique(meta[[condition_col]])
  cluster_condition_summary <- summary_df %>%
    group_by(!!sym(cluster_col), !!sym(condition_col)) %>%
    summarise(n_samples = n_distinct(!!sym(sample_col)), total_cells = sum(cell_count),
              .groups = "drop") %>%
    pivot_wider(names_from = !!sym(condition_col),
                values_from = c(n_samples, total_cells), values_fill = 0)
  
  col1 <- paste0("n_samples_", conds[1]); col2 <- paste0("n_samples_", conds[2])
  cluster_condition_summary %>%
    mutate(eligible_for_DE = !!sym(col1) >= min_samples_per_condition &
             !!sym(col2) >= min_samples_per_condition)
}

cluster_check <- check_pseudobulk_balance(seurat_object)
View(cluster_check)


##############################################################################
# 20. Functional enrichment (GO / KEGG) on pseudobulk DEGs
##############################################################################

deg_path <- output_dir_pseudobulk
padj_enrich_cutoff <- 0.05
logfc_enrich_cutoff <- 0.05

results_list_GO <- list()
results_list_KEGG <- list()

for (file in list.files(deg_path, full.names = TRUE)) {
  cluster_name <- tools::file_path_sans_ext(basename(file))
  deg <- read.csv(file)

  sig_genes <- deg %>%
    filter(adj.P.Val < padj_enrich_cutoff & abs(logFC) > logfc_enrich_cutoff) %>%
    pull(gene)

  if (length(sig_genes) == 0) {
    message("No DEGs passed thresholds for cluster ", cluster_name)
    next
  }

  go_enr <- tryCatch(
    enrichGO(gene = sig_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) NULL
  )
  if (!is.null(go_enr)) {
    go_df <- go_enr@result
    go_df$cluster <- cluster_name
    results_list_GO[[cluster_name]] <- go_df
  }

  kegg_enr <- tryCatch(
    enrichKEGG(gene = sig_genes, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05),
    error = function(e) NULL
  )
  if (!is.null(kegg_enr)) {
    kegg_df <- kegg_enr@result
    kegg_df$cluster <- cluster_name
    results_list_KEGG[[cluster_name]] <- kegg_df
  }
}

write.csv(bind_rows(results_list_GO),
          "C:/Data/scRNAseq/GO_BP_enrichment_all_clusters.csv", row.names = FALSE)
write.csv(bind_rows(results_list_KEGG),
          "C:/Data/scRNAseq/KEGG_enrichment_all_clusters.csv", row.names = FALSE)

# Hallmark gene set enrichment (MSigDB) as an alternative/complementary check
msig_h <- msigdbr(species = "Mus musculus", category = "H")
results_list_hallmark <- list()

for (file in list.files(deg_path, full.names = TRUE)) {
  cluster_name <- tools::file_path_sans_ext(basename(file))
  deg <- read.csv(file)
  sig_genes <- deg %>%
    filter(adj.P.Val < padj_enrich_cutoff & abs(logFC) > logfc_enrich_cutoff) %>%
    pull(gene)
  if (length(sig_genes) == 0) next

  enr <- enricher(sig_genes, TERM2GENE = msig_h %>% dplyr::select(gs_name, gene_symbol),
                   pvalueCutoff = 0.05)
  enr@result$cluster <- cluster_name
  results_list_hallmark[[cluster_name]] <- enr@result
}

write.csv(bind_rows(results_list_hallmark),
          "C:/Data/scRNAseq/Hallmark_enrichment_all_clusters.csv", row.names = FALSE)
