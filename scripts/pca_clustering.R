library(Seurat)
library(ggplot2)
options(ggrepel.max.overlaps = Inf) #sets maximum overlap to inf (use case:Cluster plot visualization)

#Dimentionality reduction: PCA and UMAP
seurat_qc <- RunPCA(seurat_qc)
seurat_qc <- RunUMAP(seurat_qc, dims = 1:30, reduction = "pca")
seurat_qc <- FindNeighbors(seurat_qc, dims = 1:30, reduction = "pca")
seurat_qc <- FindClusters(seurat_qc, resolution = 0.5)


# Save plots
# PCA plot (Elbow Plot)
elbow_plot <- ElbowPlot(seur, ndims = 30)
ggsave(filename = file.path(plotsDir, paste0(sampleName, "_PCAElbowPlot.png")), plot = elbow_plot, width = 6, height = 4)

# UMAP plot
umap_plot <- DimPlot(seurat_qc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP Clustering")
ggsave(filename = file.path(plotsDir, paste0(sampleName, "_UMAP.png")), plot = umap_plot, width = 6, height = 6)

# PCA Feature Plot 
pca_plot <- DimPlot(seurat_qc, reduction = "pca", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("PCA Clustering")
ggsave(filename = file.path(plotsDir, paste0(sampleName, "_PCA.png")), plot = pca_plot, width = 6, height = 6)

