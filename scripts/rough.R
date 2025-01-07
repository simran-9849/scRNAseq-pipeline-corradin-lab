library(tidyverse)


#Dimentionality reduction: PCA and UMAP
seurat_qc <- RunPCA(seurat_qc, assay = "SCT")
seurat_qc <- RunUMAP(seurat_qc, dims = 1:30, reduction = "pca")
seurat_qc <- FindNeighbors(seurat_qc, dims = 1:30, reduction = "pca")
seurat_qc <- FindClusters(seurat_qc, resolution = 0.5)
saveRDS(seurat_qc, "seuratqc-z.RDS")

#Visualize clusters
DimPlot(seurat_qc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(seurat_qc, reduction = "pca", dims= c(1,2), group.by = "seurat_clusters", label = TRUE)
ElbowPlot(seurat_qc, ndims = 50)

#Feature plot
FeaturePlot(integrated_seurat, features = c("PAPPA2", "PLS3", "PDE1A", "PRR16"),
            reduction = "umap")
FeaturePlot(seurat_qc, features = c("GAP43", "STMN2", "TUBA1A"), reduction = "pca", dims = c(1, 2))

VlnPlot(seurat_object, features = c("BCl11B", "GAD67", "PPP1R1B", "DRD1", 
                                "DRD2", "SIX3", "EBF1", "PDE10a", "FOXP1"), group.by = "orig.ident", split.by = "seurat_clusters")

# First, prepare the SCT assay for marker detection
integrated_seurat <- PrepSCTFindMarkers(integrated_seurat)

# Then run FindAllMarkers
markers <- FindAllMarkers(integrated_seurat,
                          assay = "SCT",
                          only.pos = TRUE, 
                          min.pct = 0.025, 
                          logfc.threshold = 0.025
)

top_markers <- markers %>% group_by() %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(seurat_object, features = unique(top_markers$gene)) + RotatedAxis()

write.table(doublet_counts, file = "doublet_counts.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# Visualize doublets
DimPlot(combined_seurat, group.by = combined_seurat@meta.data$DF.classifications_0.25_0.16_133)

#########################################################################

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
seurat <- PercentageFeatureSet(seurat, pattern = "^RP[SL]", col.name = "percent.rb")

#Filter genes and cells 
seurat <- subset(seurat, 
                    subset = nCount_RNA > 500 &
                      nFeature_RNA > 200 &
                      nCount_RNA < 2500 & 
                      percent.mt < 5)

object.list <- SplitObject(seurat, split.by = "orig.ident")

object.list <- lapply(X = object.list, FUN = function(x) {
  x <- SCTransform(x, 
                   vst.flavor = "v2",
                   vars.to.regress = c("percent.mt", "percent.rb"),
                   verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = object.list)
anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
integrated_seurat <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated_seurat) <- "SCT"
integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = FALSE)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

# Plot UMAP
DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

object.list <- lapply(object.list, function(x) {
  x <- RunPCA(x, assay = "SCT")
  x <- RunUMAP(x, dims = 1:30, reduction = "pca")
  x <- FindNeighbors(x, dims = 1:30, reduction = "pca")
  x <- FindClusters(x, resolution = 0.5)
  return(x)
})

umap_plots <- lapply(names(object.list), function(sample_name) {
  DimPlot(object.list[[sample_name]], reduction = "umap", group.by = "seurat_clusters") +
    ggtitle(sample_name)
})

wrap_plots(umap_plots, ncol = 2)

VlnPlot(seurat_object, features = c("BCL11B", "GAD67", "PPP1R1B", "DRD1", 
                          "DRD2", "SIX3", "EBF1", "PDE10A", "FOXP1"),
          group.by = "seurat_clusters",
        split.by = "orig.ident")

wrap_plots(plot_list, ncol = 1)

# Combine results into a single Seurat object
combined_seurat <- merge(processed_samples[[1]]$doublet_finder, 
                         y = lapply(processed_samples[-1], function(x) x$doublet_finder))

# Save the combined doublet finder results
saveRDS(combined_seurat, snakemake@output[["doublet_rds"]])

# Save plots to a PDF file for all samples
pdf(file.path(dirname(snakemake@output[["doublet_counts"]]), "doubletfinder_plots.pdf"))
for (sample in processed_samples) {
  print(sample$pk_plot)
  print(sample$doublet_plot)
}
dev.off()

# Combine and save doublet counts
doublet_counts <- do.call(rbind, lapply(processed_samples, function(x) x$doublet_table))
write.table(doublet_counts, file = snakemake@output[["doublet_counts"]], 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Filter out doublets and save the filtered Seurat object
df_cols <- grep("DF.classifications", colnames(combined_seurat@meta.data), value = TRUE)
singlet_cells <- rownames(combined_seurat@meta.data[combined_seurat@meta.data[, df_cols] == "Singlet", ])
doublet_filtered_seurat <- subset(combined_seurat, cells = singlet_cells)
saveRDS(doublet_filtered_seurat, snakemake@output[["filtered_seurat"]])

#Scatter plots for QC validation
scatter_plots <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.rb") +
  FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                 + geom_smooth(method = 'lm'))

