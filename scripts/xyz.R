library(Seurat)
library(ggplot2)

# Read Seurat object
seurat <- readRDS(snakemake@input[["seurat"]])

# Process mitochondrial and ribosomal genes
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
seurat <- PercentageFeatureSet(seurat, pattern = "^RP[SL]", col.name = "percent.rb")

# Split the Seurat object
sample_list <- SplitObject(seurat, split.by = "orig.ident")

# Process each sample separately
for (sample_name in names(sample_list)) {
  sample_seurat <- sample_list[[sample_name]]
  
  # QC plots
  p1 <- VlnPlot(sample_seurat, features = "nFeature_RNA")
  p2 <- VlnPlot(sample_seurat, features = "nCount_RNA")
  p3 <- VlnPlot(sample_seurat, features = "percent.mt")
  p4 <- VlnPlot(sample_seurat, features = "percent.rb")
  
  qc_violin <- (p1 | p2) / (p3 | p4)
  
  # Save QC violin plot
  pdf(paste0(snakemake@output[["qc_violin"]], "_", sample_name, ".pdf"), width = 12, height = 10)
  print(qc_violin + plot_annotation(title = paste0("MT and RB filtering QC - ", sample_name)))
  dev.off()
  
  # Filter genes and cells
  sample_seurat_qc <- subset(sample_seurat, 
                             subset = nCount_RNA < 200 &
                               nFeature_RNA > 2500 &
                               percent.mt < 5 &
                               percent.rb < 60)
  
  # Scatter plots for QC validation
  scatter_plots <- FeatureScatter(sample_seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    FeatureScatter(sample_seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.rb") +
    FeatureScatter(sample_seurat_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm')
  
  # Save scatter plots
  pdf(paste0(snakemake@output[["scatter_plots"]], "_", sample_name, ".pdf"), width = 10, height = 6)
  print(scatter_plots)
  dev.off()
  
  # Normalize data
  sample_seurat_qc <- SCTransform(sample_seurat_qc, vars.to.regress = c("percent.mt", "percent.rb"))
  
  # Save the processed Seurat object for each sample
  saveRDS(sample_seurat_qc, file = paste0(snakemake@output[["seurat_qc"]], "_", sample_name, ".rds"))
}
