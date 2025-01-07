library(Seurat)
library(patchwork)
library(ggplot2)

#Read seurat object
seurat <- readRDS(snakemake@input[["seurat"]])
seurat <- readRDS("/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/merged_seurat.RDS")

# Process mitochondrial and ribosomal genes
#mito_genes <- grep(pattern = "^MT-", x = rownames(x = seurat), value = TRUE, ignore.case = TRUE)
#ribo_genes <- grep(pattern = "^RP[SL]", x = rownames(x = seurat), value = TRUE, ignore.case = TRUE)
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
seurat <- PercentageFeatureSet(seurat, pattern = "^RP[SL]", col.name = "percent.rb")

#Visualize using violin plots
p1 <- VlnPlot(seurat, features = "nFeature_RNA")
p2 <- VlnPlot(seurat, features = "nCount_RNA")
p3 <- VlnPlot(seurat, features = "percent.mt")
p4 <- VlnPlot(seurat, features = "percent.rb")

# Combine plots using patchwork
qc_violin <- (p1 | p2) / (p3 | p4)

# Save plot
pdf("/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/plots/violin_QC_123.pdf", width = 12, height = 10)
print(qc_violin + plot_annotation(title = "MT and RB filtering QC"))
dev.off()

#Filter genes and cells 
seurat_qc <- subset(seurat, 
               subset = nCount_RNA < 200 &
                 nFeature_RNA > 2500 &
                 percent.mt < 5 &
                 percent.rb < 60)

#Scatter plots for QC validation
scatter_plots <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "percent.rb") +
  FeatureScatter(seurat_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

pdf("/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/plots/scatterPlot.pdf", width = 10, height = 6)
print(scatter_plots)
dev.off()

#Normalize data
seurat_qc <- SCTransform(seurat_qc, vars.to.regress = c("percent.mt", "percent.rb"))
#seurat_qc <- FindVariableFeatures(seurat_qc, selection.method = "vst", nfeatures = 2000)
#seurat_qc <- ScaleData(seurat_qc, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.rb"))

#Save plots
ggsave(filename = snakemake@output[["qc_violin"]], plot = qc_violin, width = 10, height = 6)
ggsave(filename = snakemake@output[["scatter_plots"]], plot = scatter_plots, width = 12, height = 8)
dev.off()

# Save the Seurat object
saveRDS(seurat_qc, file = snakemake@output[["seurat_qc"]])
#saveRDS(seurat_qc, file = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/seurat_qc.RDS")
