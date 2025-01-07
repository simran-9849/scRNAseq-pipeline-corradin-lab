library(DoubletFinder)
library(tidyverse)
library(ggplot2)
library(Seurat)

# Load the Seurat object
filtered_seurat <- seurat_qc

# Split the Seurat object by sample identifier
# Replace "origin.ident" with your actual sample identifier column name in the metadata
seurat_list <- SplitObject(filtered_seurat, split.by = "orig.ident")

# Function to process each sample
process_sample <- function(sample) {
  # Run paramSweep
  sweep_seurat <- paramSweep(sample, PCs = 1:30, sct = FALSE)
  summary_sweep_seurat <- summarizeSweep(sweep_seurat, GT = FALSE)
  bcmvn_seurat <- find.pK(summary_sweep_seurat)
  
  # Plot pK identification
  pk_plot <- ggplot(bcmvn_seurat, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line() +
    ggtitle(paste("pK Identification -", sample$orig.ident[1]))
  
  # Determine optimal pK
  pK <- bcmvn_seurat %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  # Estimate homotypic doublet proportion
  annotations <- sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  nExp_poi <- round(0.069 * nrow(sample@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder 
  doublet_finder <- doubletFinder(sample, 
                                  PCs = 1:30, 
                                  pN = 0.25, 
                                  pK = pK, 
                                  nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = FALSE)
  
  # Get the DoubletFinder column name
  df_col <- grep("DF.classifications", colnames(doublet_finder@meta.data), value = TRUE)
  
  # Visualize doublets
  doublet_plot <- DimPlot(doublet_finder, reduction = 'umap', group.by = df_col) +
    ggtitle(paste("Doublet Classification -", sample$orig.ident[1]))
  
  # Number of singlets and doublets
  doublet_table <- table(doublet_finder@meta.data[[df_col]])
  
  return(list(doublet_finder = doublet_finder, 
              pk_plot = pk_plot, 
              doublet_plot = doublet_plot, 
              doublet_table = doublet_table))
}

# Process all samples
processed_samples <- lapply(seurat_list, process_sample)

# Combine results into a single Seurat object
combined_seurat <- merge(processed_samples[[1]]$doublet_finder, 
                         y = lapply(processed_samples[-1], function(x) x$doublet_finder))

# Save the combined sweep results
saveRDS(lapply(processed_samples, function(x) x$doublet_finder), snakemake@output[["sweep_seurat"]])

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

