library(Seurat)

# Get input folders and sample names
input_folders <- snakemake@input[["inFolders"]]
sample_names <- snakemake@params[["sampleNames"]]

# Read files and create Seurat objects
seurat_objects <- lapply(seq_along(input_folders), function(i) {
  data <- Read10X(data.dir = input_folders[i])
  CreateSeuratObject(counts = data, project = sample_names[i])
})

# Merge Seurat objects
merged_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1], 
                       add.cell.ids = sample_names, project = "MergedProject")

# Set orig.ident
merged_seurat$orig.ident <- factor(merged_seurat$orig.ident, levels = sample_names, labels = sample_names)

# Save merged Seurat object
saveRDS(merged_seurat, file = snakemake@output[["rds"]])