library(scater)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
sce_completed <- readRDS("/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/sce.RDS")
#read sc data and create sce object
mtx <- readMM(snakemake@input[["matrix"]])
features <- read.delim(snakemake@input[["features"]], header = FALSE, 
                       stringsAsFactors = FALSE)
barcodes <- read.delim(snakemake@input[["barcodes"]], header = FALSE, 
                       stringsAsFactors = FALSE)

rownames(mtx) <- features$V1
colnames(mtx) <- barcodes$V1

sce <- SingleCellExperiment(
  assays = list(counts = mtx)
)

rowData(sce)$gene_id <- features$V1
rowData(sce)$gene_name <- features$V2

# get mitochondrial genes
is.mito <- grepl("^MT-", rowData(sce)$gene_name, ignore.case = TRUE)

# Identify spike-ins by finding rows that match the pattern (e.g., starting with "ERCC")
is_spike <- grepl(snakemake@params[[1]], rownames(sce))

# Add a new column in `rowData` to store spike-in status
rowData(sce)$is_spike <- is_spike

# calculate metrics
sce <- addPerCellQCMetrics(sce, subsets = list(Spike=is_spike, Mt=is.mito))

#Convert sce object to df to perform library size calculations
sce_df <- as.data.frame(colData(sce))
sce_df$total_features_by_counts <- colSums(assay(sce, "counts") > 0)
sce_df$detection_rate <- sce_df$total_features_by_counts / nrow(sce)
sce <- DataFrame(sce_df)

saveRDS(sce, file=snakemake@output[[1]])

