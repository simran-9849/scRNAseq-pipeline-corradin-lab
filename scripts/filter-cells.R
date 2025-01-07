#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

library(scater)

seurat_qc <- readRDS(snakemake@input[[1]])
seurat_qc <- as.SingleCellExperiment(seurat_obj)

# drop cells with too few counts or expressed features
libsize.drop <- isOutlier(seurat_qc$nCount_RNA, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(seurat_qc$nFeature_RNA, nmads=3, type="lower", log=TRUE)

# drop cells with too high proportion of mito genes or spike-ins expressed
mito.drop <- isOutlier(seurat_qc$percent.mt, nmads=3, type="higher")
#spike.drop <- isOutlier(sce$, nmads=3, type="higher")

# write stats
stats <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
     ByMito=sum(mito.drop), Remaining=ncol(sce))
write.table(stats, file=snakemake@output[["stats"]], sep='\t',
            row.names=FALSE, quote=FALSE)

# filter sce
seurat_qc <- seurat_qc[,!(libsize.drop | feature.drop | mito.drop)]

saveRDS(seurat_qc, file=snakemake@output[["rds"]])

