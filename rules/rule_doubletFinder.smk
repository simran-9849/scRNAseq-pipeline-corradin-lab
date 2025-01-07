rule run_doubletfinder:
    input:
        seurat_obj = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/seurat_qc.RDS"
    
    output:
        sweep_seurat = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/sweep_seurat.RDS",
        doublet_rds = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/doublet.RDS",
        doublet_counts = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/doublet_counts.txt",
        filtered_seurat = "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/results/doublet_filtered_seurat.RDS"
    
    script:
        "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/scripts/df_trial.R"