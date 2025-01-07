configfile: "config.yaml"

rule qc_and_process:
    input:
        seurat = config["input"]["seurat"]
    output:
        qc_violin = config["output"]["qc_violin"],
        scatter_plots = config["output"]["scatter_plots"],
        seurat_qc = config["output"]["seurat_qc"]
    script:
        "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/scripts/xyz.R"

rule all:
    input:
        rules.qc_and_process.output
