configfile: "config.yaml"

rule all:
    input:
        config["output"]["rds"]

rule create_merged_seurat_object:
    input:
        inFolders = config["samples"].values()
    output:
        rds = config["output"]["rds"]
    params:
        sampleNames = list(config["samples"].keys())
    script:
      "/lab/corradin_data/FOR_SIMRAN/scRNA_seq/scripts/seurat_cts.R"