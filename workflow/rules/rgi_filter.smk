rule filter_rgi_quality:
    input:
        rgi = "results/{sample}/rgi/{sample}.txt"
    output:
        filtered = "results/{sample}/rgi_preprocessing/{sample}_filtered.txt"
    params:
        min_identity = config["rgi"]["min_identity"],
        min_coverage = config["rgi"]["min_coverage"],
        model_types = config["rgi"]["model_types"]
    log:
        "logs/rgi_filter/{sample}.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 1
    resources:
        mem_mb = 2000
    script:
        "../scripts/filter_rgi_results.py"
