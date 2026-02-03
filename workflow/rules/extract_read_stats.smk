rule extract_nanoplot_stats:
    input:
        stats = "results/{sample}/qc/{sample}_Nanoplot/NanoStats.txt"
    output:
        "results/{sample}/stats/{sample}_read_stats.txt"
    log:
        "logs/extract_stats/{sample}.log"
    conda:
        "../envs/python_filter.yaml"
    script:
        "../scripts/extract_nanoplot_stats.py"


rule aggregate_read_stats:
    input:
        expand("results/{sample}/stats/{sample}_read_stats.txt",
               sample=config["samples_id"])
    output:
        "results/stats/sequencing_stats.tsv"
    params:
        samples = config["samples_id"]
    log:
        "logs/aggregate_stats/aggregate.log"
    conda:
        "../envs/python_filter.yaml"
    script:
        "../scripts/aggregate_read_stats.py"
