rule aggregate_rgi_samples:
    input:
        selected = expand("results/{sample}/rgi_preprocessing/{sample}_selected.txt",
                         sample=config["samples_id"])
    output:
        aggregated = "results/summary/rgi_all_samples.tsv",
        summary = "results/summary/rgi_aggr_summary.txt"
    params:
        samples = config["samples_id"]
    log:
        "logs/rgi_aggregate/aggregate.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 8
    resources:
        mem_mb = 80000 
    script:
        "../scripts/aggregate_rgi_samples.py"
