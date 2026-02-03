rule bracken_multilevel:
    input:
        report = "results/{sample}/kraken/{sample}_report.txt"
    output:
        bracken_out = "results/{sample}/bracken/{sample}_{level}.txt",
        bracken_report = "results/{sample}/bracken/{sample}_{level}_report.txt"
    params:
        db = config["kraken_db"],
        read_len = config["bracken_read_len"],
        threshold = config.get("bracken_threshold", 10)
    threads: 4
    conda:
        "../envs/bracken.yaml"
    log:
        "logs/braken/{sample}_{level}.log"
    shell:
        """
        bracken -d {params.db} \
                -i {input.report} \
                -o {output.bracken_out} \
                -w {output.bracken_report} \
                -r {params.read_len} \
                -l {wildcards.level} \
                -t {params.threshold} \
                2> {log}
        """
