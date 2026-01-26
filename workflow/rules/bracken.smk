rule bracken:
    input:
        report = "results/{sample}/kraken/{sample}_report.txt"
    output:
        bracken_out = "results/{sample}/bracken/{sample}_species.txt"
    params:
        db = config["kraken_db"],
        read_len = config["bracken_read_len"],
        level = config["bracken_level"]
    threads: 10
    conda:
        "../envs/bracken.yaml"
    shell:
        """
        bracken -d {params.db} \
                -i {input.report} \
                -o {output.bracken_out} \
                -r {params.read_len} \
                -l {params.level}
        """
