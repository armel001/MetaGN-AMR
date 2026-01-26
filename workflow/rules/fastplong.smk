rule fastplong:
    input:
        "resources/reads/{sample}.fastq.gz"
    output:
        "results/{sample}/fp_trimmed.fastq.gz"
    conda:
        "../envs/fastplong.yaml"
    threads: 16
    log:
        "logs/fastplong/{sample}.log"
    shell:

        """
        fastplong \
            -i {input} \
            -o {output} \
            -q 10 \
            -t {threads} \
            &> {log}
        """

