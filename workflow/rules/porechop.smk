
rule porechop:
    input:
        "resources/reads/{sample}.fastq.gz"
    output:
        "results/{sample}/trimmed.fastq.gz"
    conda:
        "../envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} --threads 24"
