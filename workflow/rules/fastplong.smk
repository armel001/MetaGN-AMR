rule fastplong:
    input:
        "resources/reads/{sample}.fastq.gz"
    output:
        fq="results/{sample}/fastplong/fp_trimmed.fastq.gz",
        html="results/{sample}/fastplong/fastplong.html",
        json="results/{sample}/fastplong/fastplong.json"
    conda:
        "../envs/fastplong.yaml"
    threads: 16
    log:
        "logs/fastplong/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.fq})
        mkdir -p $(dirname {log})
        fastplong \
            -i {input} \
            -o {output.fq} \
            -q 10 \
            -t {threads} \
            -h {output.html} \
            -j {output.json} \
            &> {log}
        """
