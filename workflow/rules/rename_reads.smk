rule rename_reads:
    input:
        "resources/reads/{sample}.fastq.gz"
    output:
        "resources/reads/{sample}.renamed.fastq.gz"
    shell:
        r"""
        zcat {input} \
        | awk 'NR%4==1 {{printf("@{wildcards.sample}_%d\n", ++i)}} NR%4!=1' \
        | gzip > {output}
        """

