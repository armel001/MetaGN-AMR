rule host_depletion:
    input:
        reads="results/{sample}/fp_trimmed.fastq.gz",
        host_ref="resources/references/GRCh38_no_alt.mmi"
    output:
        clean_reads="results/{sample}/clean/{sample}_noh.fq.gz"
    threads: 24
    conda:
        "../envs/host_depletion.yaml"
    shell:
        """
        minimap2 -ax map-ont {input.host_ref} {input.reads} | \
        samtools view -b -f 4 -@ {threads} | \
        samtools fastq -@ {threads} - | gzip > {output.clean_reads}
        """
