rule minimap2_align:
    input:
        reads="results/{sample}/clean/{sample}_noh.fq.gz",
        ref="reference/{ref}.fasta"
    output:
        bam="results/{sample}/align/{sample}_{ref}.bam"
    params:
        threads=16
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} | samtools sort -@ {params.threads} -o {output.bam} && samtools index {output.bam}"
