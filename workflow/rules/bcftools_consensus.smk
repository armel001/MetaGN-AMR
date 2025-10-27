
rule bcftools_consensus:
    input:
        ref="reference/{ref}.fasta",
        bam="results/{sample}/align/{sample}_{ref}.bam"
    output:
        vcf="results/{sample}/consensus/{sample}_{ref}.vcf.gz",
        consensus="results/{sample}/consensus/{sample}_{ref}_consensus.fasta"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Oz -o {output.vcf} && bcftools index {output.vcf} && cat {input.ref} | bcftools consensus {output.vcf} > {output.consensus}"
