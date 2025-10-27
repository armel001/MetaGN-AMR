
rule abricate:
    input:
        assembly="results/{sample}/assembly.fasta"
    output:
        "results/{sample}/abricate/{sample}_ARG_summary.tsv"
    conda:
        "../envs/abricate.yaml"
    shell:
        "abricate --db resfinder {input.assembly} > {output}"