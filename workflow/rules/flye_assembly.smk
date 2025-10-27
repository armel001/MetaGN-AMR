rule flye_assembly:
    input:
        "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        "results/{sample}/assembly.fasta"
    conda:
        "../envs/flye.yaml"
    shell:
        "flye --nano-raw {input} --out-dir results/{wildcards.sample}/flye --threads 24 && cp results/{wildcards.sample}/flye/assembly.fasta {output}"
