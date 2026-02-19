rule flye_assembly:
    input:
        "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        "results/{sample}/assembly.fasta"
    conda:
        "../envs/flye.yaml"
    log:
        "logs/flye_assembly/{sample}.log"
    benchmark:
        "benchmarks/flye_assembly/{sample}.txt"
    threads: 40
    resources:
        mem_mb = 128000
    shadow:
        "minimal"
    shell:
        """
        flye \
            --nano-raw {input} \
            --out-dir flye_outdir \
            --threads {threads} \
            --meta \
            --iterations 3 \
            2> {log} \
            && cp flye_outdir/assembly.fasta {output}
        """
