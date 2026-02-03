rule nanoplot:
    input:
        fastq = "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        stats = "results/{sample}/qc/{sample}_Nanoplot/NanoStats.txt",
        report = "results/{sample}/qc/{sample}_Nanoplot/NanoPlot-report.html"
    params:
        outdir = lambda wildcards: f"results/{wildcards.sample}/qc/{wildcards.sample}_Nanoplot"
    threads: 20
    resources:
        mem_mb = 16000
    conda:
        "../envs/nanoplot.yaml"
    log:
        "logs/nanoplot/{sample}.log"
    benchmark:
        "benchmarks/nanoplot/{sample}.txt"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        
        export TMPDIR={params.outdir}/tmp
        export MPLBACKEND=Agg
        export MPLCONFIGDIR={params.outdir}/.mplconfig
        
        NanoPlot \
          --fastq {input.fastq} \
          -o {params.outdir} \
          --threads {threads} \
          --huge \
          -f png \
          --verbose \
        > {log} 2>&1 || {{
            echo "Number of reads: 0" > {output.stats}
            echo "<html><body>NanoPlot failed</body></html>" > {output.report}
        }}
        
        rm -rf {params.outdir}/tmp {params.outdir}/.mplconfig
        """
