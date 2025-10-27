rule plot_ARGs_by_gene:
    input:
        summary="results/{sample}/abricate/{sample}_ARG_summary.tsv"
    output:
        plot="results/{sample}/plots/{sample}_ARGs_by_gene.png"
    params:
        sample="{sample}"
    conda:
        "../envs/python_plot.yaml"
    shell:
        """
        python /data/armel/git/MetagenAMR-main/workflow/scripts/plot_ARGs_by_gene.py {input.summary} {output.plot} {params.sample}
        """
