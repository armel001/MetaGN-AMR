rule plot_ARGs_all:
    input:
        summaries=expand("results/{sample}/abricate/{sample}_ARG_summary.tsv", sample=config["samples_id"])
    output:
        plot="results/summary/all_ARGs.png",
        table="results/summary/all_ARGs_counts.tsv"
    conda:
        "../envs/python_plot.yaml"
    shell:
        """
        mkdir -p results/summary
        python /data/armel/git/MetagenAMR-main/workflow/scripts/plot_ARGs_all.py {output.plot} {output.table} {input.summaries}
        """
