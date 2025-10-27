rule plot_ARGs_all_R:
    input:
        summaries=expand("results/{sample}/abricate/{sample}_ARG_summary.tsv", sample=config["samples_id"])
    output:
        plot="results/summary/r_all_ARGs.png",
        table="results/summary/r_all_ARGs_counts.tsv"
    conda:
        "../envs/r_env.yaml"
    shell:
        """
        mkdir -p results/summary
        Rscript /data/armel/git/MetagenAMR-main/workflow/scripts/r_plot_ARGs_all.R {output.plot} {output.table} {input.summaries}
        """
