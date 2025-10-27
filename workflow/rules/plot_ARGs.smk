
rule plot_ARGs:
    input:
        abricate="results/{sample}/abricate/{sample}_ARG_summary.tsv"
    output:
        plot="results/{sample}/plots/{sample}_resistance_classes.png"
    conda:
        "../envs/python_plot.yaml"
    script:
        "../scripts/plot_resistance_classes.py"
