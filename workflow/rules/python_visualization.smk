rule generate_arg_visualizations:
    input:
        arg_matrix_normalized = "results/r_analysis/arg_matrix_normalized.tsv",
        arg_matrix_counts = "results/r_analysis/arg_matrix_counts.tsv",
        arg_matrix_presence = "results/r_analysis/arg_matrix_presence.tsv",
        drug_class = "results/r_analysis/drug_class_abundance.tsv",
        mechanism = "results/r_analysis/mechanism_abundance.tsv",
        family = "results/r_analysis/family_abundance.tsv",
        arg_counts = "results/r_analysis/arg_counts.tsv"
    output:
        fig1 = "results/figures/1_drug_classes_distribution.png",
        fig2 = "results/figures/2_alpha_diversity.png",
        fig3 = "results/figures/3_heatmap_top30.png",
        fig4 = "results/figures/4_resistance_mechanisms.png",
        fig5 = "results/figures/5_rarefaction_curves.png",
        fig6 = "results/figures/6_pcoa_beta_diversity.png",
        fig7 = "results/figures/7_top20_args.png",
        colab_notebook = "results/figures/ARG_Visualization_Colab.ipynb"
    log:
        "logs/python_visualization/generate_figures.log"
    conda:
        "../envs/python_viz.yaml"
    threads: 1
    script:
        "../scripts/generate_arg_visualizations.py"
