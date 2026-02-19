rule normalize_and_matrices:
    input:
        rgi_aggregated = "results/summary/rgi_all_samples.tsv",  # ← Corrigé
        sequencing_stats = "results/stats/sequencing_stats.tsv"
    output:
        arg_counts = "results/r_analysis/arg_counts.tsv",
        arg_relative = "results/r_analysis/arg_relative.tsv",
        arg_normalized = "results/r_analysis/arg_normalized.tsv",
        arg_presence = "results/r_analysis/arg_presence.tsv",
        arg_matrix_counts = "results/r_analysis/arg_matrix_counts.tsv",
        arg_matrix_relative = "results/r_analysis/arg_matrix_relative.tsv",
        arg_matrix_normalized = "results/r_analysis/arg_matrix_normalized.tsv",
        arg_matrix_presence = "results/r_analysis/arg_matrix_presence.tsv",
        drug_class_abundance = "results/r_analysis/drug_class_abundance.tsv",
        mechanism_abundance = "results/r_analysis/mechanism_abundance.tsv",
        family_abundance = "results/r_analysis/family_abundance.tsv"
    log:
        "logs/r_analysis/normalize_and_matrices.log"
    conda:
        "../envs/r_analysis.yaml"
    threads: 1
    script:
        "../scripts/normalize_and_matrices.R"
