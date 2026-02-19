rule plot_taxonomy_figures:
    input:
        genus_matrix   = fig_genus_matrix,    # ← défini dans Snakefile selon stratégie
        species_matrix = fig_species_matrix,  # ← défini dans Snakefile selon stratégie
        alpha_div      = fig_alpha_div        # ← défini dans Snakefile selon stratégie
    output:
        composition_genus   = "results/figures/taxonomy/composition_genus.png",
        composition_species = "results/figures/taxonomy/composition_species.png",
        alpha_diversity     = "results/figures/taxonomy/alpha_diversity.png",
        heatmap_species     = "results/figures/taxonomy/heatmap_species.png",
        pathogens           = "results/figures/taxonomy/pathogens.png"
    params:
        top_genera  = config.get("fig_top_genera",    15),
        top_species = config.get("fig_top_species",   20),
        top_heatmap = config.get("fig_heatmap_top_n", 30),
        dpi         = config.get("fig_dpi", 300),
        pathogens   = config.get("priority_pathogens", [
            "Escherichia coli",
            "Klebsiella pneumoniae",
            "Pseudomonas aeruginosa",
            "Clostridioides difficile",
            "Clostridium perfringens",
            "Acinetobacter baumannii",
            "Enterococcus faecium"
        ])
    conda:
        "../envs/python_viz.yaml"
    log:
        "logs/figures/taxonomy_figures.log"
    script:
        "../scripts/plot_taxonomy.py"
