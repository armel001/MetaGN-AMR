# ============================================================================
# Règles pour l'analyse taxonomique multi-niveaux
# ============================================================================

# Niveaux taxonomiques
TAXONOMIC_LEVELS = ["S", "G", "F", "P"]

# ============================================================================
# Combinaison des échantillons
# ============================================================================

rule combine_bracken_samples:
    input:
        expand("results/{sample}/bracken/{sample}_{{level}}.txt", 
               sample=config["samples_id"])
    output:
        abundance_matrix = "results/taxonomy/abundance_matrix_{level}.tsv",
        metadata = "results/taxonomy/sample_metadata_{level}.tsv"
    params:
        samples = config["samples_id"]
    wildcard_constraints:
        level = "S|G|F|P"  # ← AJOUT : Contrainte stricte sur le niveau
    conda:
        "../envs/python_filter.yaml"
    log:
        "logs/combine/combine_{level}.log"
    script:
        "../scripts/combine_bracken.py"

# ============================================================================
# Filtrage
# ============================================================================

rule filter_low_abundance:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}.tsv"
    output:
        filtered = "results/taxonomy/abundance_matrix_{level}_filtered.tsv",
        removed = "results/taxonomy/removed_taxa_{level}.txt"
    params:
        min_reads = config.get("min_reads_per_taxon", 100),
        min_samples = config.get("min_samples_present", 2),
        min_abundance = config.get("min_relative_abundance", 0.001)
    wildcard_constraints:
        level = "S|G|F|P"  # ← AJOUT : Contrainte stricte sur le niveau
    conda:
        "../envs/python_filter.yaml"
    log:
        "logs/filter/filter_{level}.log"
    script:
        "../scripts/filter_abundance.py"

# ============================================================================
# Normalisation
# ============================================================================

rule normalize_abundance:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}_filtered.tsv"
    output:
        relative = "results/taxonomy/abundance_matrix_{level}_relative.tsv",
        clr = "results/taxonomy/abundance_matrix_{level}_clr.tsv"
    wildcard_constraints:
        level = "S|G|F|P"  # ← AJOUT : Contrainte stricte sur le niveau
    conda:
        "../envs/python_filter.yaml"
    log:
        "logs/normalize/normalize_{level}.log"
    script:
        "../scripts/normalize_abundance.py"

# ============================================================================
# Diversité alpha
# ============================================================================

rule calculate_alpha_diversity:
    input:
        matrix = "results/taxonomy/abundance_matrix_S.tsv"
    output:
        diversity = "results/taxonomy/alpha_diversity.tsv"
    conda:
        "../envs/python_filter.yaml"
    log:
        "logs/diversity/alpha_diversity.log"
    script:
        "../scripts/alpha_diversity.py"

# ============================================================================
# Résumé
# ============================================================================

rule summary_statistics:
    input:
        kraken_reports = expand("results/{sample}/kraken/{sample}_report.txt", 
                               sample=config["samples_id"]),
        diversity = "results/taxonomy/alpha_diversity.tsv",
        matrices = expand("results/taxonomy/abundance_matrix_{level}_filtered.tsv",
                         level=TAXONOMIC_LEVELS)
    output:
        summary = "results/summary/taxonomy_summary.txt"
    conda:
        "../envs/python_filter.yaml"
    log:
        "logs/summary/taxonomy_summary.log"
    script:
        "../scripts/taxo_summary.py"
