#!/usr/bin/env bash

################################################################################
# AMR Pipeline for Metagenomic Analysis
# Author       : Thibaut Armel Chérif GNIMADI
# Affiliation  : CERFIG
# Description  : Pipeline d'analyse métagénomique avec Kraken2/Bracken
# Version      : 2.0
# Date         : 2026-02-11
################################################################################

set -euo pipefail  # Arrêt si erreur, variables non définies, ou erreur dans pipe

# ==============================================================================
# CONFIGURATION
# ==============================================================================

CORES=22
CONFIG_FILE="config/config.yaml"
SNAKEFILE="workflow/Snakefile"
RESULTS_DIR="results"
SUMMARY_DIR="${RESULTS_DIR}/summary"

echo "===================================="
echo "  EXÉCUTION DU PIPELINE"
echo "===================================="
echo "Configuration : ${CONFIG_FILE}"
echo "Cœurs         : ${CORES}"
echo "Snakefile     : ${SNAKEFILE}"
echo ""

snakemake \
    --snakefile "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --cores "${CORES}" \
    --use-conda \
    --rerun-triggers mtime 

echo ""
echo "✓ Pipeline terminé avec succès"
echo ""

# Rapport HTML
echo "Génération du rapport HTML..."
snakemake \
    --snakefile "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --report "${SUMMARY_DIR}/workflow_report.html" \
    2>/dev/null || echo "⚠ Rapport HTML non généré (optionnel)"

# DAG (Directed Acyclic Graph)
echo "Génération du graphique du workflow..."
if command -v dot &> /dev/null; then
    snakemake \
        --snakefile "${SNAKEFILE}" \
        --configfile "${CONFIG_FILE}" \
        --dag | dot -Tpng > "${RESULTS_DIR}/workflow_dag.png"
    echo "✓ DAG sauvegardé: ${RESULTS_DIR}/workflow_dag.png"
else
    echo "⚠ Graphviz (dot) non installé, DAG non généré"
    echo "  Installation: conda install graphviz"
fi

echo ""
echo "===================================="
echo "  PIPELINE TERMINÉ"
echo "===================================="
echo "Résultats    : ${RESULTS_DIR}/"
echo "Résumés      : ${SUMMARY_DIR}/"
echo "===================================="
