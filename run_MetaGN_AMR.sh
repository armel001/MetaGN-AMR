#!/usr/bin/env bash

################################################################################
# AMR Pipeline for Metagenomic Analysis
# Author       : Thibaut Armel Chérif GNIMADI
# Affiliation  : CERFIG
# Description  : Pipeline d'analyse métagénomique avec Kraken2/Bracken
# Version      : 2.1
# Date         : 2026-02-11
################################################################################

set -euo pipefail

# ==============================================================================
# CONFIGURATION PAR DÉFAUT
# ==============================================================================

CORES=${CORES:-22}
CONFIG_FILE=${CONFIG_FILE:-"config/config.yaml"}
SNAKEFILE=${SNAKEFILE:-"workflow/Snakefile"}
RESULTS_DIR=${RESULTS_DIR:-"results"}
SUMMARY_DIR="${RESULTS_DIR}/summary"
LOGS_DIR="${RESULTS_DIR}/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/pipeline_${TIMESTAMP}.log"

# Options par défaut
DRY_RUN=false
VERBOSE=false
UNLOCK=false
CLUSTER_MODE=false
GENERATE_REPORTS=true

# ==============================================================================
# FONCTIONS UTILITAIRES
# ==============================================================================

# Affichage avec couleurs
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "${LOG_FILE}"
}

log_success() {
    echo -e "${GREEN}[✓]${NC} $1" | tee -a "${LOG_FILE}"
}

log_warning() {
    echo -e "${YELLOW}[⚠]${NC} $1" | tee -a "${LOG_FILE}"
}

log_error() {
    echo -e "${RED}[✗]${NC} $1" | tee -a "${LOG_FILE}"
}

# Afficher l'usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

AMR Pipeline for Metagenomic Analysis

OPTIONS:
    -h, --help              Afficher cette aide
    -n, --dry-run           Exécution à sec (voir ce qui serait fait)
    -c, --cores N           Nombre de cœurs (défaut: ${CORES})
    -v, --verbose           Mode verbeux
    -u, --unlock            Déverrouiller le répertoire de travail
    --no-reports            Ne pas générer les rapports finaux
    --cluster               Mode cluster (SLURM)
    --config FILE           Fichier de configuration (défaut: ${CONFIG_FILE})

EXEMPLES:
    $0                                  # Exécution standard
    $0 -n                               # Dry-run pour vérifier
    $0 -c 32 --verbose                  # 32 cœurs en mode verbeux
    $0 --unlock                         # Déverrouiller après crash
    $0 --cluster                        # Exécution sur cluster SLURM

EOF
    exit 0
}

# ==============================================================================
# PARSING DES ARGUMENTS
# ==============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        --no-reports)
            GENERATE_REPORTS=false
            shift
            ;;
        --cluster)
            CLUSTER_MODE=true
            shift
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        *)
            log_error "Option inconnue: $1"
            usage
            ;;
    esac
done

# ==============================================================================
# INITIALISATION
# ==============================================================================

# Créer les dossiers nécessaires
mkdir -p "${RESULTS_DIR}" "${SUMMARY_DIR}" "${LOGS_DIR}"

# Banner
cat << "EOF"
╔════════════════════════════════════════════════════════════════╗
║                                                                ║
║        AMR METAGENOMIC ANALYSIS PIPELINE - CERFIG             ║
║                                                                ║
╚════════════════════════════════════════════════════════════════╝
EOF

log_info "Démarrage du pipeline: $(date)"
log_info "Fichier de log: ${LOG_FILE}"
echo ""

# ==============================================================================
# VÉRIFICATIONS PRÉLIMINAIRES
# ==============================================================================

log_info "Vérifications préliminaires..."

# Vérifier Snakemake
if ! command -v snakemake &> /dev/null; then
    log_error "Snakemake n'est pas installé"
    log_info "Installation: conda install -c bioconda snakemake"
    exit 1
fi
SNAKEMAKE_VERSION=$(snakemake --version)
log_success "Snakemake ${SNAKEMAKE_VERSION} détecté"

# Vérifier le Snakefile
if [[ ! -f "${SNAKEFILE}" ]]; then
    log_error "Snakefile introuvable: ${SNAKEFILE}"
    exit 1
fi
log_success "Snakefile trouvé: ${SNAKEFILE}"

# Vérifier le fichier de config
if [[ ! -f "${CONFIG_FILE}" ]]; then
    log_error "Fichier de configuration introuvable: ${CONFIG_FILE}"
    exit 1
fi
log_success "Configuration trouvée: ${CONFIG_FILE}"

# Vérifier conda
if ! command -v conda &> /dev/null; then
    log_warning "Conda non détecté (environnements ne seront pas créés automatiquement)"
else
    log_success "Conda détecté"
fi

echo ""

# ==============================================================================
# DÉVERROUILLAGE (si demandé)
# ==============================================================================

if [[ "${UNLOCK}" == true ]]; then
    log_info "Déverrouillage du répertoire de travail..."
    snakemake \
        --snakefile "${SNAKEFILE}" \
        --configfile "${CONFIG_FILE}" \
        --unlock
    log_success "Répertoire déverrouillé"
    exit 0
fi

# ==============================================================================
# CONFIGURATION DE L'EXÉCUTION
# ==============================================================================

SNAKEMAKE_OPTS=(
    --snakefile "${SNAKEFILE}"
    --configfile "${CONFIG_FILE}"
    --cores "${CORES}"
    --use-conda
    --rerun-triggers mtime
    --printshellcmds
    --reason
    --keep-going  # Continuer malgré les échecs non critiques
)

if [[ "${VERBOSE}" == true ]]; then
    SNAKEMAKE_OPTS+=(--verbose)
fi

if [[ "${DRY_RUN}" == true ]]; then
    SNAKEMAKE_OPTS+=(--dry-run)
    log_info "MODE DRY-RUN: Aucune commande ne sera exécutée"
fi

if [[ "${CLUSTER_MODE}" == true ]]; then
    log_info "Mode cluster activé (SLURM)"
    SNAKEMAKE_OPTS+=(
        --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M --time=24:00:00"
        --jobs 10
        --latency-wait 60
    )
fi

# ==============================================================================
# AFFICHAGE DE LA CONFIGURATION
# ==============================================================================

log_info "Configuration:"
log_info "  Cœurs              : ${CORES}"
log_info "  Configuration      : ${CONFIG_FILE}"
log_info "  Snakefile          : ${SNAKEFILE}"
log_info "  Résultats          : ${RESULTS_DIR}"
log_info "  Mode dry-run       : ${DRY_RUN}"
log_info "  Mode cluster       : ${CLUSTER_MODE}"
log_info "  Générer rapports   : ${GENERATE_REPORTS}"
echo ""

# ==============================================================================
# EXÉCUTION DU PIPELINE
# ==============================================================================

log_info "Démarrage de l'exécution du pipeline..."
echo ""

START_TIME=$(date +%s)

snakemake "${SNAKEMAKE_OPTS[@]}" 2>&1 | tee -a "${LOG_FILE}"

PIPELINE_EXIT_CODE=${PIPESTATUS[0]}
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
if [[ ${PIPELINE_EXIT_CODE} -ne 0 ]]; then
    log_error "Le pipeline a échoué (code: ${PIPELINE_EXIT_CODE})"
    log_info "Temps écoulé: ${ELAPSED}s"
    log_info "Consultez le fichier de log: ${LOG_FILE}"
    exit ${PIPELINE_EXIT_CODE}
fi

log_success "Pipeline terminé avec succès"
log_info "Temps d'exécution: ${ELAPSED}s ($(date -ud "@${ELAPSED}" +"%H:%M:%S"))"
echo ""

# ==============================================================================
# GÉNÉRATION DES RAPPORTS (si activé et non dry-run)
# ==============================================================================

if [[ "${GENERATE_REPORTS}" == true ]] && [[ "${DRY_RUN}" == false ]]; then
    log_info "Génération des rapports..."
    
    # Rapport HTML
    log_info "  → Rapport HTML..."
    if snakemake \
        --snakefile "${SNAKEFILE}" \
        --configfile "${CONFIG_FILE}" \
        --report "${SUMMARY_DIR}/workflow_report.html" \
        >> "${LOG_FILE}" 2>&1; then
        log_success "    Rapport HTML: ${SUMMARY_DIR}/workflow_report.html"
    else
        log_warning "    Rapport HTML non généré (optionnel)"
    fi
    
    # DAG
    log_info "  → Graphique du workflow (DAG)..."
    if command -v dot &> /dev/null; then
        if snakemake \
            --snakefile "${SNAKEFILE}" \
            --configfile "${CONFIG_FILE}" \
            --dag | dot -Tpng > "${RESULTS_DIR}/workflow_dag.png" 2>> "${LOG_FILE}"; then
            log_success "    DAG: ${RESULTS_DIR}/workflow_dag.png"
        else
            log_warning "    DAG non généré"
        fi
    else
        log_warning "    Graphviz (dot) non installé"
        log_info "    Installation: conda install graphviz"
    fi
    
    # Rulegraph (graphe simplifié)
    log_info "  → Graphique des règles..."
    if command -v dot &> /dev/null; then
        if snakemake \
            --snakefile "${SNAKEFILE}" \
            --configfile "${CONFIG_FILE}" \
            --rulegraph | dot -Tpng > "${RESULTS_DIR}/rulegraph.png" 2>> "${LOG_FILE}"; then
            log_success "    Rulegraph: ${RESULTS_DIR}/rulegraph.png"
        fi
    fi
    
    echo ""
fi

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

cat << EOF
╔════════════════════════════════════════════════════════════════╗
║                    PIPELINE TERMINÉ                            ║
╚════════════════════════════════════════════════════════════════╝

📁 Résultats principaux:
   • Résultats          : ${RESULTS_DIR}/
   • Résumés            : ${SUMMARY_DIR}/
   • Logs               : ${LOGS_DIR}/
   • Log de cette exec. : ${LOG_FILE}

📊 Fichiers de sortie clés:
   • Diversité alpha    : ${RESULTS_DIR}/taxonomy/alpha_diversity*.tsv
   • Matrices filtrées  : ${RESULTS_DIR}/taxonomy/abundance_matrix_*_filtered*.tsv
   • Résumé taxonomique : ${RESULTS_DIR}/summary/taxonomy_summary*.txt

⏱  Temps d'exécution: ${ELAPSED}s
📅 Fin: $(date)

EOF

# ==============================================================================
# STATISTIQUES FINALES
# ==============================================================================

if [[ "${DRY_RUN}" == false ]]; then
    log_info "Statistiques:"
    
    # Compter les fichiers générés
    N_OUTPUTS=$(find "${RESULTS_DIR}" -type f 2>/dev/null | wc -l)
    log_info "  Fichiers générés: ${N_OUTPUTS}"
    
    # Taille totale
    TOTAL_SIZE=$(du -sh "${RESULTS_DIR}" 2>/dev/null | cut -f1)
    log_info "  Taille totale: ${TOTAL_SIZE}"
fi

log_success "Analyse terminée ! ✨"
