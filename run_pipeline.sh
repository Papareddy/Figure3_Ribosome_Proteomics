#!/bin/zsh
set -e

# Get the directory where the script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# --- DEFAULTS ---
TOP_N=10; TERM_MAX=300; PVAL=0.01; FOLD=5

# --- ARGUMENT PARSER ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            echo "Usage: $(basename "$0") [--top_n <int>] [--term_max <int>] [--pval <float>] [--fold <float>]"
            echo ""
            echo "This pipeline analyzes ribosome-associated proteome data."
            echo ""
            echo "Options:"
            echo "  --top_n <int>   : Number of top GO terms to display per cluster in dotplots (default: 10)."
            echo "  --term_max <int>: Maximum size of GO terms to consider (default: 300)."
            echo "  --pval <float>  : P-value threshold for GO terms (default: 0.01)."
            echo "  --fold <float>  : Fold enrichment threshold for GO terms (default: 5)."
            echo "  -h, --help      : Display this help message and exit."
            exit 0
            ;;
        --top_n) TOP_N="$2"; shift ;;
        --term_max) TERM_MAX="$2"; shift ;;
        --pval) PVAL="$2"; shift ;;
        --fold) FOLD="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# --- PATHS ---
PYTHON=$(which python)
RSCRIPT=$(which Rscript)

DATA="${SCRIPT_DIR}/data/Ribosome_Associated_proteome.tsv"
RES_DIR="${SCRIPT_DIR}/results"
mkdir -p $RES_DIR

echo "------------------------------------------------------------"
echo ">> RUNNING PROJECT PIPELINE"
echo ">> Results target: $RES_DIR"
echo "------------------------------------------------------------"

# STEP 1: Volcano (Internal output path changed to results/ via the script)
echo ">> Running Volcano Plots..."
$RSCRIPT "${SCRIPT_DIR}/src/Curved_volcanoplots_Figure3.R" "$DATA"

# STEP 2 & 3: GO Enrichment
mkdir -p $RES_DIR/GO_UP $RES_DIR/GO_DOWN
$PYTHON "${SCRIPT_DIR}/src/GOterms_by_gprofiler.py" "$DATA" --filter "_up" --outdir "$RES_DIR/GO_UP"
$PYTHON "${SCRIPT_DIR}/src/GOterms_by_gprofiler.py" "$DATA" --filter "_down" --outdir "$RES_DIR/GO_DOWN"

# STEP 4: DotPlots
echo ">> Running GO DotPlots..."
$RSCRIPT "${SCRIPT_DIR}/src/plot_GOterm_clusters.R" "$RES_DIR/GO_UP" "$RES_DIR/Figure3_GO_UP.pdf" "UP Clusters" $TOP_N $TERM_MAX $PVAL $FOLD
$RSCRIPT "${SCRIPT_DIR}/src/plot_GOterm_clusters.R" "$RES_DIR/GO_DOWN" "$RES_DIR/Figure3_GO_DOWN.pdf" "DOWN Clusters" $TOP_N $TERM_MAX $PVAL $FOLD

echo ">> PIPELINE COMPLETE. See results/ directory."