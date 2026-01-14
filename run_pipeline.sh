#!/bin/zsh
set -e

# --- DEFAULTS ---
TOP_N=10; TERM_MAX=300; PVAL=0.01; FOLD=5

# --- ARGUMENT PARSER ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --top_n) TOP_N="$2"; shift ;;
        --term_max) TERM_MAX="$2"; shift ;;
        --pval) PVAL="$2"; shift ;;
        --fold) FOLD="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# --- PATHS ---
# Using current environment binaries
PYTHON=$(which python3)
RSCRIPT=$(which Rscript)

DATA="data/Ribosome_Associated_proteome.tsv"
RES_DIR="results"
mkdir -p $RES_DIR

echo "------------------------------------------------------------"
echo ">> RUNNING PROJECT PIPELINE"
echo ">> Results target: $RES_DIR"
echo "------------------------------------------------------------"

# STEP 1: Volcano (Internal output path changed to results/ via the script)
echo ">> Running Volcano Plots..."
$RSCRIPT src/Curved_volcanoplots_Figure3.R

# STEP 2 & 3: GO Enrichment
mkdir -p $RES_DIR/GO_UP $RES_DIR/GO_DOWN
$PYTHON src/simplifiedGO.py "$DATA" --filter "_up" --outdir "$RES_DIR/GO_UP"
$PYTHON src/simplifiedGO.py "$DATA" --filter "_down" --outdir "$RES_DIR/GO_DOWN"

# STEP 4: DotPlots
echo ">> Running GO DotPlots..."
$RSCRIPT src/SimplifiedGO_ClusterBased.R "$RES_DIR/GO_UP" "$RES_DIR/Figure3_GO_UP.pdf" "UP Clusters" $TOP_N $TERM_MAX $PVAL $FOLD
$RSCRIPT src/SimplifiedGO_ClusterBased.R "$RES_DIR/GO_DOWN" "$RES_DIR/Figure3_GO_DOWN.pdf" "DOWN Clusters" $TOP_N $TERM_MAX $PVAL $FOLD

echo ">> PIPELINE COMPLETE. See results/ directory."
