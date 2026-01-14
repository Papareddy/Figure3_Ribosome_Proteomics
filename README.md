# Figure 3: Ribosome-Associated Proteome Analysis Pipeline

This repository contains the scripts and data to reproduce Figure 3, which involves analyzing ribosome-associated proteome data to generate curved volcano plots, perform GO enrichment analysis, and visualize the results as dotplots.

### Directory Structure

```
.
├── data/
│   └── Ribosome_Associated_proteome.tsv    # Raw proteome data
├── results/                                # Output directory for all generated figures and enrichment results
├── src/
│   ├── Curved_volcanoplots_Figure3.R       # R script for generating curved volcano plots
│   ├── SimplifiedGO_ClusterBased.R         # R script for generating GO enrichment dotplots
│   └── simplifiedGO.py                     # Python script for performing GO enrichment analysis
├── .vscode/                                # VS Code settings
├── Heatmap_figure3.R                       # Supplementary R script for heatmap generation (not part of main Figure 3 pipeline)
├── Readexcels_and_Barfplot_figure3.R       # Supplementary R script for barplot generation (not part of main Figure 3 pipeline)
├── RibosomeAssociatedProteome.yml          # Conda/Mamba environment definition
├── run_pipeline.sh                         # Main script to execute the entire analysis pipeline
└── README.md                               # This file
```

### Environment Setup

The analysis relies on a specific Conda/Mamba environment. You can set it up using the provided `RibosomeAssociatedProteome.yml` file:

```bash
mamba env create -f RibosomeAssociatedProteome.yml
conda activate RibosomeAssociatedProteome
```

### Pipeline Overview

The main analysis pipeline is orchestrated by the `run_pipeline.sh` script, which performs the following steps:

1.  **Volcano Plot Generation:** Creates curved volcano plots highlighting significant protein changes.
2.  **GO Enrichment Analysis (UP and DOWN):** Identifies over-represented Gene Ontology terms for "up-regulated" and "down-regulated" protein clusters.
3.  **GO Enrichment Dotplot Visualization:** Generates dotplots to visualize the most significant GO terms.

### Usage

To run the entire pipeline, ensure your `RibosomeAssociatedProteome` environment is activated, then execute the `run_pipeline.sh` script:

```bash
# Activate the environment
conda activate RibosomeAssociatedProteome

# Run the pipeline
./run_pipeline.sh
```

You can customize certain parameters of the pipeline using command-line arguments for `run_pipeline.sh`:

```bash
./run_pipeline.sh --top_n 10 --term_max 300 --pval 0.01 --fold 5
```
however the above parameters were used in the publication

**Available arguments for `run_pipeline.sh`:**
*   `--top_n <int>`: Number of top GO terms to display per cluster in dotplots (default: 10).
*   `--term_max <int>`: Maximum size of GO terms to consider (default: 300).
*   `--pval <float>`: P-value threshold for GO terms (default: 0.01).
*   `--fold <float>`: Fold enrichment threshold for GO terms (default: 5).

All generated output files, including PDF figures and GO enrichment TSV files, will be saved in the `results/` directory.

### Script Descriptions

#### `src/Curved_volcanoplots_Figure3.R`

*   **Purpose:** Generates a multi-panel PDF (`results/Figure3_Proteome_Volcanoes.pdf`) containing curved volcano plots for different ribosome fractions. It uses a custom hyperbolic threshold to define significance and colors points based on their significance (up/down/NS).
*   **Input Data:** Reads `data/Ribosome_Associated_proteome.tsv` (hardcoded path).
*   **Key Features:** Custom hyperbolic significance threshold, highlighting of specific genes.
*   **Dependencies:** `dplyr`.

#### `src/simplifiedGO.py`

*   **Purpose:** Performs Gene Ontology (GO) enrichment analysis on the proteome data. It takes the main data file, filters it based on cluster suffixes (`_up` or `_down`), and then uses g:Profiler to find enriched GO terms for each cluster.
*   **Input Data:** Expects a path to the main proteome TSV file (e.g., `data/Ribosome_Associated_proteome.tsv`).
*   **Output:** Generates a `cluster_sizes.tsv` file, `cluster_sizes.png` (barplot), `go_terms_all.tsv`, `go_terms_bp.tsv`, `go_terms_mf.tsv`, `go_terms_cc.tsv`, and category-specific heatmaps (`heatmap_bp.png`, `heatmap_mf.png`, `heatmap_cc.png`) within the specified output directory (e.g., `results/GO_UP` or `results/GO_DOWN`).
*   **Arguments:**
    *   `<input_file>`: Path to the input proteome TSV file.
    *   `--outdir <path>`: Output directory for GO enrichment results (default: `GO_Results`).
    *   `--filter <suffix>`: Suffix to filter clusters (e.g., `_up` or `_down`, default: `_up`).
*   **Organism:** `athaliana` (Arabidopsis thaliana).
*   **Dependencies:** `pandas`, `numpy`, `matplotlib`, `seaborn`, `gprofiler-official`.

#### `src/SimplifiedGO_ClusterBased.R`

*   **Purpose:** Visualizes the results of the GO enrichment analysis as multi-page PDF dotplots. It filters GO terms based on p-value, fold enrichment, and term size, and displays the top N terms for each cluster across different GO categories.
*   **Input Data:** Reads `go_terms_bp.tsv`, `go_terms_mf.tsv`, and `go_terms_cc.tsv` from a specified results directory (e.g., `results/GO_UP`).
*   **Output:** A multi-page PDF file containing dotplots (e.g., `results/Figure3_GO_UP.pdf`).
*   **Arguments:**
    *   `<results_dir>`: Directory containing GO enrichment TSV files.
    *   `<output_file>`: Path for the output PDF dotplot.
    *   `<plot_title>`: Title for the plot.
    *   `<top_n>`: Number of top GO terms.
    *   `<term_max>`: Maximum term size.
    *   `<pval>`: P-value threshold.
    *   `<fold>`: Fold enrichment threshold.
*   **Dependencies:** `ggplot2`, `tidyr`, `forcats`, `ggpubr`, `dplyr`, `readr`.

#### `Heatmap_figure3.R` and `Readexcels_and_Barfplot_figure3.R`

These are supplementary R scripts found in the root directory. They are not directly integrated into the `run_pipeline.sh` but might be used for additional analysis or visualization related to Figure 3. Their specific functions would require individual examination.

### Data

#### `data/Ribosome_Associated_proteome.tsv`

This file contains the raw ribosome-associated proteome data used as input for the analysis pipeline. It is expected to be a tab-separated values file with columns such as `file`, `log2fc_wt`, `pval_wt`, `gene_id`, and `cluster`, among others, as processed by the R and Python scripts.