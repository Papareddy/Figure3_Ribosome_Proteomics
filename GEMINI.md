# Ribosome associate proteome GO Enrichment & Differential Expression

## Project Goal
Visualize functional enrichment and expression trends for splicing-related subsets.

## Core Scripts
- `simplifiedGO.py`: Python logic for clustering/processing GO terms.
- `SimplifiedGO_ClusterBased.R`: R script for plotting clustered GO results.
- `Readexcels_and_Barfplot_figure3.R`: Data ingestion and initial bar chart generation.
- `Curved_volcanoplots_Figure3.R`: High-quality volcano plots with custom thresholds.
- `Heatmap_figure3.R`: Heatmap generation for expression signatures.

## Workflow
1. Run `simplifiedGO.py` to process raw enrichment data.
2. Use `SimplifiedGO_ClusterBased.R` to visualize functional groups.
3. Generate volcano and heatmaps for specific contrast groups.
