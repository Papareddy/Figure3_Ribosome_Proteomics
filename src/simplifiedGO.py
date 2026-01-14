import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler
import argparse

# Step 0: Parse arguments
parser = argparse.ArgumentParser(description="GO enrichment for Ribosome Associated Proteome.")
parser.add_argument("input_file", help="Path to Ribosome_Associated_proteome.tsv")
parser.add_argument("--outdir", default="GO_Results", help="Output directory")
parser.add_argument("--filter", default="_up", help="Suffix to filter clusters (e.g., _up or _down)")
args = parser.parse_args()

input_file = args.input_file
output_dir = args.outdir
os.makedirs(output_dir, exist_ok=True)

print(f"\n=== GO Enrichment Pipeline ===")
print(f"Input: {input_file} | Filter: {args.filter} | Output: {output_dir}\n")

# Step 1: Load data
data = pd.read_csv(input_file, sep='\t')

# Standardize column names to lowercase to avoid KeyErrors
data.columns = [c.lower() for c in data.columns]

# Step 2: Filter by cluster suffix
# FIX: added 'na=False' to handle the NaN values causing the ValueError
data_filtered = data[data['cluster'].str.endswith(args.filter, na=False)].copy()

if data_filtered.empty:
    print(f"No clusters found ending with '{args.filter}'. Check your TSV column names and content.")
    exit()

# Summarize cluster sizes
cluster_counts = data_filtered['cluster'].value_counts().sort_index()
cluster_counts.to_csv(os.path.join(output_dir, "cluster_sizes.tsv"), sep='\t')

# Plot cluster sizes
plt.figure(figsize=(10, 5))
sns.barplot(x=cluster_counts.index, y=cluster_counts.values, hue=cluster_counts.index, palette="tab10", legend=False)
plt.xticks(rotation=45, ha='right')
plt.title(f"Gene Counts: {args.filter} Clusters")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "cluster_sizes.png"), dpi=300)
plt.close()

# Step 3: GO term enrichment
gp = GProfiler(return_dataframe=True)
cluster_results = []

for cluster in data_filtered['cluster'].unique():
    print(f"Analyzing {cluster}...")
    # Using 'gene_id' based on your provided file header
    genes = data_filtered[data_filtered['cluster'] == cluster]['gene_id'].tolist()
    
    if not genes:
        continue
        
    result = gp.profile(
        organism='athaliana',
        query=genes,
        significance_threshold_method='fdr'
    )
    
    if not result.empty:
        result['Cluster'] = cluster
        cluster_results.append(result)

if cluster_results:
    go_df = pd.concat(cluster_results)
    go_df['-log10(p_value)'] = -np.log10(go_df['p_value'])
    go_df.to_csv(os.path.join(output_dir, "go_terms_all.tsv"), sep='\t', index=False)

    # Split and save by category
    for category, label in {"GO:BP": "bp", "GO:MF": "mf", "GO:CC": "cc"}.items():
        cat_df = go_df[go_df['source'] == category]
        if not cat_df.empty:
            cat_df.to_csv(os.path.join(output_dir, f"go_terms_{label}.tsv"), sep='\t', index=False)
            
            # Category Heatmap (Top 10 per cluster)
            top_cat = cat_df.groupby('Cluster').head(10)
            if len(top_cat['Cluster'].unique()) > 1:
                plt.figure(figsize=(12, 10))
                pivot_df = top_cat.pivot(index="name", columns="Cluster", values="-log10(p_value)").fillna(0)
                sns.heatmap(pivot_df, cmap="YlGnBu", cbar_kws={'label': '-log10(p-value)'})
                plt.title(f"Top GO {label.upper()} Terms")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"heatmap_{label}.png"), dpi=300)
                plt.close()
    print(f"\nSuccess. Results saved in {output_dir}")
else:
    print("No significant GO terms found.")