import pandas as pd
import numpy as np
import re

# Define file names
gprofiler_file = "/Users/zainabbutul/Downloads/gProfiler_zmays_10-12-2025_9-19-09 PM__intersections.csv"
degs_file = "/Users/zainabbutul/maize_rnaseq_project/DESeq2_analysis/DEGs_padj0.05_absLFC1.csv"

# --- 1. Load and Prepare Data ---

# Load gProfiler results
gprof_df = pd.read_csv(gprofiler_file)

# Load DEG results
degs_df = pd.read_csv(degs_file)

# Rename the first column (gene ID) and convert IDs to uppercase
degs_df.rename(columns={degs_df.columns[0]: 'id'}, inplace=True)
degs_df['id'] = degs_df['id'].str.upper()

# --- 2. Generate the Edge List (GO Term <--> Gene ID) ---

edge_list_rows = []
for index, row in gprof_df.iterrows():
    go_term = row['term_id']
    # Split the comma-separated string of genes, handling optional spaces
    genes = re.split(r',\s*', str(row['intersections']))

    for gene_id in genes:
        if gene_id:
            edge_list_rows.append({'Source': go_term, 'Target': gene_id})

cytoscape_edges_df = pd.DataFrame(edge_list_rows)

# --- 3. Generate the Node Attributes File ---

# A. GO Term Nodes
go_nodes_df = gprof_df[[
    'term_id', 'term_name', 'adjusted_p_value', 'source', 'negative_log10_of_adjusted_p_value'
]].copy()
go_nodes_df.rename(columns={
    'term_id': 'id',
    'negative_log10_of_adjusted_p_value': 'negLog10_padj'
}, inplace=True)
go_nodes_df['type'] = 'GO_term'
# Re-index columns to prepare for concatenation with gene data (setting gene attributes to NaN)
go_nodes_df = go_nodes_df.reindex(columns=[
    'id', 'type', 'term_name', 'adjusted_p_value', 'source', 'negLog10_padj',
    'log2FoldChange', 'padj', 'abs_log2FC'
])

# B. Gene Nodes
gene_nodes_df = degs_df[[
    'id', 'log2FoldChange', 'padj', 'abs_log2FC'
]].copy()
gene_nodes_df['type'] = 'gene'
# Re-index columns to prepare for concatenation with GO data (setting GO attributes to NaN)
gene_nodes_df = gene_nodes_df.reindex(columns=[
    'id', 'type', 'term_name', 'adjusted_p_value', 'source', 'negLog10_padj',
    'log2FoldChange', 'padj', 'abs_log2FC'
])

# C. Combine and Finalize Node Attributes
cytoscape_nodes_df = pd.concat([go_nodes_df, gene_nodes_df], ignore_index=True)
cytoscape_nodes_df.fillna('', inplace=True)

# --- 4. Write Cytoscape-Ready Files ---
cytoscape_nodes_df.to_csv("cytoscape_nodes_output.csv", index=False)
cytoscape_edges_df.to_csv("cytoscape_edges_output.csv", index=False)