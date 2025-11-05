# Maize-RNAseq-drought-analysis
# Decoding Maize's Drought Response through RNA-seq  
*Differential Expression • Functional Enrichment • Network Visualization*

---

##  1. Overview / Snapshot
A comprehensive RNA-seq analysis of **maize (Zea mays)** under drought stress.  
This project decodes how maize reprograms its genes across **leaf, ear, and tassel** tissues to adapt and survive water deficit.  

Pipeline integrates:  
**Kallisto · DESeq2 · g:Profiler · STRING · Cytoscape**

> From raw FASTQ → Differential Expression → Functional Enrichment → Network Visualization

---

## 2. Why This Matters
Drought stress is one of the most severe global challenges impacting crop yield.  
Understanding maize’s molecular response under drought provides insights for:  
- Developing **drought-tolerant cultivars**  
- Identifying **key stress-regulatory genes**  
- Enhancing **agricultural resilience** in a changing climate  

> “Every plant faces a choice under drought — to grow or to survive.”

This project investigates how maize makes that molecular decision.

---

## 3. Data & Dataset Information
- **Organism:** *Zea mays* (inbred line B73)  
- **Conditions:** Drought vs. Control  
- **Tissues:** Leaf, Ear, Tassel  
- **Developmental Stages:** V12 and R1  
- **Replicates:** 3 per tissue-condition combination  
- **Total Samples:** 36 RNA-seq runs  
- **Source:** [NCBI BioProject PRJNA291919](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA291919)

---

## 4. Methods / Workflow Overview

### Pipeline Stages
1. **Data Acquisition** — Download 36 RNA-seq samples using SRA Toolkit.  
2. **QC & Trimming** — Evaluate and clean reads (FastQC, fastp).  
3. **Quantification** — Estimate transcript abundances (Kallisto).  
4. **Transcript–Gene Mapping** — Generate mapping file for DESeq2 (`tx2gene.csv`).  
5. **Differential Expression Analysis** — Identify DEGs (DESeq2).  
6. **Functional Enrichment** — Interpret DEGs biologically (g:Profiler).  
7. **Network Visualization** — Map interactions (STRING, Cytoscape).

### Tools Used
| Category | Tools / Packages |
|-----------|------------------|
| QC & Trimming | FastQC, fastp |
| Quantification | Kallisto |
| Differential Expression | DESeq2, tximport, ggplot2 |
| Enrichment | g:Profiler, clusterProfiler |
| Network Analysis | STRING, Cytoscape, UniProt Biomart |
| Languages | Bash, Python, R |

---

## 5. Results & Biological Insight

Each visualization unveils part of maize’s adaptive logic to drought stress.

- **PCA Plots:** Show distinct clustering between drought and control samples, strongest in leaf tissue.  
- **MA Plots:** Reveal transcriptional flexibility — leaf and ear adapt dynamically; tassel shows downregulation.  
- **Volcano Plots:** Highlight DEG magnitude and direction; leaf > ear > tassel in reprogramming intensity.  
- **Heatmaps:** Display co-expression blocks; leaf shows diverse modules, tassel shows uniform repression.  
- **GO Enrichment:** Enriched in “response to heat,” “ROS detoxification,” and “MAPK signaling.”  
- **STRING Network:** Reveals **HSP** and **ABA** signaling hubs coordinating proteostasis and hormonal control.  
- **Cytoscape Visualization:** Integrates STRING and enrichment results into a unified, interactive network — showing maize’s coordinated, multi-system drought defense strategy.

> “In Cytoscape, the maize of genes becomes a map — showing maize’s coordinated survival plan, where heat shock proteins and ABA signals connect to defend and conserve.”

📖 *For full plots, interpretations, and visuals → [Full Project Story] https://github.com/zainabbutul14/Maize-RNAseq-drought-analysis/blob/main/results/full_project_story.md *

---

## 6. Project Structure
```
Maize-RNAseq-drought-analysis/

├── data
│   ├── count_matrix.csv
│   ├── gProfiler_zmays_10-12-2025_9-19-09 PM__intersections.csv
│   ├── mart_export.txt
│   ├── Project-1_metadata.csv
│   ├── QC_stats.csv
│   ├── string_interactions_short.tsv
│   └── tx2gene_corrected.csv
├── results
│   ├── Plots
│   │   ├── CYTOSCAPE.pdf
│   │   ├── GO_enrichment_barplot.png
│   │   ├── GO_enrichment_dotplot.png
│   │   ├── Heatmap_ALL.png
│   │   ├── Heatmap_ear.png
│   │   ├── Heatmap_leaf.png
│   │   ├── Heatmap_tassel.png
│   │   ├── MA_ALL.png
│   │   ├── MA_ear.png
│   │   ├── MA_leaf.png
│   │   ├── MA_tassel.png
│   │   ├── PCA_ALL.png
│   │   ├── PCA_ear.png
│   │   ├── PCA_leaf.png
│   │   ├── PCA_tassel.png
│   │   ├── STRING NETWORK.png
│   │   ├── Top20_GO_KEGG.png
│   │   ├── Volcano_ALL.png
│   │   ├── Volcano_ear.png
│   │   ├── Volcano_leaf.png
│   │   └── Volcano_tassel.png
│   └── Tables
│       ├── DEG_summary_by_tissue.csv
│       ├── DEGs_DOWN_ALL.csv
│       ├── DEGs_for_enrichment.txt
│       ├── DEGs_padj0.05_absLFC1.csv
│       ├── DEGs_UP_ALL.csv
│       ├── DESeq2_results_ALL.csv
│       ├── DESeq2_results_ear.csv
│       ├── DESeq2_results_leaf.csv
│       ├── DESeq2_results_tassel.csv
│       └── normalized_counts_all.csv
└── scripts
    ├── cytoscape.py
    ├── final_deseq2_analysis.R
    ├── final_enrichment.R
    ├── gtf.py
    ├── qc_trimming.py
    ├── quanitifcation.py
    └── run_rnaseq_pipeline.sh
```

## Reproducibility
Run full RNA-seq pipeline using 
```
bash scripts/run_rnaseq_pipeline.sh
```

Differential expression
```
Rscript scripts/final_deseq2_analysis.R
```
Functional enrichment
```
Rscript scripts/final_enrichment.R
```
Cytoscape & STRING integration
```
python3 cytoscape.py
```
Then, used autoannotate in cytoscape to cluster and label modules
DEGs are mapped to **UniProt IDs** (via biomaRt in R) for STRING-based interaction analysis.

**STRING Database**  provides protein–protein networks that integrate directly into Cytoscape.

This step unifies functional enrichment and protein-interaction data — creating a systems-level drought-response network.

All outputs (plots, DEG tables, and enrichment results) are saved in results/.


# License

This repository is licensed under the MIT License.

“Drought teaches plants resilience.
This project taught me how data reveals that resilience — gene by gene.”
