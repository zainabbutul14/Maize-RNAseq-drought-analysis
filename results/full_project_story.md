# 📖 Full Project Story — Decoding Maize's Drought Response through RNA-seq

## 🌾 Abstract
Maize (*Zea mays*) faces a critical trade-off under drought: survival versus growth.  
Using RNA-seq data from NCBI BioProject PRJNA291919, this study maps the genome-wide transcriptional landscape across three tissues (leaf, ear, tassel) under drought and well-watered conditions.  
Through an integrated pipeline—Kallisto → DESeq2 → g:Profiler → STRING → Cytoscape—thousands of differentially expressed genes were identified and functionally characterized.  
Results reveal a hierarchical drought-response network: global repression of growth and metabolic genes, tissue-specific activation of stress-defense modules, and cross-tissue coordination between heat-shock (HSP) and ABA-signaling pathways.  
Together, these findings expose the molecular logic behind maize resilience and highlight key genetic targets for breeding drought-tolerant cultivars.

---

## 🌱 Background

Drought leaves plants at a crossroads — to invest in growth or to survive.  
This study explores the genetic logic behind that decision in *Zea mays* (maize), revealing how it reprograms its genes across tissues to withstand water stress.

> “Every plant faces a choice under drought — to grow or to survive.”

By integrating differential expression, enrichment, and protein-interaction network analysis, this project uncovers how maize coordinates molecular defense and energy-conservation strategies during drought.

Understanding these transcriptional shifts not only explains how maize survives drought but also provides a foundation for engineering crops that can thrive under climate stress.

---
## 🧬 Pipeline Overview

1. **Data Acquisition** — 36 RNA-seq samples (3 tissues × 2 stages × 2 conditions × 3 replicates) from [NCBI BioProject PRJNA291919](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA291919).  
2. **Quality Control & Trimming** — FastQC + fastp.  
3. **Quantification** — Transcript abundance estimation with **Kallisto**.  
4. **Transcript–Gene Mapping** — `tx2gene.csv` for DESeq2.  
5. **Differential Expression** — DESeq2 (FDR < 0.05 | |log2FC| > 1).  
6. **Functional Enrichment** — g:Profiler + clusterProfiler.  
7. **Network Visualization** — STRING + Cytoscape integration.

---
## 🧩 Biologically Insightful Visualization Outcomes

Each visualization reveals a different layer of maize’s molecular drought adaptation — from global expression reprogramming to functional network organization.

---
### 🧭 PCA Plots — *Global Expression Patterns Across Tissues*

**What it shows:** Principal Component Analysis summarizes expression variation across thousands of genes.  
**What it reveals:**  
- **Leaf:** Strong drought–control separation — major transcriptional reprogramming.  
- **Ear:** Moderate separation — partial adjustment.  
- **Tassel:** Overlapping clusters — limited plasticity.  
**Key takeaway:** Drought drives the strongest divergence in leaf, moderate in ear, minimal in tassel.

<img width="1800" height="1500" alt="PCA_tassel" src="https://github.com/user-attachments/assets/d719e4f2-f1b8-48b8-afcb-a818db5c6a7e" />
<img width="1800" height="1500" alt="PCA_leaf" src="https://github.com/user-attachments/assets/bfb654bb-6adc-403c-a7ca-4acc6631e7c8" />
<img width="1800" height="1500" alt="PCA_ear" src="https://github.com/user-attachments/assets/ffe6ef9a-50e9-44ac-b1ca-6c022a4df92a" />

*Figure 1. PCA plots by tissue — leaf shows the clearest drought–control separation.*

---
### 🌋 Volcano Plots — *Differential Expression Significance Across Tissues*

**What it shows:** log₂ fold-change vs −log₁₀ p-value for every gene.  
**What it reveals:**  
- **Leaf:** Thousands of DEGs — strong activation + repression.  
- **Ear:** Balanced regulation — survival vs growth trade-off.  
- **Tassel:** Sparse DEGs — global slowdown.  
**Key takeaway:** Drought-response gradient: **Leaf ≫ Ear ≫ Tassel.**

<img width="1200" height="900" alt="Volcano_tassel" src="https://github.com/user-attachments/assets/9564e278-009d-4da2-81c1-5a4fbcfb2d80" />
<img width="1200" height="900" alt="Volcano_leaf" src="https://github.com/user-attachments/assets/965a1f16-24b1-449f-be28-3a044f71cd03" />
<img width="1200" height="900" alt="Volcano_ear" src="https://github.com/user-attachments/assets/6ea112ab-5876-466a-a6d3-f3ccd8efe712" />
 
*Figure 2. Volcano plots show the progressive decline in DEG number and magnitude from leaf to tassel.*

---

### 📈 MA Plots — *Expression Balance and Fold Changes Across Tissues*

**What it shows:** Mean expression (A) vs log fold-change (M).  
**What it reveals:**  
- **Leaf:** Broad, balanced spread — flexible regulation.  
- **Ear:** Moderate spread — partial adaptation.  
- **Tassel:** Narrow spread — global repression.  
**Key takeaway:** The broader the spread, the greater the adaptive flexibility — leaf > ear > tassel.

<img width="800" height="600" alt="MA_tassel" src="https://github.com/user-attachments/assets/d5f9c077-054f-409d-89d8-f927fb8e229a" />
<img width="800" height="600" alt="MA_leaf" src="https://github.com/user-attachments/assets/59d3d25e-4261-4b66-a013-184ee3cfdac1" />
<img width="800" height="600" alt="MA_ear" src="https://github.com/user-attachments/assets/6767181e-b9a1-4776-9a30-24a1a0fbe60f" />

*Figure 3. MA plots illustrating tissue-specific transcriptional responsiveness.*

---

### 🔥 Heatmaps — *Tissue-specific Co-expression Patterns*

**What it shows:** Normalized expression clustering of DEGs across samples.  
**What it reveals:**  
- **Leaf:** Mixed activation + repression of large gene sets.  
- **Ear:** Two dominant clusters — upregulated stress vs downregulated growth.  
- **Tassel:** Uniform blue — broad repression.  
**Key takeaway:** Distinct tissue strategies but a shared energy-saving drought module.

<img width="2100" height="2100" alt="Heatmap_leaf" src="https://github.com/user-attachments/assets/1f82adc0-b92d-4b43-ab21-583f2ce468df" />
<img width="2100" height="2100" alt="Heatmap_ear" src="https://github.com/user-attachments/assets/518ac547-43c3-4020-b7d9-09fc6503afcd" />
<img width="2100" height="2100" alt="Heatmap_tassel" src="https://github.com/user-attachments/assets/773f4146-c25d-4cd3-bba9-0c5642081611" />

*Figure 4. Heatmaps showing distinct co-expression signatures per tissue and a conserved downregulated core.*

---

## 🔎 Combined (ALL Tissues) Visualizations — Shared & Core Responses

Combined plots summarize **cross-tissue trends** and highlight globally conserved drought-response signatures.

### 🧭 PCA_ALL — *Global sample structure across tissues*
**What it shows:** PCA using all samples together.  
**What it reveals:** Tissue identity drives most variance, yet drought–control separation is also visible across tissues.  
**Key takeaway:** A common drought axis overlays tissue-specific expression space.

<img width="1800" height="1500" alt="PCA_ALL" src="https://github.com/user-attachments/assets/80152c07-aac4-435e-a294-6fdf4f588018" />
  
*Figure 5. Combined PCA across all tissues showing both tissue clustering and shared drought-driven separation.*

---

### 📈 MA_ALL — *Overall Differential-Expression Landscape*
**What it shows:** Average expression vs log₂ fold-change for all samples combined.  
**What it reveals:** Dense DEG bands mark consistently regulated genes; outliers represent tissue-specific shifts.  
**Key takeaway:** A core set of genes is consistently regulated across tissues — universal drought responders.

<img width="800" height="600" alt="MA_ALL" src="https://github.com/user-attachments/assets/1eff6ffc-5e9d-4cf6-a7c2-8416627fded6" />

*Figure 6. Combined MA plot highlighting genes with consistent regulation across tissues.*

---

### 🌋 Volcano_ALL — *Shared Significant DEGs Across Tissues*
**What it shows:** Significance vs fold-change aggregated across tissues.  
**What it reveals:** High-confidence DEGs that remain significant in multiple tissues.  
**Key takeaway:** Identifies robust, cross-tissue drought-responsive genes for follow-up studies.

<img width="1200" height="900" alt="Volcano_ALL" src="https://github.com/user-attachments/assets/a3cdf06a-25e0-496d-94bf-a0b29f85c505" />

*Figure 7. Combined volcano isolating strong DEGs conserved across tissue contexts.*

---

### 🔥 Heatmap_ALL — *Core Co-expression Modules Across the Whole Dataset*
**What it shows:** Expression of top variable/core DEGs across all samples.  
**What it reveals:** Conserved clusters — energy-saving genes repressed everywhere, defense genes activated broadly.  
**Key takeaway:** Reveals a compact “core drought module” shared across tissues.

  <img width="2100" height="2100" alt="Heatmap_ALL" src="https://github.com/user-attachments/assets/37234920-5b6c-435f-b6ed-a1e83fe84825" />

*Figure 8. Combined heatmap revealing shared co-expression modules underlying universal drought adaptation.*

---

### 🧠 GO Enrichment — *Dominant Functional Pathways*

**What it shows:** Enriched biological processes among DEGs.  
**What it reveals:** Activation of stress-signaling, antioxidant defense, and protein-folding mechanisms.  
**Key takeaway:** Heat, ROS detoxification, and MAPK signaling dominate drought adaptation.

<img width="2400" height="1800" alt="GO_enrichment_barplot" src="https://github.com/user-attachments/assets/8be6b782-734a-40d1-8782-6ebb7a835514" />

*Figure 9. GO enrichment analysis identifying stress and signaling pathways as key drought mechanisms.*

---

### 🕸️ STRING Network — *Protein–Protein Interaction Modules*

**What it shows:** Functional relationships among drought-responsive proteins.  
**What it reveals:** Two major hubs — **HSP proteostasis** and **ABA signaling**.  
**Key takeaway:** Drought triggers coordinated modules that control protein stability and water regulation.

 <img width="654" height="551" alt="STRING NETWORK" src="https://github.com/user-attachments/assets/fa3ccd01-c975-4d64-b45c-34ae1ae02707" />

*Figure 10. STRING network depicting HSP and ABA signaling hubs as central drought-response regulators.*

---

### 🧬 Cytoscape Visualization — *Functional Network Mapping*

**What it shows:** Cytoscape integrates STRING + enrichment data into an annotated interaction map.  
**What it reveals:** Clusters of co-functional genes tied to proteostasis, oxidative stress, and hormonal signaling.  
**Key takeaway:** Confirms maize’s drought response operates as an interconnected, multi-system network.

[CYTOSCAPE.pdf](https://github.com/user-attachments/files/23037559/CYTOSCAPE.pdf)

*Figure 11. Cytoscape map integrating enrichment and STRING results — AutoAnnotate clustering reveals Heat Shock and ABA modules.*

🧠 *AutoAnnotate identified two key hubs:*  
1. **Heat Shock Cluster:** HSP82 / HSP26 — protein protection.  
2. **ABA Signaling Cluster:** PYL7–PP2C–SnRK2 axis — osmotic control.

---

## 🌾 Biological Conclusion

Maize displays a **hierarchical drought-response network**:  
- **Global:** Growth repression and energy conservation.  
- **Leaf:** Balances photosynthesis and defense.  
- **Ear:** Activates protection while repressing reproduction.  
- **Tassel:** Broad shutdown to preserve resources.

🧠 Integration of HSP and ABA modules shows a two-pronged strategy:  
1. Maintain protein stability through proteostasis networks.  
2. Regulate water balance via ABA-mediated signaling.

🌿 *Together, these findings highlight maize’s molecular resilience and identify targets for future drought-tolerant crop improvement.*

---

## 💧 Interpreting Drought vs. Well-Watered Differentiation

The multi-tissue transcriptome comparison clearly distinguishes **drought-stressed** and **well-watered** maize plants at multiple biological levels:

- **Expression Profiles (PCA, MA, Volcano, Heatmaps):**  
  Under drought, the **leaf** shows distinct transcriptional separation from controls, driven by repression of photosynthetic genes and activation of defense and ABA-responsive genes.  
  The **ear** demonstrates moderate differentiation — balancing reproductive growth with stress survival.  
  The **tassel** remains transcriptionally conservative, showing minimal deviation from well-watered controls.

- **Functional Pathways (GO, STRING):**  
  Drought stress activates pathways associated with **heat response**, **ROS detoxification**, **protein folding**, and **hormonal signaling**, while downregulating those linked to **cell growth** and **metabolism**.  
  Well-watered plants maintain high expression of **growth-promoting** and **metabolic** genes that are suppressed under drought.

- **Molecular Networks (Cytoscape):**  
  Network topology reveals that drought-stressed plants form new regulatory clusters — notably **HSP-based proteostasis modules** and **ABA signaling nodes** — absent or weakly connected in well-watered plants.  
  This shift marks a systemic reorganization from **growth optimization** to **stress defense coordination**.

---
:

🌻 This integration of computational genomics and plant biology reveals how big data can illuminate the invisible resilience of crops — and guide us toward more sustainable agriculture.

## 🌿 Translational & Treatment Insight

Understanding these transcriptional differentiations opens potential routes for improving drought tolerance:

1. **Molecular Breeding:**  
   Target upregulated drought-responsive genes (e.g., *HSP82*, *HSP26*, *PYL7*, *PP2C*) as biomarkers for selecting stress-resilient genotypes.

2. **Genetic Engineering:**  
   Overexpress key drought-defense genes (from HSP or ABA clusters) to strengthen proteostasis and water-use regulation under field droughts.

3. **Agronomic Strategies:**  
   Integrate physiological monitoring with transcript-level markers — early activation of ABA- and ROS-related genes could guide optimized irrigation or stress-priming schedules.

4. **Predictive Genomics:**  
   The conserved “core drought module” identified across tissues can be used to build predictive models of stress tolerance for different maize cultivars.

---

🧠 *In essence, drought-stressed maize reprograms itself at the molecular level — suppressing growth to survive, while well-watered plants sustain developmental expansion. Understanding this switch enables targeted interventions to engineer or breed maize that can both grow and endure.*
