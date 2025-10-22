## ==== Libraries ====
suppressPackageStartupMessages({
  library(readr)
  library(tximport)
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(ggplot2)
})

## ==== Paths ====
proj_dir  <- "/Users/zainabbutul/maize_rnaseq_project"
quant_dir <- file.path(proj_dir, "QUANT")
tx2g_path <- file.path(proj_dir, "reference", "tx2gene_corrected.csv")
meta_path <- file.path(proj_dir, "Project-1_metadata.csv")
out_dir   <- file.path(proj_dir, "DESeq2_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ==== Metadata ====
metadata <- read.csv(meta_path, row.names = 1, check.names = FALSE)
metadata$tissue <- as.character(metadata$tissue)
metadata$Developmental_stage <- as.character(metadata$Developmental_stage)
metadata$growth_condition <- as.character(metadata$growth_condition)
metadata$condition <- ifelse(tolower(metadata$growth_condition) %in% c("well-watered","control"),
                             "control", "drought")
metadata$condition <- factor(metadata$condition, levels = c("control","drought"))
samples <- rownames(metadata)

## ==== Tx2gene mapping ====
tx2gene <- read_csv(tx2g_path, show_col_types = FALSE)
stopifnot(all(c("transcript","gene") %in% colnames(tx2gene)))

## ==== Kallisto files ====
files <- file.path(quant_dir, samples, "abundance.h5")
names(files) <- samples
missing <- files[!file.exists(files)]
if (length(missing)) {
  message("Missing kallisto files: ", paste(names(missing), collapse=", "))
  keep <- file.exists(files)
  files <- files[keep]
  metadata <- metadata[keep, , drop = FALSE]
  samples <- rownames(metadata)
}
stopifnot(length(files) > 0)

## ==== Import data ====
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

## ==== DESeq2 dataset ====
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)

## ==== Filter low counts globally ====
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

## ==== Save normalized counts ====
write.csv(counts(dds, normalized = TRUE), file.path(out_dir, "normalized_counts_all.csv"))

## ==== Global analysis ====
res_all <- results(dds, contrast = c("condition","drought","control"))
res_all <- lfcShrink(dds, coef = "condition_drought_vs_control", type = "apeglm")
write.csv(as.data.frame(res_all), file.path(out_dir, "DESeq2_results_ALL.csv"))

### PCA plot
vsd_all <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd_all, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "PCA_ALL.png"), plot = p, width = 6, height = 5)

### MA plot
png(file.path(out_dir, "MA_ALL.png"), 800, 600)
plotMA(res_all, ylim=c(-5,5))
dev.off()

### Volcano plot (overall)
png(file.path(out_dir, "Volcano_ALL.png"), 1200, 900)
EnhancedVolcano(as.data.frame(res_all),
                lab = rownames(res_all),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = "Drought vs Control (ALL)")
dev.off()

### Heatmap
topVarGenes <- head(order(rowVars(assay(vsd_all)), decreasing = TRUE), 50)
mat <- assay(vsd_all)[topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = metadata, filename = file.path(out_dir, "Heatmap_ALL.png"))

### Upregulated and Downregulated DEG lists (overall)
up_genes_all <- subset(as.data.frame(res_all), padj < pCut & log2FoldChange > FCcut)
down_genes_all <- subset(as.data.frame(res_all), padj < pCut & log2FoldChange < -FCcut)
write.csv(up_genes_all, file.path(out_dir, "DEGs_UP_ALL.csv"), row.names = TRUE)
write.csv(down_genes_all, file.path(out_dir, "DEGs_DOWN_ALL.csv"), row.names = TRUE)

message("✅ Overall upregulated and downregulated DEG lists saved")

## ==== Per tissue analysis  ====
tissues <- unique(metadata$tissue)
summ <- data.frame()

for (t in tissues) {
  message("🔍 Running tissue: ", t)
  
  idx <- metadata$tissue == t
  dds_sub <- dds[, idx]
  dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 1, drop = FALSE]
  
  if (length(unique(colData(dds_sub)$condition)) < 2) {
    message("Skipping ", t, " (not enough conditions)")
    next
  }
  
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition","drought","control"))
  res <- lfcShrink(dds_sub, coef = "condition_drought_vs_control", type = "apeglm")
  
  write.csv(as.data.frame(res), file.path(out_dir, paste0("DESeq2_results_", t, ".csv")))
  
  ### PCA
  vsd <- vst(dds_sub, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw(base_size = 14)
  ggsave(file.path(out_dir, paste0("PCA_", t, ".png")), plot = p, width = 6, height = 5)
  
  ### Volcano plot 
  res$padj[is.na(res$padj)] <- 1
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  
  sig_count <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1.5)
  pCut <- if (sig_count == 0) 0.1 else 0.05
  FCcut <- if (sig_count == 0) 0.5 else 1.5
  
  print(
    EnhancedVolcano(as.data.frame(res),
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = pCut,
                    FCcutoff = FCcut,
                    title = paste("Drought vs Control (", t, ")", sep = ""),
                    subtitle = paste("pCutoff:", pCut, "FCcutoff:", FCcut))
  )
  
  png(file.path(out_dir, paste0("Volcano_", t, ".png")), width = 1200, height = 900)
  print(
    EnhancedVolcano(as.data.frame(res),
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = pCut,
                    FCcutoff = FCcut,
                    title = paste("Drought vs Control (", t, ")", sep = ""),
                    subtitle = paste("pCutoff:", pCut, "FCcutoff:", FCcut))
  )
  dev.off()
  
  ### Heatmap
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  mat <- assay(vsd)[topVarGenes, ]
  mat <- mat - rowMeans(mat)
  pheatmap(mat, annotation_col = metadata[idx,,drop=FALSE],
           filename = file.path(out_dir, paste0("Heatmap_", t, ".png")))
  
  ### MA plot
  png(file.path(out_dir, paste0("MA_", t, ".png")), 800, 600)
  plotMA(res, ylim=c(-5,5))
  dev.off()
  
  ### DEG summary
  up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  tot <- sum(!is.na(res$padj))
  summ <- rbind(summ, data.frame(tissue = t,
                                 n_samples = sum(idx),
                                 control = sum(metadata$condition[idx] == "control"),
                                 drought = sum(metadata$condition[idx] == "drought"),
                                 genes_tested = tot,
                                 up = up,
                                 down = down))
}

write.csv(summ, file.path(out_dir, "DEG_summary_by_tissue.csv"))

message("✅ All analyses complete. Results saved in: ", out_dir)

# Filtering of DESEQ2 RESULTS & saving
input_path <- "/Users/zainabbutul/maize_rnaseq_project/DESeq2_analysis/DESeq2_results_ALL.csv"
output_dir <- "/Users/zainabbutul/maize_rnaseq_project/DESeq2_analysis"
output_path <- file.path(output_dir, "DEGs_padj0.05_absLFC1.csv")

# Load DESeq2 results
deseq <- read_csv(input_path)

# 4. Filter DEGs
filtered_DEGs <- deseq %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  mutate(
    abs_log2FC = abs(log2FoldChange),
    negLog10_padj = -log10(padj)
  )

# Summary 
total_genes <- nrow(deseq)
sig_genes <- nrow(filtered_DEGs)
upregulated <- sum(filtered_DEGs$log2FoldChange > 0, na.rm = TRUE)
downregulated <- sum(filtered_DEGs$log2FoldChange < 0, na.rm = TRUE)

cat("Total genes:", total_genes, "\n")
cat("Significant DEGs (padj < 0.05 & |log2FC| > 1):", sig_genes, "\n")
cat(" - Upregulated:", upregulated, "\n")
cat(" - Downregulated:", downregulated, "\n")

#  Save filtered table
write_csv(filtered_DEGs, output_path)
cat("\n✅ Filtered DEGs saved to:", output_path, "\n")

# preview first few rows
print(head(filtered_DEGs, 10))
