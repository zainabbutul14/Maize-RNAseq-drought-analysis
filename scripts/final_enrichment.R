# ==== Libraries ====
library(dplyr)
library(readr)
library(ggplot2)

# ==== Parameters ====
padj_cutoff <- 0.05
logfc_cutoff <- 1.5
topN <- 200  # fallback if too few DEGs

# ==== Paths ====
proj_dir <- "/Users/zainabbutul/maize_rnaseq_project"
deg_results_path <- file.path(proj_dir, "DESeq2_analysis", "DESeq2_results_ALL.csv")
gp_results_path <- "/Users/zainabbutul/Downloads/gProfiler_zmays_10-12-2025_9-19-09 PM__intersections.csv"

# ==== Load DESeq2 results ====
res <- read.csv(deg_results_path, row.names = 1)

# ==== Extract DEGs ====
deg_up <- res %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > logfc_cutoff) %>%
  rownames()

deg_down <- res %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -logfc_cutoff) %>%
  rownames()

deg_all <- c(deg_up, deg_down)

# ==== Fallback: Take top N by |log2FC| if no DEGs ====
if (length(deg_all) == 0) {
  deg_all <- rownames(res)[order(abs(res$log2FoldChange), decreasing = TRUE)][1:topN]
}

# ==== Save DEG list for enrichment analysis ====
write.table(deg_all,
            file.path(proj_dir, "DESeq2_analysis", "DEGs_for_enrichment.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
message("✅ DEG list saved: DEGs_for_enrichment.txt")

# ==== Load g:Profiler results ====
gp <- read_csv(gp_results_path)

# ==== Keep top 20 enriched terms by adjusted p-value ====
top_gp <- gp %>%
  arrange(adjusted_p_value) %>%
  head(20)

# ==== Barplot: -log10(adj p-value) ====
barplot <- ggplot(top_gp,
                  aes(x = reorder(term_name, -negative_log10_of_adjusted_p_value),
                      y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Top 20 Enriched GO Terms") +
  theme_minimal()

ggsave(file.path(proj_dir, "DESeq2_analysis", "GO_enrichment_barplot.png"),
       plot = barplot, width = 8, height = 6)
message("✅ GO enrichment barplot saved: GO_enrichment_barplot.png")

# ==== Dotplot of GO enrichment ====
dotplot <- ggplot(top_gp,
                  aes(x = reorder(term_name, -negative_log10_of_adjusted_p_value),
                      y = negative_log10_of_adjusted_p_value,
                      size = intersection_size,
                      color = term_size)) +
  geom_point() +
  coord_flip() +
  xlab("") +
  ylab("-log10(adjusted p-value)") +
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("GO Enrichment Dotplot") +
  theme_minimal()

ggsave(file.path(proj_dir, "DESeq2_analysis", "GO_enrichment_dotplot.png"),
       plot = dotplot, width = 8, height = 6)
message("✅ GO enrichment dotplot saved: GO_enrichment_dotplot.png")
