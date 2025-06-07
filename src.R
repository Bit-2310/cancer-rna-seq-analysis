# Gene Expression Analysis in Cancer using RNA-Seq
# Author: Pranava Upparlapalli
# Date: Sys.Date()

# Load required libraries
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(pROC)
library(caret)
library(parallel)
library(mltools)
library(MASS)
library(ggfortify)
library(stats)
library(data.table)
library(matrixStats)
library(cluster)
library(factoextra)

# Load data and merge with labels
data <- read.csv("data/data.csv", header = TRUE)
labels <- read.csv("data/labels.csv", header = TRUE)
merged_data <- merge(data, labels, by = "X", all.x = TRUE)
merged_data <- na.omit(merged_data)
merged_data$Class <- as.factor(merged_data$Class)
gc()

# Visualize sample distribution across cancer types
class_data <- as.data.frame(table(merged_data$Class))
colnames(class_data) <- c("Cancer_Type", "Count")
class_data$Count <- sort(class_data$Count, decreasing = TRUE)

ggplot(class_data, aes(x = Cancer_Type, y = Count, fill = Cancer_Type)) +
  geom_bar(stat = "identity") +
  theme_gray() +
  labs(title = "Distribution of Samples Across Cancer Types",
       x = "Cancer Type",
       y = "Number of Samples")
ggsave("plots/sample_distribution.png", width = 8, height = 5)

# Prepare DESeq2 inputs
gene_data <- merged_data[, !names(merged_data) %in% c("X", "Class")]
rownames(gene_data) <- merged_data$X
gene_expr_matrix <- round(t(gene_data))

sample_info <- data.frame(Class = merged_data$Class)
rownames(sample_info) <- merged_data$X

stopifnot(all(colnames(gene_expr_matrix) == rownames(sample_info)))
dds <- DESeqDataSetFromMatrix(countData = gene_expr_matrix,
                              colData = sample_info,
                              design = ~ Class)
dds <- DESeq(dds)
res <- results(dds)

# Filter significant genes
res_sig <- subset(res, padj < 0.05 & !is.na(padj))
sig_gene_names <- rownames(res_sig)
sig_genes_data <- as.data.frame(counts(dds, normalized = TRUE)[sig_gene_names, ])
sig_genes_data$Gene <- rownames(sig_genes_data)

# Filter and reshape expression data for clustering
gene_names <- sig_genes_data$Gene
numeric_data <- sig_genes_data[, -which(names(sig_genes_data) == "Gene")]
rownames(numeric_data) <- gene_names

is_row_valid <- apply(numeric_data, 1, function(row) {
  all(is.finite(row)) && var(row, na.rm = TRUE) > 0
})

gene_matrix_for_clustering <- numeric_data[is_row_valid, , drop = FALSE]

if (nrow(gene_matrix_for_clustering) >= 2) {
  gene_hclust <- hclust(dist(gene_matrix_for_clustering))
  ordered_gene_names <- rownames(gene_matrix_for_clustering)[gene_hclust$order]

  sig_genes_data_clean <- sig_genes_data[is_row_valid, ]
  sig_genes_data_long <- pivot_longer(sig_genes_data_clean,
                                      cols = -Gene,
                                      names_to = "Sample",
                                      values_to = "Expression")
  sig_genes_data_long$Gene <- factor(sig_genes_data_long$Gene,
                                     levels = ordered_gene_names)
} else {
  sig_genes_data_long <- data.frame()
}
gc()

# Volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("lightcoral", "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_gray()
ggsave("plots/volcano_plot.png", width = 8, height = 6)

# MA plot
ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("lightcoral", "darkgray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "MA Plot of Differentially Expressed Genes",
       x = "Average Expression (baseMean)", y = "Log2 Fold Change") +
  theme_gray()
ggsave("plots/ma_plot.png", width = 8, height = 6)

# PCA
vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_assay <- assay(vst_data)
vst_assay_filtered <- vst_assay[rowVars(vst_assay) > 0, ]

pca_results <- prcomp(t(vst_assay_filtered), center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_results$x)
pca_scores$Class <- colData(dds)$Class

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 1.5) +
  labs(title = "PCA of Cancer Samples", x = "PC1", y = "PC2") +
  theme_bw()

# K-means clustering
kmeans_result <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 5, nstart = 25)
pca_scores$K_Cluster <- factor(kmeans_result$cluster)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = K_Cluster)) +
  geom_point(size = 1.5) +
  labs(title = "K-means Clustering of Cancer Samples", x = "PC1", y = "PC2") +
  theme_bw()

# Silhouette plot
silhouette_analysis <- silhouette(kmeans_result$cluster, dist(pca_scores[, c("PC1", "PC2")]))
fviz_silhouette(silhouette_analysis) +
  theme_minimal() +
  labs(title = "Silhouette Plot for K-means Clustering")
ggsave("plots/silhouette_plot.png", width = 8, height = 6)

# Heatmap of top 50 variable genes
top_var_genes <- head(order(rowVars(assay(vst_data)), decreasing = TRUE), 50)
heatmap_matrix <- assay(vst_data)[top_var_genes, ]

annotation_col <- data.frame(CancerType = colData(dds)$Class)
rownames(annotation_col) <- colnames(heatmap_matrix)

pheatmap(heatmap_matrix,
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Top 50 Variable Genes",
         filename = "plots/heatmap_expression.png",
         width = 8,
         height = 6)

# Print session info
sessionInfo()
