if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")  
BiocManager::install("DESeq2")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")


library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)


expr_raw <- read.csv("C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/rna_seq_diff_exp/rna_seq_diff_exp/raw_counts_ici_samples.tsv", sep='\t', row.names = 1)

meta = data.frame(RNA = colnames(expr_raw))
data = expr_raw
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ 1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/data/normalized_counts_ici_samples.tsv", sep="\t", quote=F, col.names=NA)


sample_metadata <- read.csv('C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/rna_seq_diff_exp/rna_seq_diff_exp/meta_responses.tsv', row.names = 1, sep='\t')
sample_metadata <- sample_metadata %>%
  filter(X0 %in% c('R', 'NR'))
rownames(sample_metadata) <- gsub("-", ".", rownames(sample_metadata))

normalized_counts <- read.csv('C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/data/normalized_counts_ici_samples.tsv', sep='\t', row.names = 1)

normalized_counts <- normalized_counts[row.names(sample_metadata)]

normalized_counts <- round(normalized_counts)

dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = sample_metadata,
                              design = ~ X0)

# Стабилизирующее преобразование дисперсии (для PCA)
vsd <- vst(dds, blind = TRUE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "X0", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = X0, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, vjust = 1.5) +
  labs(
    title = "PCA: Responders vs Non-responders",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("R" = "blue", "NR" = "red"), name = "Response")


# удаляем выброс
sample_metadata <- sample_metadata[rownames(sample_metadata) != "LuC_56_S19_R1_001ReadsPerGene", , drop = FALSE]
normalized_counts <- normalized_counts[, rownames(sample_metadata)]


dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = sample_metadata,
                              design = ~ X0)

#Clustermap
install.packages("pheatmap")
vsd <- vst(dds, blind = TRUE)
library(pheatmap)

sample_dist <- as.matrix(dist(t(assay(vsd))))
annotation_df <- data.frame(Response = colData(dds)$X0)
rownames(annotation_df) <- colnames(dds)

pheatmap(sample_dist,
         clustering_distance_rows = as.dist(sample_dist),
         clustering_distance_cols = as.dist(sample_dist),
         annotation_col = annotation_df,
         annotation_row = annotation_df,
         main = "Clustermap образцов",
         fontsize_row = 8,
         fontsize_col = 8)

dds$X0 <- relevel(dds$X0, ref = "NR")
dds <- DESeq(dds)

res <- results(dds)
summary(res)

# volcano plot
res_df <- as.data.frame(res)

res_df$significant <- ifelse(
  !is.na(res_df$padj) &           # отфильтровываем NA
    res_df$padj < 0.05 & 
    abs(res_df$log2FoldChange) > 1,
  "Yes", "No"
)

res_df$log10padj <- -log10(res_df$padj)

res_df$log10padj[!is.na(res_df$log10padj) & res_df$log10padj > 20] <- 20

plot_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$log10padj), ]

top_genes <- rownames(plot_df[order(plot_df$padj), ])[1:10]
top_genes_df <- plot_df[rownames(plot_df) %in% top_genes, ]

ggplot(plot_df, aes(x = log2FoldChange, y = log10padj, color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(
    title = "Volcano Plot: Responders vs Non-responders",
    x = "Log2 Fold Change (Positive = higher in R)",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(data = top_genes_df, 
                  aes(label = rownames(top_genes_df)), 
                  size = 3, 
                  max.overlaps = 10)




# далее код с семинара

plotMA(res, main="DESeq2 MA Plot")

pvalue_threshold <- 10**(-2)
log2fc_threshold <- 3

res$significance <- ifelse(res$padj < pvalue_threshold & abs(res$log2FoldChange) > log2fc_threshold, 
                           "Significant", "Not Significant")

res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$padj <- as.numeric(res$padj)


hgnc_table <- read.csv('C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/rna_seq_diff_exp/rna_seq_diff_exp/hgnc_complete_set.txt', row.names = 1, sep='\t')
symbol_map <- setNames(hgnc_table$symbol, hgnc_table$ensembl_gene_id)

rownames(res) <- sapply(rownames(res), function(id) {
  if (id %in% names(symbol_map) && !is.na(symbol_map[id])) {
    symbol_map[id]
  } else {
    id
  }
})

significant_genes <- subset(res, padj < pvalue_threshold & abs(log2FoldChange) > log2fc_threshold)

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +  # Non-significant genes in grey, significant in red
  labs(title = "Control Subculture", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(data = significant_genes, 
                  aes(label = rownames(significant_genes)), size = 3, max.overlaps = 10) +
  xlim(-10, 10)

write.csv(as.data.frame(res), "C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/data/differential_expression_results_ICI.csv")

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  
  # ДОБАВЬТЕ ЭТИ ТРИ ЛИНИИ:
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
             linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = -log10(pvalue_threshold), 
             linetype = "dashed", color = "blue", alpha = 0.7) +
  
  labs(title = "Control Subculture", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(data = significant_genes, 
                  aes(label = rownames(significant_genes)), 
                  size = 3, max.overlaps = 10) +
  xlim(-10, 10)
