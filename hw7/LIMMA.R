BiocManager::install("limma")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")

library(limma)
library(Biobase)

geo_id <- 'GSE63885'
exp = read.csv(file = "C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/rna_seq_diff_exp/rna_seq_diff_exp/GSE63885/expression_for_limma2.csv", header = TRUE, row.names = 'Gene.Symbol')
ann = read.csv(file = "C:/Users/Дарья Алексеевна/Documents/МФТИ учёба/Bioinf_2026/hw7/rna_seq_diff_exp/rna_seq_diff_exp/GSE63885/annotation_for_limma.csv", header = TRUE, row.names = 'X')


exp_set <- ExpressionSet(assayData = as.matrix(exp), phenoData = AnnotatedDataFrame(ann))

slope <- factor(ann$clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1, levels = c("pCR",
                                                                                                                                                                     "pNC"), labels = c(1, 0))

pCR=as.integer(as.vector(slope))

iterse <- rep(1, length(slope))

design <- cbind(npCR=iterse,pCR=pCR)

fit <- lmFit(exp_set, design)
fit <- eBayes(fit)
top <-topTable(fit, coef="pCR", adjust="BH", n = Inf)

# добавление столбцов для LogFC = 1, 2, 3

results <- top

results$significant_FC1 <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
  "Yes", "No"
)

results$significant_FC2 <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 2,
  "Yes", "No"
)

results$significant_FC3 <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 3,
  "Yes", "No"
)

cat("количество значимых генов\n")
cat("|logFC| > 1:", sum(results$significant_FC1 == "Yes"), "\n")
cat("|logFC| > 2:", sum(results$significant_FC2 == "Yes"), "\n")
cat("|logFC| > 3:", sum(results$significant_FC3 == "Yes"), "\n")

write.csv(results, "LIMMA_results_all_thresholds.csv")

library(ggplot2)

counts_df <- data.frame(
  LogFC_threshold = c(1, 2, 3),
  Significant_genes = c(
    sum(results$significant_FC1 == "Yes"),
    sum(results$significant_FC2 == "Yes"),
    sum(results$significant_FC3 == "Yes")
  )
)

# Столбчатая диаграмма
p_bar <- ggplot(counts_df, aes(x = factor(LogFC_threshold), y = Significant_genes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Significant_genes), vjust = -0.5, size = 5) +
  labs(
    title = "LIMMA: Number of Significant Genes by LogFC Threshold",
    subtitle = "GSE63885: pCR vs pNC",
    x = "Log2 Fold Change Threshold",
    y = "Number of Significant Genes (adj.P.Val < 0.05)"
  ) +
  theme_minimal()

print(p_bar)

# Функция для построения volcano plot
plot_volcano_limma <- function(results, fc_threshold, title_extra = "") {
  results$log10padj <- -log10(results$adj.P.Val)
  
  results$log10padj[results$log10padj > 20] <- 20
  
  results$significant <- ifelse(
    results$adj.P.Val < 0.05 & abs(results$logFC) > fc_threshold,
    "Significant", "Not Significant"
  )
  
  top_genes <- rownames(results[order(results$adj.P.Val), ])[1:10]
  top_df <- results[rownames(results) %in% top_genes, ]
  
  ggplot(results, aes(x = logFC, y = log10padj, color = significant)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("grey", "red")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "blue", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed", color = "blue", alpha = 0.7) +
    labs(
      title = paste("Volcano Plot: pCR vs pNC (|LogFC| >", fc_threshold, ")", title_extra),
      x = "Log2 Fold Change (positive = higher in pCR)",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_df, 
                    aes(label = rownames(top_df)), 
                    size = 2.5, 
                    max.overlaps = 10)
}

library(ggrepel)

p1 <- plot_volcano_limma(results, 1)
p2 <- plot_volcano_limma(results, 2)
p3 <- plot_volcano_limma(results, 3)

print(p1)
print(p2)
print(p3)
