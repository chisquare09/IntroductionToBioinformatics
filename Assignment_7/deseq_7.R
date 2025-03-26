library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)
countData <- read.csv('pasilla-countData.csv',header=TRUE,sep = ",")

metaData <- read.csv('pasilla-colData.csv',header=TRUE,sep = ",")

metaData$condition <- as.factor(metaData$condition)

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData = metaData,
                              design=~condition,
                              tidy = TRUE)
dds <- DESeq(dds)

res <- results(dds)

res_df <- as.data.frame(res)
res_df$neg_log10_padj <- -log10(res_df$padj)

# No threshold
total_DEGs <- sum(!is.na(res$padj) & res$padj < 0.05)
cat("Total differentially expressed genes before applying fold change threshold:", total_DEGs, "\n")

p0 <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = (padj < 0.05)), alpha = 0.6) +
  scale_color_manual(name="Type",
                     values = c("grey", "#B5EAD7"),
                     labels = c("Not significant","Significant")) +
  theme_minimal() +
  labs(title = "Volcano Plot (Before Applying Threshold)",
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# first threshold 
fc_threshold_1 <- log2(2)  # log2(2) = 1 (since 2-fold change is log2FC = ±1)
padj_threshold_1 <- 0.05

# Filter differentially expressed genes
deg_1 <- res_df %>%
  filter(abs(log2FoldChange) >= fc_threshold_1 & padj < padj_threshold_1)

# Count number of differentially expressed genes
num_deg_1 <- nrow(deg_1)
print(paste("Number of differentially expressed genes (Fold Change ≥ 2, p-adj < 0.05):", num_deg_1))


p1 <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = (abs(log2FoldChange) >= fc_threshold_1 & padj < padj_threshold_1)), alpha = 0.6) +
  scale_color_manual(name = "Type",
                     values = c("grey", "#FFB3BA"),
                     labels = c("Not significant","Significant")) +
  theme_minimal() +
  labs(title = "Volcano Plot (Fold Change ≥ 2, p-value < 0.05)",
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(padj_threshold_1), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold_1, fc_threshold_1), linetype = "dashed", color = "blue")

fc_threshold_2 <- log2(1.5)  # log2(1.5) ≈ 0.58
padj_threshold_2 <- 0.01

# second threshold
deg_2 <- res_df %>%
  filter(abs(log2FoldChange) >= fc_threshold_2 & padj < padj_threshold_2)


num_deg_2 <- nrow(deg_2)
print(paste("Number of differentially expressed genes (Fold Change ≥ 1.5, p-value < 0.01):", num_deg_2))

p2<-ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = (abs(log2FoldChange) >= fc_threshold_2 & padj < padj_threshold_2)), alpha = 0.6) +
  scale_color_manual(name = "Type",
                     values = c("grey", "#B3CDE3"),
                     labels = c("Not significant","Significant")) +
  theme_minimal() +
  labs(title = "Volcano Plot (Fold Change ≥ 1.5, padj < 0.01)",
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(padj_threshold_2), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-fc_threshold_2, fc_threshold_2), linetype = "dashed", color = "red")

plot <- p0/p1/p2
print(plot)

