#&=================================================================================================*
#&  WRITTEN BY : Asmaa Yassin                                                                      *
#&  DATE       : April 17, 2025                                                                    *
#&-------------------------------------------------------------------------------------------------*
#& DESCRIPTION: This script performs differential expression analysis on RNA-seq data specifically *
#& DU145 and PC3 prostate cancer cell lines from GSE140440 using DESeq2.                           *
#& Volcano plots of up and downregulated genes are generated. GO and KEGG enrichment analyses      *
#& are performed for each cell line separately. Differential expression outputs, enrichment tables,*
#& and plots, are saved as CSV or image files.                                                     *
#&=================================================================================================*

# Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))
if(!require(ggplot2)) install.packages("ggplot2")

# Load libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Set file paths
counts_path <- "path/to/GSE140440_merged_matrix.csv"
metadata_path <- "path/to/GSE140440 metadata.csvv"

# Load merged raw counts and metadata
counts <- read.csv(counts_path, row.names = 1)
metadata <- read.csv(metadata_path)

colnames(metadata) <- c("Sample.ID", "Condition") # Ensure metadata column names are correct
metadata <- metadata[, !duplicated(colnames(metadata))] # Remove duplicate condition columns if they exist
rownames(metadata) <- metadata$Sample.ID #Set row names in metadata
metadata$Condition <- factor(metadata$Condition)
counts <- round(counts) # Ensure counts are integers (since DESeq2 requires raw integer counts)
metadata$Condition <- gsub("_\\d+$", "", metadata$Condition) # Ensure the condition column is properly grouped
metadata$Condition <- factor(metadata$Condition)

##############################
###   DU145 ANALYSIS       ###
##############################

# Subset DU145 samples
du145_metadata <- metadata[grep("DU145", metadata$Condition), ]
du145_counts <- counts[, rownames(du145_metadata)]

du145_metadata$Condition <- gsub("DU145_", "", du145_metadata$Condition) # Remove cell line prefix from Condition (so it becomes "res" or "sen")
du145_metadata$Condition <- factor(du145_metadata$Condition, levels = c("Sen", "Res")) # Set "Sen" (sensitive) as the reference

# Create DESeq2 dataset for DU145 and run DEA
dds_du145 <- DESeqDataSetFromMatrix(countData = du145_counts, colData = du145_metadata, design = ~ Condition)
dds_du145 <- DESeq(dds_du145)
res_du145 <- results(dds_du145)

# Save DU145 results
write.csv(res_du145, "path/to/GSE140440_deseq2_results_DU145.csv")

sig_res_du145 <- subset(res_du145, padj < 0.05 & abs(log2FoldChange) > 1) # Filter for significant genes (padj < 0.05 & |log2FoldChange| > 1)
write.csv(sig_res_du145, "path/to/GSE140440_significant_genes_DU145.csv")

# Volcano Plot for DU145
res_du145_df <- as.data.frame(res_du145)
res_du145_df$pvalue[is.na(res_du145_df$pvalue)] <- 1
res_du145_df$pvalue[res_du145_df$pvalue == 0] <- 1e-300
res_du145_df$negLogP <- -log10(res_du145_df$pvalue)
res_du145_df$log2FoldChange[is.na(res_du145_df$log2FoldChange)] <- 0
res_du145_df$Significant <- "NS"
res_du145_df$Significant[res_du145_df$padj < 0.05 & abs(res_du145_df$log2FoldChange) > 1] <- "Significant"
res_du145_df <- na.omit(res_du145_df)


ggplot(res_du145_df, aes(x = log2FoldChange, y = negLogP, color = Significant)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("gray", "red")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 10)) +
  theme_minimal() +
  labs(title = "Volcano Plot - DU145", x = "Log2 Fold Change", y = "-log10 p-value")
ggsave("path/to/GSE140440_volcano_plot_DU145.jpg")

# Merge volcano results with DESeq2 table to include gene names for DU145
deseq2_results_du145 <- read.csv("path/to/GSE140440_deseq2_results_DU145.csv", row.names = 1)
volcano_results_du145 <- res_du145_df[res_du145_df$Significant == "Significant", ]
deseq2_results_du145$Gene <- rownames(deseq2_results_du145)
merged_results_du145 <- merge(volcano_results_du145, deseq2_results_du145[, c("log2FoldChange", "padj", "Gene")], 
                              by = c("log2FoldChange", "padj"), all.x = TRUE)
write.csv(merged_results_du145, "path/to/GSE140440_volcano_DU145.csv", row.names = FALSE)


##############################
###   PC3 ANALYSIS         ###
##############################

# Subset PC3 samples
pc3_metadata <- metadata[grep("PC3", metadata$Condition), ]
pc3_counts <- counts[, rownames(pc3_metadata)]

pc3_metadata$Condition <- gsub("PC3_", "", pc3_metadata$Condition) # Remove cell line prefix from Condition (so it becomes "res" or "sen")
pc3_metadata$Condition <- factor(pc3_metadata$Condition, levels = c("Sen", "Res")) # Set "Sen" (sensitive) as the reference

# Run DESeq2 for PC3
dds_pc3 <- DESeqDataSetFromMatrix(countData = pc3_counts, colData = pc3_metadata, design = ~ Condition)
dds_pc3 <- DESeq(dds_pc3)
res_pc3 <- results(dds_pc3)

# Save PC3 results
write.csv(res_pc3, "path/to/GSE140440_deseq2_results_PC3.csv")

sig_res_pc3 <- subset(res_pc3, padj < 0.05 & abs(log2FoldChange) > 1) # Extract significant genes (padj < 0.05 & |log2FoldChange| > 1)
write.csv(sig_res_pc3, "path/to/GSE140440_significant_genes_PC3.csv")

# Volcano Plot for PC3
res_pc3_df <- as.data.frame(res_pc3)
res_pc3_df$pvalue[is.na(res_pc3_df$pvalue)] <- 1
res_pc3_df$pvalue[res_pc3_df$pvalue == 0] <- 1e-300
res_pc3_df$negLogP <- -log10(res_pc3_df$pvalue)
res_pc3_df$log2FoldChange[is.na(res_pc3_df$log2FoldChange)] <- 0
res_pc3_df$Significant <- "NS"
res_pc3_df$Significant[res_pc3_df$padj < 0.05 & abs(res_pc3_df$log2FoldChange) > 1] <- "Significant"
res_pc3_df <- na.omit(res_pc3_df)

ggplot(res_pc3_df, aes(x = log2FoldChange, y = negLogP, color = Significant)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("gray", "red")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 10)) +
  theme_minimal() +
  labs(title = "Volcano Plot - PC3", x = "Log2 Fold Change", y = "-log10 p-value")
ggsave("path/to/GSE140440_volcano_plot_PC3.jpg")

# Merge volcano results with gene names for PC3
deseq2_results_pc3 <- read.csv("path/to/GSE140440_deseq2_results_PC3.csv", row.names = 1)
volcano_results_pc3 <- res_pc3_df[res_pc3_df$Significant == "Significant", ]
deseq2_results_pc3$Gene <- rownames(deseq2_results_pc3)
merged_results_pc3 <- merge(volcano_results_pc3, deseq2_results_pc3[, c("log2FoldChange", "padj", "Gene")], 
                            by = c("log2FoldChange", "padj"), all.x = TRUE)
write.csv(merged_results_pc3, "path/to/GSE140440_volcano_PC3.csv", row.names = FALSE)


print("DESeq2 Analysis Complete for both DU145 and PC3! All results and plots have been saved.")


##############################
### ENRICHMENT ANALYSIS    ###
##############################


#### DU145 Analysis ###
sig_genes_du145 <- read.csv("path/to/GSE140440_significant_genes_DU145.csv", row.names = 1) # Load DU145 significant genes file
ensembl_ids_du145 <- gsub("\\..*", "", rownames(sig_genes_du145)) # Extract Ensembl IDs from row names and remove version numbers

# Convert cleaned Ensembl IDs to Entrez IDs
entrez_ids_du145 <- mapIds(org.Hs.eg.db, 
                           keys = ensembl_ids_du145, 
                           column = "ENTREZID", 
                           keytype = "ENSEMBL", 
                           multiVals = "first")

entrez_ids_du145 <- na.omit(entrez_ids_du145) # Remove missing values

# Print first few mapped Entrez IDs to verify
print(head(entrez_ids_du145))


# Gene Ontology (GO) Enrichment Analysis for DU145
go_du145 <- enrichGO(gene = entrez_ids_du145, 
                             OrgDb = org.Hs.eg.db, 
                             keyType = "ENTREZID",
                             ont = "BP",        # Biological Process
                             pAdjustMethod = "BH", 
                             pvalueCutoff = 0.05)

# Plot the GO enrichment results
dotplot(go_du145, showCategory = 15, title = "GO Biological Processes - DU145")

# Save GO enrichment results
write.csv(go_du145@result, 
          "path/to/GSE140440_GO_results_DU145.csv", 
          row.names = FALSE)

# KEGG Pathway Enrichment Analysis for DU145
kegg_du145 <- enrichKEGG(gene = entrez_ids_du145, 
                                 organism = "hsa", 
                                 pvalueCutoff = 0.05)

# Plot the KEGG pathway enrichment results
dotplot(kegg_du145, showCategory = 15, title = "KEGG Pathways - DU145")

# Save KEGG enrichment results
write.csv(kegg_du145@result, 
          "path/to/GSE140440_KEGG_results_DU145.csv", 
          row.names = FALSE)


### PC3 Analysis ###


# Load PC3 significant genes file (cell-line specific)
sig_genes_pc3 <- read.csv("path/to/GSE140440_significant_genes_PC3.csv", 
                          row.names = 1)

ensembl_ids_pc3 <- gsub("\\..*", "", rownames(sig_genes_pc3)) # Extract Ensembl IDs from row names and remove version numbers

# Convert Ensembl IDs to Entrez IDs
entrez_ids_pc3 <- mapIds(org.Hs.eg.db, 
                         keys = ensembl_ids_pc3, 
                         column = "ENTREZID", 
                         keytype = "ENSEMBL", 
                         multiVals = "first")

entrez_ids_pc3 <- na.omit(entrez_ids_pc3) # Remove missing values

# Print first few mapped Entrez IDs to verify
print(head(entrez_ids_pc3))

# Gene Ontology (GO) Enrichment Analysis for PC3
go_pc3 <- enrichGO(gene = entrez_ids_pc3, 
                           OrgDb = org.Hs.eg.db, 
                           keyType = "ENTREZID",
                           ont = "BP",        # Biological Process
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05)

# Plot the GO enrichment results
dotplot(go_pc3, showCategory = 15, title = "GO Biological Processes - PC3")

# Save GO enrichment results
write.csv(go_pc3@result, 
          "path/to/GSE140440_GO_results_PC3.csv", 
          row.names = FALSE)

# ---------------------------
# KEGG Pathway Enrichment Analysis for PC3
# ---------------------------
kegg_pc3 <- enrichKEGG(gene = entrez_ids_pc3, 
                               organism = "hsa", 
                               pvalueCutoff = 0.05)

# Plot the KEGG pathway enrichment results
dotplot(kegg_pc3, showCategory = 15, title = "KEGG Pathways - PC3")

# Save KEGG enrichment results
write.csv(kegg_pc3@result, 
          "path/to/GSE140440_KEGG_results_PC3.csv", 
          row.names = FALSE)
