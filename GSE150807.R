#&=================================================================================================*
#&  WRITTEN BY : Asmaa Yassin                                                                      *
#&  DATE       : April 17, 2025                                                                    *
#&-------------------------------------------------------------------------------------------------*
#& DESCRIPTION: This script performs differential expression analysis on RNA-seq data              *
#& from the GSE150807 using DESeq2. Volcano plot is generated using identified up and downregulated*
#& genes. GO and KEGG enrichment analyses are performed on the same genes. Results are saved as    *
#& CSV files and volcano plot as an image file                                                     *
#&=================================================================================================*

#Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "pathview"))
if (!require(ggplot2)) install.packages("ggplot2")

# Load libraries
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(pathview)


#set input files' paths
counts_path <- "path/to/GSE150807_raw_counts_filtered.csv"
metadata_path <- "path/to/GSE150807_metadata.csv"

# Load raw count matrix and metadata
raw_data <- read.csv(counts_path, row.names = 1)
metadata <- read.csv(metadata_path)

# Set row names and ensure sample alignment
rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition, levels = c("sensitive", "resistant"))

# Ensure the columns of raw data match the rownames of metadata
all(colnames(raw_data) == rownames(metadata))  # should return TRUE

# Round counts
raw_data <- round(raw_data)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_data, colData = metadata, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)

# Save full DESeq2 results
write.csv(res, "path/to/GSE150807_deseq2_results.csv")

# Filter significant genes (padj < 0.05 and abs(log2FC) > 1)
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_res, "path/to/GSE150807_significant_genes.csv")


# Volcano Plot

# Convert rownames into a column
res_df <- as.data.frame(res)
res_df$EntrezID <- rownames(res_df)


# Handle NAs
res_df$pvalue[is.na(res_df$pvalue)] <- 1
res_df$pvalue[res_df$pvalue == 0] <- 1e-300
res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
res_df$negLogP <- -log10(res_df$pvalue) # Calculate negative log10 p-values for volcano plot

# Label significance
res_df$Significant <- "NS"
res_df$Significant[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant"

# Annotate with gene symbols
res_df$GeneSymbol <- mapIds(org.Hs.eg.db,
                            keys = res_df$EntrezID,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
# Volcano plot
ggplot(res_df, aes(x=log2FoldChange, y=negLogP, color=Significant)) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  scale_color_manual(values=c("gray", "red")) +
  coord_cartesian(xlim=c(-5,5), ylim=c(0, 10)) +
  theme_minimal() +
  labs(title="Volcano Plot - GSE150807", x="Log2 Fold Change", y="-log10 p-value")
  
# Save result
write.csv(subset(res_df, Significant == "Significant"), 
          "path/to/GSE150807_volcano.csv.", 
          row.names = FALSE)

### Enrichment Analysis ###


# Separate upregulated and downregulated genes
up_sig <- subset(sig_res, log2FoldChange > 1)
down_sig <- subset(sig_res, log2FoldChange < -1)


# Convert gene IDs to Entrez IDs for upregulated & downregulated genes
up_entrez <- na.omit(mapIds(org.Hs.eg.db, 
                    keys = rownames(up_sig), 
                    column = "ENTREZID", 
                    keytype = "ENTREZID", 
                    multiVals = "first"))

down_entrez <- na.omit(mapIds(org.Hs.eg.db, 
                      keys = rownames(down_sig), 
                      column = "ENTREZID", 
                      keytype = "ENTREZID", 
                      multiVals = "first"))

# GO enrichment analysis for upregulated genes
go_up <- enrichGO(gene = up_entrez,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

#  GO enrichment analysis for downregulated genes
go_down <- enrichGO(gene = down_entrez,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)

# Plot the GO enrichment results
dotplot(go_up, showCategory = 15, title = "GO: Upregulated Genes (Biological Process)")
dotplot(go_down, showCategory = 15, title = "GO: Downregulated Genes (Biological Process)")

# Convert enrichment results to data frames and save to CSV
go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)
write.csv(go_up_df, "path/to/GSE150807_GO_upregulated_results.csv", row.names = FALSE)
write.csv(go_down_df, "path/to/GSE150807_GO_downregulated_results.csv", row.names = FALSE)



# KEGG Pathway Enrichment Analysis


# KEGG for upregulated genes 
kegg_up <- enrichKEGG(gene = up_entrez,
                     organism = 'hsa',
                     pvalueCutoff = 0.05)

# KEGG for downregulated genes 
kegg_down <- enrichKEGG(gene = down_entrez,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)

# Make the KEGG results more readable by mapping them back to gene symbols
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


# Plot the KEGG enrichment for upregulated and downregulated genes
dotplot(kegg_up, showCategory = 15, title = "KEGG: Upregulated Genes")
dotplot(kegg_down, showCategory = 15, title = "KEGG: Downregulated Genes")


# Save KEGG upregulated and downregulated results to CSV
write.csv(as.data.frame(kegg_up), "path/to/GSE150807_kegg_upregulated_results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_down), "path/to/GSE150807_kegg_downregulated_results.csv", row.names = FALSE)

cat("GSE150807 DESeq2 and enrichment analysis complete!\n")
