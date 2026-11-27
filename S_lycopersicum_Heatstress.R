#Bioinformatics (transcriptomic analysis) ####
#Comparative transcriptome analysis of Solanumlycopersicum across temperature regimes
## Packages####
library(DESeq2)
library(Matrix)
library(dplyr)
library(tidyverse) 
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)  
library(ggplot2)
library(car)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

#### 1. LOAD DATA ####
# Transposed counts from Kallisto
counts <- t(readMM("output/counts_unfiltered/cells_x_genes.mtx"))
rownames(counts) <- readLines("output/counts_unfiltered/cells_x_genes.genes.names.txt")
colnames(counts) <- readLines("output/counts_unfiltered/cells_x_genes.barcodes.txt")
counts <- as.matrix(counts)

# Rename columns (order must match how samples were fed to Kallisto)
colnames(counts) <- c("LyC_1","LyC_2","LyC_3",
                      "PeC_1","PeC_2","PeC_3",
                      "LyW_1","LyW_2","LyW_3",
                      "PeW_1","PeW_2","PeW_3")

# Subset to S. lycopersicum only
ly_samples <- c("LyC_1", "LyC_2", "LyC_3", "LyW_1", "LyW_2", "LyW_3")
counts_ly  <- counts[, ly_samples]

# Cooler = 20/12∘C (day/night)
# Warm = 24/20∘C (day/night)

coldata_ly <- data.frame(
  condition = factor(rep(c("Cooler", "Warm"), each = 3),
                     levels = c("Cooler", "Warm")),
  row.names = ly_samples
)

# Verify
dim(counts_ly)
head(coldata_ly)

#### 2. BUILD DESeqDataSet ####
dds <- DESeqDataSetFromMatrix(
  countData = counts_ly,
  colData   = coldata_ly,
  design    = ~ condition
)

#### 3. FILTER LOW-COUNT GENES ####
keep <- rowSums(counts(dds)) > 200
dds  <- dds[keep, ]
cat("Genes retained:", nrow(dds), "\n")

#### 4. PCA ####
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
PCS_SLy <- ggplot(pcaData, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pct[1], "%")) +
  ylab(paste0("PC2: ", pct[2], "%")) +
  coord_fixed() + theme_minimal()

png("pca_lycopersicum.png", width = 10, height = 5, units = "in", res = 600)
print(PCS_SLy)
dev.off() 

#### 5. RUN DESeq ####
dds <- DESeq(dds)

#### 6. DISPERSION PLOT ####
plotDispEsts(dds)
png("dispersion_lycopersicum.png", width = 10, height = 5, units = "in", res = 600)
plotDispEsts(dds, main = "S. lycopersicum: Dispersion Estimates")
dev.off()

#### 7. CONTRAST ####
resultsNames(dds)
# "Intercept"
# "condition_Warm_vs_Cooler"

res <- results(dds, name = "condition_Warm_vs_Cooler",
               alpha = 1e-3, pAdjustMethod = "BH")

# Summary
table(res$padj < 0.05 & abs(res$log2FoldChange) > 1)

#### 8. MA PLOT ####
plotMA(res, main = "S. lycopersicum: Warm vs Cooler")
png("MA_lycopersicum.png", width = 10, height = 5, units = "in", res = 600)
plotMA(res, main = "S. lycopersicum: Warm vs Cooler")
dev.off()
#### 9. VOLCANO PLOT ####
library(EnhancedVolcano)
top_genes <- rownames(res)[order(res$pvalue)][1:30]

vol <- EnhancedVolcano(
  res,
  lab = rownames(res),
  selectLab = top_genes,
  title = "S. lycopersicum",
  subtitle = "+ Up in Warm\n- Up in Cooler",
  col = c("grey30", "gold", "royalblue", "red2"),
  x = "log2FoldChange", y = "pvalue",
  pCutoff = 0.05, FCcutoff = 1
)

png("volcano_lycopersicum.png", width = 10, height = 10, units = "in", res = 600)
print(vol)
dev.off()   

#### 10. HEATMAP ####
library(ComplexHeatmap)
library(RColorBrewer)

# Top 100 DEGs
deg_sig <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])
if (length(deg_sig) > 50) deg_sig <- deg_sig[order(res[deg_sig, "padj"])]
if (length(deg_sig) < 30)  deg_sig <- rownames(vsd)[order(rowVars(assay(vsd)), decreasing = TRUE)][1:100]
deg_sig <- head(deg_sig, 50)

tmp1 <- assay(vsd)[deg_sig, ]
tmp2 <- t(tmp1)
tmp3 <- scale(tmp2)
expr_scaled <- t(tmp3)

dim(expr_scaled)
# Should say: 50  6
# If it says NULL → the bug still fired → restart R and try again   

ha <- HeatmapAnnotation(
  condition = coldata_ly$condition,
  col = list(condition = c(Cooler = "skyblue", Warm = "darkorange"))
)

ht <- Heatmap(
  expr_scaled,
  name = "Expression (z-scored)",
  col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 9),
  border = TRUE,
  top_annotation = ha
)

png("heatmap_lycopersicum.png", width = 8, height = 10, units = "in", res = 600)
draw(ht)
dev.off()    

#### 11. ENRICHMENT ####
#"Which biological processes are overrepresented among the 
# genes that go UP when S. lycopersicum shifts from cooler to warmer?" 
library(clusterProfiler)
# Install (one-time)
if (!requireNamespace("AnnotationHub", quietly = TRUE))
  BiocManager::install("AnnotationHub")

library(AnnotationHub)
ah <- AnnotationHub()

# Find it
query(ah, 'org.Solanum_lycopersicum.eg.sqlite')   
# Should show: AH85737 | org.Solanum_lycopersicum.eg.sqlite

# Load it (downloads on first use, ~200 MB)
org.Sl.eg.db <- ah[["AH119847"]]   


# ID types it supports?
keytypes(org.Sl.eg.db)
# [1] "ACCNUM" "ALIAS"  "ENTREZID" "EVIDENCE" "EVIDENCEALL" "GENENAME"
# [7] "GID"    "GO"     "GOALL"    "ONTOLOGY"  "ONTOLOGYALL" "PMID"
# [13] "REFSEQ" "SYMBOL" "UNIGENE"   

library(clusterProfiler)

# Define DEGs
deg_up   <- rownames(res[res$padj < 0.05 & res$log2FoldChange > 1, ])
deg_down <- rownames(res[res$padj < 0.05 & res$log2FoldChange < -1, ])

# --- ID conversion (handle both LOC and symbol IDs) ---

# 1. LOC IDs → strip prefix → Entrez
loc_genes  <- deg_up[grepl("^LOC", deg_up)]
loc_entrez <- gsub("^LOC", "", loc_genes)

mapped_loc <- bitr(loc_entrez,
                   fromType = "ENTREZID",
                   toType   = c("ENTREZID", "SYMBOL"),
                   OrgDb    = org.Sl.eg.db)

# 2. Non-LOC IDs → try as SYMBOL or ALIAS
sym_genes  <- deg_up[!grepl("^LOC", deg_up)]

mapped_sym <- bitr(sym_genes,
                   fromType = "SYMBOL",
                   toType   = c("ENTREZID", "SYMBOL"),
                   OrgDb    = org.Sl.eg.db)

# If SYMBOL fails, try ALIAS
if (nrow(mapped_sym) == 0) {
  mapped_sym <- bitr(sym_genes,
                     fromType = "ALIAS",
                     toType   = c("ENTREZID", "SYMBOL"),
                     OrgDb    = org.Sl.eg.db)
}

# 3. Combine
mapped <- rbind(mapped_loc, mapped_sym)
mapped <- mapped[!duplicated(mapped$ENTREZID), ]

cat("Total up DEGs:", length(deg_up), "\n")
cat("Mapped to Entrez:", nrow(mapped), "\n")
cat("Unmapped:", length(deg_up) - nrow(mapped), "\n")

# --- Same for the universe (all genes) ---
#universe_loc  <- gsub("^LOC", "", rownames(res)[grepl("^LOC", rownames(res))])
universe_sym  <- rownames(res)[!grepl("^LOC", rownames(res))]

# Universe: LOC genes (already Entrez IDs after stripping)
universe_loc <- gsub("^LOC", "", rownames(res)[grepl("^LOC", rownames(res))])

# Keep only those that exist in the OrgDb
all_entrez <- keys(org.Sl.eg.db, keytype = "ENTREZID")
uni_loc <- universe_loc[universe_loc %in% all_entrez]

# Universe: non-LOC genes (symbols)
universe_sym <- rownames(res)[!grepl("^LOC", rownames(res))]
uni_sym <- bitr(universe_sym,
                fromType = "SYMBOL",
                toType   = c("ENTREZID", "SYMBOL"),
                OrgDb    = org.Sl.eg.db)

universe_entrez <- unique(c(uni_loc, uni_sym$ENTREZID))
cat("Universe size:", length(universe_entrez), "\n")   

ego_up <- enrichGO(gene      = mapped$ENTREZID,
                   universe  = universe_entrez,
                   OrgDb     = org.Sl.eg.db,
                   keyType   = "ENTREZID",
                   ont       = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

dotplot(ego_up, showCategory = 15) + ggtitle("Up in Warm — GO BP")   
png("GO_BP_up_in_warm.png", width = 10, height = 8, units = "in", res = 600)
print(dotplot(ego_up, showCategory = 15) + ggtitle("Up in Warm — GO BP"))
dev.off()   

ego_up$Description   
write.csv(as.data.frame(ego_up), "GO_BP_up_in_warm.csv", row.names = FALSE) 


