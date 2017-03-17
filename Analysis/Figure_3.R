#############################################################
#### Code to generate Figure 3 ##############################
#############################################################
library(Rtsne)
library(pheatmap)
library(ggplot2)
library(BASiCS)
setwd("/Users/eling01/GitHub/ImmuneAging2017/")

# Load gene names
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

# Load data
exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")

### Figure 3A
### tSNE of B6 and CAST naive/active
Data <- exp.data[,which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)) |
                                 grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)) |
                                 grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)) |
                                 grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]

# Calculate tSNE 
tsne <- Rtsne(t(log10(Data + 1)), perplexity = 30)

# Color data points by activation state and species
state <- vector(length = ncol(Data))
state[which(grepl("SS51_naive", colnames(Data)) | grepl("SS52_naive", colnames(Data)))] <- "B6_naive" 
state[which(grepl("SS51_active", colnames(Data)) | grepl("SS52_active", colnames(Data)))] <- "B6_active"
state[which(grepl("SS19_naive", colnames(Data)) | grepl("SS25_naive", colnames(Data)))] <- "CAST_naive" 
state[which(grepl("SS19_active", colnames(Data)) | grepl("SS25_active", colnames(Data)))] <- "CAST_active" 

tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], State = state)

# Plot tSNE
ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 10, aes(col = state), shape = "-") +
  scale_colour_manual(values = c("black", "coral4", "tan3", "red")) + 
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

### Figure 3B
# Find shared activation genes
# Load BASiCS Test objects for following comparisons
# B6 young naive vs B6 young active
DE.B6.young <- Test_logFC2_B6
DE.B6.young.active <- DE.B6.young[DE.B6.young$ResultDiffExp == "SS51_active_SS52_active+",]

# CAST young naive vs CAST young active
DE.CAST.young <- Test_logFC2_CAST
DE.CAST.young.active <- DE.CAST.young[DE.CAST.young$ResultDiffExp == "SS19_active_SS25_active+",]

# Collect shared activation genes
shared.a <- intersect(as.character(DE.B6.young.active$GeneNames), as.character(DE.CAST.young.active$GeneNames))

# To find species specific activation genes, we only test the genes that are upregulated in either B6 or CAST
# Load this specifc test: B6 young active vs CAST young active
spec.test <- Test_logFC2_active
B6.spec.a <- as.character(spec.test$GeneNames[spec.test$ResultDiffExp == "SS51_active_SS52_active+"])

CAST.spec.a <- as.character(spec.test$GeneNames[spec.test$ResultDiffExp == "SS19_active_SS25_active+"])

### We remove DE genes that could be cause by mapping artifatcs
# We test differences in mapping comparing B6 young active mapped to B6 vs mapped to CAST
# and CAST young active mapped to B6 vs mapped to CAST
# Again, we only test the genes that are upregulated during activation
DE.B6.map <- Test_logFC2_activationGenes_CAST
DE.CAST.map <- Test_logFC2_activationGenes_B6

genes <- c(DE.B6.map$GeneNames[DE.B6.map$ResultDiffExp != "NoDiff" & DE.B6.map$ResultDiffExp != "ExcludedByUser"], DE.CAST.map$GeneNames[DE.CAST.map$ResultDiffExp != "NoDiff" & DE.CAST.map$ResultDiffExp != "ExcludedByUser"])
genes <- unique(genes)

# Species specifc genes
B6.spec.a <- B6.spec.a[-which(!is.na(match(B6.spec.a, genes)))]
CAST.spec.a <- CAST.spec.a[-which(!is.na(match(CAST.spec.a, genes)))]

# Plot heatmaps
# Downsample to 30 cells per condition 
Data <- exp.data[,c(which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))[sample(1:length(which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))), 30)],
                    which(grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)))[sample(1:length(which(grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)))), 30)],
                    which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))[sample(1:length(which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))), 30)],
                    which(grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))[sample(1:length(which(grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))), 30)]
)]

# Proportionally downsample genes
genes.shared <- shared.a[sample(1:length(shared.a), 100)]
genes.B6 <- B6.spec.a[sample(1:length(B6.spec.a), 38)]
genes.CAST <- CAST.spec.a[sample(1:length(CAST.spec.a), 34)]

# plot heatmaps
pheatmap(log10(Data[genes.shared[order(rowMeans(Data[genes.shared,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.shared,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))
pheatmap(log10(Data[genes.B6[order(rowMeans(Data[genes.B6,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.B6,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))
pheatmap(log10(Data[genes.CAST[order(rowMeans(Data[genes.CAST,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.CAST,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))

### Figure 3C
# Visualize GO categories
# GO anlysis was perfomed by using the online tool DAVID based on the shared and species specific activation genes
# The results can be read in and plotted
results <- read.table("path/to/GOfile", header = TRUE, sep = '\t')
results <- results[results$Bonferroni < 0.1,]

par(mar = c(16,4,1,4))
barplot(height = -log10(results$Bonferroni), col = "white", names.arg = results$Term, las = 2, cex.names = 0.5, ylim = c(0,5))
abline(a = 1, b = 0, col = "red", lwd = 2)

# For similar GO categories only the more significant was presented

### Figure 3D
# Frequencies of genes expression - downsample to 70 genes
genes.shared <- shared.a[sample(1:length(shared.a), 70)]
genes.B6 <- B6.spec.a[sample(1:length(B6.spec.a), 70)]
genes.CAST <- CAST.spec.a[sample(1:length(CAST.spec.a), 70)]

shared.cells <- exp.data[genes.shared,
                         which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)) | grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]
B6.cells <- exp.data[genes.B6,
                     which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)) | grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]
CAST.cells <- exp.data[genes.CAST,
                       which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)) | grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]


plot(hist(apply(shared.cells[sample(1:nrow(shared.cells), 70),], 1, function(n){length(which(n > 0))/ncol(shared.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(shared.cells, 1, function(n){length(which(n > 0))/ncol(shared.cells)})))

plot(hist(apply(B6.cells[sample(1:nrow(B6.cells), 70),], 1, function(n){length(which(n > 0))/ncol(B6.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(B6.cells, 1, function(n){length(which(n > 0))/ncol(B6.cells)})))

plot(hist(apply(CAST.cells[sample(1:nrow(CAST.cells), 70),], 1, function(n){length(which(n > 0))/ncol(CAST.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(CAST.cells, 1, function(n){length(which(n > 0))/ncol(CAST.cells)})))

