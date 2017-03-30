#############################################################
#### Code to generate Figure 1 ##############################
#############################################################
library(Rtsne)
library(pheatmap)
library(ggplot2)
library(BASiCS)
library(caroline)
setwd("/Users/eling01/GitHub/ImmuneAging2017/")

# Load gene names
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

### Figure 1
# Fig 1B
exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")

Data <- exp.data[,which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)) |
                           grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)))]

# Colour by species
species <- vector(length = ncol(Data))
species[which(grepl("SS51_naive", colnames(Data)) | grepl("SS52_naive", colnames(Data)))] <- "B6" 
species[which(grepl("SS19_naive", colnames(Data)) | grepl("SS25_naive", colnames(Data)))] <- "CAST" 

# Calculated tSNE
tsne <- Rtsne(t(log10(Data + 1)), perplexity = 10)
tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], species = species)

# Plot tSNE
ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 4, aes(bg = species), pch = 21) +
  scale_fill_manual(values = c("black", "coral4", "tan3")) + 
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

# Caluclate PCA
pca <- prcomp(t(log10(Data + 1)))
pca.df <- data.frame(pca.1 = pca$x[,1], pca.2 = pca$x[,2], species = species)

# Plot PCA
ggplot(data = pca.df, aes(pca.1, pca.2)) + 
  geom_point(size = 4, aes(bg = species), pch = 21) +
  scale_fill_manual(values = c("black", "coral4", "tan3")) + 
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

# Fig 1C

# Load BASiCS Test file comparing B6 young naive and CAST young naive
load("/path/to/Test_logFC2_file/")
B6_CAST <- Test_logFC2$Table

### Find genes that are differentially expressed
# B6 specific
B6 <- B6_CAST[B6_CAST$ResultDiffExp == "SS51_naive_SS52_naive+",]
B6.spec <- as.character(B6$GeneNames)

# CAST specific
CAST <- B6_CAST[B6_CAST$ResultDiffExp == "SS19_naive_SS25_naive+",]
CAST.spec <- as.character(CAST$GeneNames)

### Remove genes that are caused by mapping artifacts - See quality control
# B6 and CAST samples mapped to CAST genome
load("/path/to/Test_logFC2_file_CAST/")
DE.CAST.map <- Test_logFC2$Table
# B6 and CAST samples mapped to B6 genome
load("/path/to/Test_logFC2_file_B6/")
DE.B6.map <- Test_logFC2$Table

# Collect genes that are differenttially mapped between the two Genomes - for B6 and CAST
genes <- c(DE.B6.map$GeneNames[DE.B6.map$ResultDiffExp != "NoDiff" & DE.B6.map$ResultDiffExp != "ExcludedByUser"], 
           DE.CAST.map$GeneNames[DE.CAST.map$ResultDiffExp != "NoDiff" & DE.CAST.map$ResultDiffExp != "ExcludedByUser"])
genes <- unique(genes)

# Remove these genes from B6 specifc genes ... 
B6.spec <- B6.spec[-which(B6$GeneNames %in% genes)]

# ... and CAST specifc genes
CAST.spec <- CAST.spec[-which(CAST$GeneNames %in% genes)]

# Collect cells from B6 young naive and CAST young naive for plotting
Data <- exp.data[,c(which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data))),
                    which(grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data))))]

### Plot specific genes
# For visualization purposes sample 30 genes and 30 cells 
B6 <- colnames(Data)[which(grepl("SS51", colnames(Data)) | grepl("SS52", colnames(Data)))]
B6 <- B6[sample(1:length(B6), 30)]
B6.genes <- B6.spec[sample(1:length(B6.spec), 30)]
CAST <- colnames(Data)[which(grepl("SS19", colnames(Data)) | grepl("SS25", colnames(Data)))]
CAST <- CAST[sample(1:length(CAST), 30)]
CAST.genes <- CAST.spec[sample(1:length(CAST.spec), 30)]
cells <- c(B6, CAST)

# Plot species specifc heatmaps
pheatmap(log10(Data[B6.genes,cells] + 1), cellwidth = 2, cellheight = 2, fontsize = 8, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(30))
pheatmap(log10(Data[CAST.genes,cells] + 1), cellwidth = 2, cellheight = 2, fontsize = 8, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(30))


### Figure 1D
### Plot vioplots for selected genes 
genes <- c("Klf13","Ets1", "Ly6e", "Rhoh", "Il21r", "Tmed7",
           "Gm2a", "Serinc3", "Ctse", "Arl5a", "Mccc2", "Hdgfrp3")

gene.ids <- rownames(genenames)[match(genes, as.character(genenames[,2]))]

### Look at variability of the selected genes - in increasing order
B6_CAST[match(gene.ids[1:6], B6_CAST$GeneNames),11]
B6_CAST[match(genes[7:12], as.character(genenames[as.character(B6_CAST$GeneNames),2])),10]


# Plot violin plots for B6 and CAST
violins(x = as.data.frame(log10(t(exp.data[rev(gene.ids),which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))])+ 1)),
        horizontal = TRUE, drawRect = FALSE, deciles = FALSE, connect = FALSE, ylim = c(0,6),
        main = gene)

violins(x = as.data.frame(log10(t(exp.data[rev(gene.ids),which(grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)))])+ 1)),
        horizontal = TRUE, drawRect = FALSE, deciles = FALSE, connect = FALSE, ylim = c(0,6),
        main = gene)

