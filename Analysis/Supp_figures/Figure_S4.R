################################################################################
#### Code to generate Figure S4 ################################################
################################################################################
library(pheatmap)
library(ggplot2)

exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

### Fig S4A
Data <- exp.data[,c(which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data))),  
                    which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data))))]

pca <- prcomp(log10(t(Data) + 1))
pca.df <- data.frame(pca.1 = pca$x[,1], pca.2 = pca$x[,2], activation = c(rep("+", 53), rep("-", 93)))

ggplot(data = pca.df, aes(pca.1, pca.2)) + 
  geom_point(size = 8, aes(colour = activation),shape = "-") +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

### Fig S4B
# Load a BASiCS Test dataframe tested on logFC=2 (Test_logFC2) between naive and activated cells
DE <- Test_logFC2$Table
genes.naive <- as.character(DE$GeneNames[DE$ResultDiffExp == "SS51_naive_SS52_naive+"])
genes.active <- as.character(DE$GeneNames[DE$ResultDiffExp == "SS51_active_SS52_active+"])

# Plot heatmap ordered by mean expression
pheatmap(log10(Data[c(as.character(genenames[genes.naive[order(rowMeans(Data[as.character(genenames[genes.naive,2]),]), decreasing = TRUE)],2]), 
                      as.character(genenames[genes.active[order(rowMeans(Data[as.character(genenames[genes.active,2]),]), decreasing = TRUE)],2])),] + 1), cluster_cols = FALSE, cluster_rows = FALSE,
         col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), show_rownames = FALSE, show_colnames = FALSE)

### Fig S4C
# Fing genes with high variability in the naive state
# Load BASiCS DV object (Test_logFC0)
DV <- Test_logFC0$Table
DV <- DV[DV$ResultDiffExp == "NoDiff",]

DV.genes <- DV$GeneNames[DV$ResultDiffOverDisp == "SS51_naive_SS52_naive+"]
DV.genes <- as.character(genenames[DV.genes,2])

# Collect naive and activated cells
naive.cells <- exp.data[DV.genes,which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]
active.cells <- exp.data[DV.genes,which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]

# Downsample cells
naive.cells <- naive.cells[,sample(1:93, 50)]
active.cells <- active.cells[,sample(1:53, 50)]

# Only consider genes that are expressed in both conditions
naive.cells <- naive.cells[rowMeans(naive.cells) > 50 & rowMeans(active.cells) > 50,]
active.cells <- active.cells[rownames(naive.cells),]

# Compute pearson correlations between genes
cors.naive <- cor(t(naive.cells), method = "pearson")
cors.active <- cor(t(active.cells), method = "pearson")

# Set diagonal to 0
diag(cors.active) <- 0
diag(cors.naive) <- 0

# Collect highly correlating genes
r.a <- rownames(cors.active)[apply(cors.active, 1, function(n){ifelse(length(which(n > 0.8)) > 0, TRUE, FALSE)})]
r.n <- rownames(cors.naive)[apply(cors.naive, 1, function(n){ifelse(length(which(n > 0.8)) > 0, TRUE, FALSE)})]

# Plot most correlated genes
pheatmap(log10(naive.cells[r.n,] + 1), col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(naive.cells[r.n,]))/log10(max(naive.cells[r.n,]))*100, 0)], border_color = NA, cellwidth = 4, cellheight = 4, fontsize = 5, clustering_distance_rows = as.dist(1 - cor(log10(t(naive.cells[r.n,]) + 1), method = "spearman")))
pheatmap(log10(active.cells[r.n,] + 1), col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(naive.cells[r.n,]))/log10(max(naive.cells[r.n,]))*100, 0)], border_color = NA, cellwidth = 4, cellheight = 4, fontsize = 5, clustering_distance_rows = as.dist(1 - cor(log10(t(naive.cells[r.n,]) + 1), method = "spearman")))

### Fig S4D
# Calculate the fraction of cells expressing up-regulated genes in:
# active cells
active.cells.active <- exp.data[which(!is.na(match(rownames(exp.data), as.character(genenames[genes.active,2])))),
                                which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]

# naive cells
naive.cells.active <- exp.data[which(!is.na(match(rownames(exp.data), as.character(genenames[genes.active,2])))),
                               which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]

# down-regulated genes in:
# naive cells
naive.cells.naive <- exp.data[which(!is.na(match(rownames(exp.data), as.character(genenames[genes.naive,2])))),
                              which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]

# active cells
active.cells.naive <- exp.data[which(!is.na(match(rownames(exp.data), as.character(genenames[genes.naive,2])))),
                               which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]

# Downsample to 600 genes
plot(hist(apply(naive.cells.naive[sample(1:nrow(naive.cells.naive), 600),], 1, function(n){length(which(n > 0))/ncol(naive.cells.naive)}), breaks = 50), ylim = c(0,250), xlim = c(0,1))
abline(v = median(apply(naive.cells.naive, 1, function(n){length(which(n > 0))/ncol(naive.cells.naive)})))

plot(hist(apply(active.cells.active[sample(1:nrow(active.cells.active), 600),], 1, function(n){length(which(n > 0))/ncol(active.cells.active)}), breaks = 50), ylim = c(0,250), xlim = c(0,1))
abline(v = median(apply(active.cells.active, 1, function(n){length(which(n > 0))/ncol(active.cells.active)})))

plot(hist(apply(active.cells.naive[sample(1:nrow(active.cells.naive), 600),], 1, function(n){length(which(n > 0))/ncol(active.cells.naive)}), breaks = 50), ylim = c(0,250), xlim = c(0,1))
abline(v = median(apply(active.cells.naive, 1, function(n){length(which(n > 0))/ncol(active.cells.naive)})))

plot(hist(apply(naive.cells.active[sample(1:nrow(naive.cells.active), 600),], 1, function(n){length(which(n > 0))/ncol(naive.cells.active)}), breaks = 50), ylim = c(0,250), xlim = c(0,1))
abline(v = median(apply(naive.cells.active, 1, function(n){length(which(n > 0))/ncol(naive.cells.active)})))

### Fig S4E
# See Fig 3C

### Fig S4F+G
# Use scLVM to correct for translation effect
input <- read.table("Data/raw_data.txt")

dat <- input[which(grepl("ENS", rownames(input))),which(grepl("SS51_active", colnames(input)) | grepl("SS52_active", colnames(input)))]
ERCC <- input[grepl("ERCC-", rownames(input)),colnames(dat)]
dat <- dat[rowSums(dat) > 0,]

# Read in GO categories for genes up-regulated during activation
GOs <- read.table("/path/to/GO_categories.txt", header = TRUE, sep = "\t")

# Select translation related GO categories
genes <- as.character(GOs[1:6, 6])
genes <- as.character(unlist(sapply(genes, function(n){unlist(strsplit(n, ", "))})))
genes <- unique(genes)

# Prepare for scLVM correction
Y = t(log10(dat+1)) 
geneID = rownames(dat) 
techNoise <- scLVM::fitTechnicalNoise(nCountsEndo = dat, nCountsERCC = ERCC, use_ERCC = TRUE, fit_type = "counts", plot = TRUE)
tech_noise = as.vector(techNoise$techNoiseLog) #technical noise

# Highly variable genes
is_het = scLVM::getVariableGenes(dat, techNoise$fit, method = "fit", 
                                 threshold = 0.1, fit_type="counts")

# Initiate scLVM
sclvm = new("scLVM")
sclvm = scLVM::init(sclvm,Y=Y,tech_noise = tech_noise)

# Find numbers of latent factors
TranslationARD = fitFactor(sclvm,geneSet = genes, k=20,use_ard = TRUE)

plot(seq(1, length(TranslationARD$X_ard)), TranslationARD$X_ard, xlab = '# Factor', ylab = 'Variance explained')
title('Variance explained by latent factors')

# Use one factor
Translation = fitFactor(sclvm,geneSet = genes, k=1,)

Kt = Translation$K
Xt = Translation$X

image(Kt,xaxt = "n", yaxt = "n", xlab = 'cells', ylab = 'cells')
title('Similarity matrix based on cell cycle')

idx_het = which(is_het)

# fit the model for variable genes
sclvm = varianceDecomposition(sclvm, K=Kt, idx = idx_het)

results_var = getVarianceComponents(sclvm)

var_filtered = results_var$var[results_var$conv,] 

# filter out genes for which vd has not converged
head(var_filtered)

# correct expression values
Ycorr = getCorrectedExpression(sclvm)
dim(Ycorr)

# Look at PCA for uncorrected and corrected values
### Fig S4F
pca <- prcomp(log10(Ycorr + 1))
plot(pca$x[,1], pca$x[,2], pch = 16, pch = "|")

### Fig S4G
# normalised data
dat.norm <- exp.data[which(grepl("ENS", rownames(exp.data))),which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]
pca <- prcomp(t(log10(dat.norm + 1)))
plot(pca$x[,1], pca$x[,2], pch = 16, pch = "|")

