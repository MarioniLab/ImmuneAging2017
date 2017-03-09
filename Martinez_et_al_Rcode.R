###############################################################
#### Differential expression and variability analysis #########
###############################################################

library(BASiCS)

#### Generation of MCMC objects ###############################

# Read in data
input <- read.table("/path/to/file", header = TRUE, sep = "\t")

# Split biological genes and spike-ins
bio.g <- input[which(grepl("ENSMUS", rownames(input))),]
ERCC <- input[which(grepl("ERCC", rownames(input))),]

# Read in ERCC concentrations
cur.url <- "https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt"
download.file(url = cur.url, destfile = basename(cur.url), mode = "w")
ERCC.conc <- read.table("path/to/cms_095046.txt", header=TRUE, sep = "\t", fill = TRUE)

# Calculate total number of ERCC molecules per well
ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/50000 # dilution factor
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

# Pre-scale reads to adjust for differences in library size
rpm <- (bio.g/colSums(bio.g))*1000000

# Filtering on genes
cell.count <- apply(rpm, 1, function(n){length(which(n > 20))})
bio.g.1 <- bio.g[names(which(cell.count > 2)),]
ERCC.1 <- ERCC[rowSums(ERCC) > 0, ]

# Forming the Data object
Counts <- rbind(bio.g.1, ERCC.1)
Tech <- c(rep(FALSE, nrow(bio.g.1)), rep(TRUE, nrow(ERCC.1)))
SpikeInput <- ERCC.num.final[rownames(ERCC.1),1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data <- newBASiCS_Data(as.matrix(Counts), Tech, SpikeInput.1)

# Starting the MCMC - log-normal prior on delta if genes have zero counts
MCMC_Output <- BASiCS_MCMC(Data, N = 40000, Thin = 20, 
                           Burn = 20000, PrintProgress = TRUE, 
                           PriorDelta = "log-normal")

# Check for Chain convergence
plot(MCMC_Summary, Param = "delta") #### cell-to-cell heterogeneity regarding to mRNA content
plot(MCMC_Summary, Param = "s") #### possible evidence of cell-to-cell differeces in capture efficiency 
plot(MCMC_Output, Param = "delta", Gene = 1)
plot(MCMC_Output, Param = "mu", Gene = 1)
plot(MCMC_Output, Param = "s", Cell = 20)

#### Retrieving normalized counts ###############################
DenoisedCounts = BASiCS_DenoisedCounts(Data = Data, Chain = MCMC_Output)
DenoisedRates = BASiCS_DenoisedRates(Data = Data, Chain = MCMC_Output)

#### Differential expression and variability testing ############
### For testing 2 MCMC objects are needed
MCMC.A <- MCMC_output.1
MCMC.B <- MCMC_output.2
Data.A <- Data.1
Data.B <- Data.2

# Function for correcting differences in population wide RNA content
OffSetCorrection <- function(MCMC1, MCMC2){
  median(rowSums(MCMC1@mu)/rowSums(MCMC2@mu)) 
}

# Forming the Data_D object
CountsTest = Data.B@Counts
CountsRef = Data.A@Counts
Tech = as.logical(c(rep("FALSE", length(which(grepl("ENS", rownames(CountsRef))))), 
                    rep("TRUE", length(which(grepl("ERCC", rownames(CountsRef)))))))

SpikeInput <- ERCC.num.final[rownames(CountsRef)[which(grepl("ERCC", rownames(CountsRef)))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data_DV <- newBASiCS_D_Data(CountsTest = CountsTest, GeneNames = rownames(CountsTest),
                            CountsRef = CountsRef,
                            Tech = Tech,
                            SpikeInputTest = SpikeInput.1[,2],
                            SpikeInputRef = SpikeInput.1[,2]
                            )

# Offset correction
OffSet = OffSetCorrection(MCMC.A, MCMC.B)

# Forming the MCMC_D object
MCMC_DV <- newBASiCS_D_Chain(muTest = MCMC.B@mu,
                             muRef = MCMC.A@mu / OffSet,
                             deltaTest = MCMC.B@delta,
                             delta = MCMC.A@delta,
                             phi = cbind(MCMC.B@phi, MCMC.A@phi * OffSet), s = cbind(MCMC.B@s, MCMC.A@s),
                             nu = cbind(MCMC.B@nu, MCMC.A@nu),
                             thetaTest = MCMC.B@theta,
                             thetaRef = MCMC.A@theta)

### Genes to exclude from testing
excl.genes <- colnames(MCMC_DV@muTest)[which(colMeans(MCMC_DV@muTest) < 50 & 
                                               colMeans(MCMC_DV@muRef) < 50)]
sel <- rep(TRUE, length(genes))
sel[match(excl.genes, genes)] <- FALSE

# Testing for Differential expression 
Test_logFC2 <- BASiCS_D_TestDE(Data_DV, MCMC_DV, GeneNames = genes,
                        EpsilonM = 2, EpsilonD = 0.4, GenesSelect = sel, 
                        EFDR_M = 0.05, EFDR_D = 0.05,
                        OrderVariable = "GeneNames", GroupLabelRef = cur_names[1], 
                        GroupLabelTest = cur_names[2])

# Testing for Differential variability
Test_logFC0 <- BASiCS_D_TestDE(Data_DV, MCMC_DV, GeneNames = genes,
                        EpsilonM = 0, EpsilonD = 0.4, GenesSelect = sel, 
                        EFDR_M = 0.05, EFDR_D = 0.05,
                        OrderVariable = "GeneNames", GroupLabelRef = cur_names[1], 
                        GroupLabelTest = cur_names[2])

###############################################################
#### Code to reproduce Figures ################################
###############################################################

library(Rtsne)
library(pheatmap)
library(ggplot2)
library(BASiCS)
library(RColorBrewer)
library(caroline)

#### Figure 1 #################################################
exp.data <- read.table("/path/to/normalized_data/")
load("/Users/eling01/Google Drive/scMouse_Immun/Analysis/Variance/DV_DE_all_B6CAST.RData")

### Figure 1B

# Select cell from B6 young naive and CAST young naive
Data <- exp.data[,which(grepl("B6_young_naive", colnames(exp.data)) |
                           grepl("CAST_young_naive", colnames(exp.data)))]

# Color factor for plotting
species <- vector(length = ncol(Data))
species[which(grepl("B6_young_naive", colnames(Data)))] <- "B6" 
species[which(grepl("CAST_young_naive", colnames(Data)))] <- "CAST" 

# Calculate tSNE
tsne <- Rtsne(t(log10(Data + 1)), perplexity = 40)
tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], species = species)

# Plot tSNE
ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 4, aes(bg = species), pch = 21) +
  scale_fill_manual(values = c("black", "coral4", "tan3")) + 
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

### Figure 1C

# DE test between B6 young naive and CAST young naive - see above
B6_CAST <- Test_logFC2_B6_CAST

# Find genes that are differentially expressed
# B6 specific
DE.B6_CAST <- B6_CAST[B6_CAST$ResultDiffExp == "B6_young_naive+",]
B6.spec <- as.character(genenames[as.character(DE.B6_CAST$GeneNames),2])

# CAST specific
DE.CAST_B6 <- B6_CAST[B6_CAST$ResultDiffExp == "CAST_young_naive+",]
CAST.spec <- as.character(genenames[as.character(DE.CAST_B6$GeneNames),2])

# Heatmaps of 30 randomly selected species-specifc genes and cells
B6 <- colnames(Data)[which(grepl("B6_young_naive", colnames(Data)))]
B6 <- B6[sample(1:length(B6), 30)]
B6.genes <- B6.spec[sample(1:length(B6.spec), 30)]

CAST <- colnames(Data)[which(grepl("CAST_young_naive", colnames(Data)))]
CAST <- CAST[sample(1:length(CAST), 30)]
CAST.genes <- CAST.spec[sample(1:length(CAST.spec), 30)]

cells <- c(B6, CAST)

pheatmap(log10(Data[B6.genes,cells] + 1), cellwidth = 2, cellheight = 2, fontsize = 8, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(30, 60))

pheatmap(log10(Data[CAST.genes,cells] + 1), cellwidth = 2, cellheight = 2, fontsize = 8, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(30, 60))

### Figure 1D

# Plot violinplots of hand picked genes with different levels of variabilty
genes <- c("Klf13","Ets1", "Ly6e", "Rhoh", "Il21r", "Tmed7",
           "Gm2a", "Serinc3", "Ctse", "Arl5a", "Mccc2", "Hdgfrp3")

# Look at variability measures
B6_CAST[match(genes[1:6], as.character(genenames[as.character(B6_CAST$GeneNames),2])),11]
B6_CAST[match(genes[7:12], as.character(genenames[as.character(B6_CAST$GeneNames),2])),10]

violins(x = as.data.frame(log10(t(exp.data[rev(genes),which(grepl("B6_young_naive", colnames(exp.data)) )])+ 1)),
        horizontal = TRUE, drawRect = FALSE, deciles = FALSE, connect = FALSE, ylim = c(0,6))

violins(x = as.data.frame(log10(t(exp.data[rev(genes),which(grepl("CAST_young_naive", colnames(exp.data)))])+ 1)),
        horizontal = TRUE, drawRect = FALSE, deciles = FALSE, connect = FALSE, ylim = c(0,6))

#### Figure 2 #################################################

### Figure 2A

# Select naive and activated cells from young B6 animals
Data <- exp.data[,c(which(grepl("B6_young_naive", colnames(exp.data))),
                       which(grepl("B6_young_active", colnames(exp.data))))]

# Compute tSNE
tsne <- Rtsne(t(log10(Data + 1)), perplexity = 10)
tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], activation = c(rep("-", 93), rep("+", 53)))

ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 8, aes(colour = activation),shape = "-") +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

### Figure 2B
# Analysis of DE between naive and activated cells from young B6
# For this analysis a DE and DV Test object is needed - Test_logFC2 and Test_logFC0

Test_logFC2$cols <- ifelse(Test_logFC2$ResultDiffExp == "B6_young_active+", "Active", ifelse(Test_logFC2$ResultDiffExp == "B6_young_naive+", "Naive", ifelse(Test_logFC0$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("dark red", "black", "dark blue", "grey")) + 
  theme_minimal(base_size = 20)

### Figure 2C
# For this analysis only consider genes that do not change in mean expression - Test_logFC0

var.gene <- data.frame(exp = log10(c(Test_logFC0$OverDispRef[which(Test_logFC0$ResultDiffExp == "NoDiff")], Test_logFC0$OverDispTest[which(Test_logFC0$ResultDiffExp == "NoDiff")])),
                       cond = factor(c(rep("Naive", length(Test_logFC0$OverDispRef[which(Test_logFC0$ResultDiffExp == "NoDiff")])),
                                       rep("Active", length(Test_logFC0$OverDispTest[which(Test_logFC0$ResultDiffExp == "NoDiff")]))), levels = c("Naive", "Active"))
)

ggplot(var.gene, aes(cond, exp)) + 
  geom_boxplot(aes(fill = cond), outlier.shape = NA) +
  theme_minimal(base_size = 30) +
  ylim(0,2.5) +
  xlab("State") + ylab("log10[Over-dispersion]") + 
  scale_fill_manual(values = c("red", "blue"))

# Test for difference in overall varibility measures
wilcox.test(log10(Test_logFC0$OverDispRef[which(Test_logFC0$ResultDiffExp == "NoDiff")]),
            log10(Test_logFC0$ExpTest[which(Test_logFC0$ResultDiffExp == "NoDiff")]))

### Figure 2D
# Plot the frequency of gene expression for genes up-regulated and down-regulated during activation

active.cells <- Data[as.character(genenames[Test_logFC2$GeneNames[Test_logFC2$ResultDiffExp == "B6_young_active+"],2]),
                         which(grepl("B6_young_active", colnames(exp.data)))]
naive.cells <- Data[as.character(genenames[Test_logFC2$GeneNames[Test_logFC2$ResultDiffExp == "B6_young_naive+"],2]),
                        which(grepl("B6_young_naive", colnames(exp.data)))]

### Histogram of cells detected
plot(hist(apply(naive.cells[sample(1:nrow(naive.cells), 600),], 1, function(n){length(which(n > 0))/ncol(naive.cells)}), breaks = 50), xlim = c(0,1),ylim = c(0,50))
abline(v = median(apply(naive.cells[sample(1:nrow(naive.cells), 600),], 1, function(n){length(which(n > 0))/ncol(naive.cells)})))

plot(hist(apply(active.cells[sample(1:nrow(active.cells), 600),], 1, function(n){length(which(n > 0))/ncol(active.cells)}), breaks = 50), xlim = c(0,1), ylim = c(0,50))
abline(v = median(apply(active.cells[sample(1:nrow(active.cells), 600),], 1, function(n){length(which(n > 0))/ncol(active.cells)})))

### Figure 2EFG
# Plot single-cell transcript counts for selected genes

# Downsample to 50 cells for visualization purposes
cells.naive <- colnames(Data)[which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]
cells.active <- colnames(Data)[which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]

# Genes are: Sell, Il2ra and Eif1
gene <-"ENSMUSG00000026581"

barplot(log10(Data[gene,cells.naive[sample(1:length(cells.naive), 50)]] + 1), ylim = c(0,6))
abline(h = log10(Test_logFC2$ExpRef[match(gene, Test_logFC2$GeneNames)]))

barplot(log10(Data[gene,cells.active[sample(1:length(cells.active), 50)]] + 1), ylim = c(0,6))
abline(h = log10(Test_logFC2$ExpTest[match(gene, Test_logFC2$GeneNames)]))

#### Figure 3 #################################################

### Figure 3A
# tSNE of naive and activated cells from CAST and B6
Data <- exp.data[,which(grepl("B6_young_naive", colnames(exp.data)) |
                                 grepl("B6_young_active", colnames(exp.data)) | 
                                 grepl("CAST_young_naive", colnames(exp.data)) | 
                                 grepl("CAST_young_active", colnames(exp.data)))]

tsne <- Rtsne(t(log10(Data + 1)), perplexity = 10)

state <- vector(length = ncol(Data))
state[which(grepl("B6_young_naive", colnames(Data)))] <- "B6_naive" 
state[which(grepl("B6_young_active", colnames(Data)))] <- "B6_active"
state[which(grepl("CAST_young_naive", colnames(Data)))] <- "CAST_naive" 
state[which(grepl("CAST_young_active", colnames(Data)))] <- "CAST_active" 

tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], State = state)

ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 10, aes(col = state), shape = "-") +
  scale_colour_manual(values = c("black", "coral4", "tan3", "red")) + 
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 

### Figure 3B

# Detection of shared activation program - three Test objects are needed: 
# B6 young naive vs active: Test_logFC2_B6_nva
# CAST young naive vs active: Test_logFC2_CAST_nva
# Young active B6 vs CAST: Test_logFC2_a_B6vCAST
DE.B6.young.active <- Test_logFC2_B6_nva[Test_logFC2_B6_nva$ResultDiffExp == "B6_young_active+",]

DE.CAST.young.active <- Test_logFC2_CAST_nva[Test_logFC2_CAST_nva$ResultDiffExp == "CAST_young_active+",]

B6_CAST_active.CAST <- Test_logFC2_a_B6vCAST[Test_logFC2_a_B6vCAST$ResultDiffExp == "CAST_young_active+",]
B6_CAST_active.B6 <- Test_logFC2_a_B6vCAST[Test_logFC2_a_B6vCAST$ResultDiffExp == "B6_young_active+",]

# The shared activation program are genes that are similarly up-regulated in B6 and CAST
shared.a <- match(as.character(DE.B6.young.active$GeneNames), as.character(DE.CAST.young.active$GeneNames))
shared.a <- as.character(DE.CAST.young.active$GeneNames)[shared.a[!is.na(shared.a)]]

# To find species specifc up-regulated genes we test a specific subset of all genes
# Only the 1208 genes that are up-regulated in either B6 or CAST
# Therefore we want to find genes that are up-regulated but DE between the two species
spec.test <- Test_log2FC_up_B6vs_CAST
B6.spec.a <- as.character(spec.test$GeneNames[spec.test$ResultDiffExp == "B6_young_active+"])

CAST.spec.a <- as.character(spec.test$GeneNames[spec.test$ResultDiffExp == "CAST_young_active+"])

# Heatmap of shared activation process - proportionally downsampled genes - cells downsampled to 30
Data <- exp.data[,c(which(grepl("B6_young_naive", colnames(exp.data)))[sample(1:length(which(grepl("B6_young_naive", colnames(exp.data)))), 30)],
                    which(grepl("CAST_young_naive", colnames(exp.data)))[sample(1:length(which(grepl("CAST_young_naive", colnames(exp.data)))), 30)],
                    which(grepl("B6_young_active", colnames(exp.data)))[sample(1:length(which(grepl("B6_young_active", colnames(exp.data)))), 30)],
                    which(grepl("CAST_young_active", colnames(exp.data)))[sample(1:length(which(grepl("CAST_young_active", colnames(exp.data)))), 30)]
)]

genes.shared <- shared.a[sample(1:length(shared.a), 100)]
genes.B6 <- B6.spec.a[sample(1:length(B6.spec.a), 38)]
genes.CAST <- CAST.spec.a[sample(1:length(CAST.spec.a), 34)]

pheatmap(log10(Data[genes.shared[order(rowMeans(Data[genes.shared,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.shared,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))

pheatmap(log10(Data[genes.B6[order(rowMeans(Data[genes.B6,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.B6,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))

pheatmap(log10(Data[genes.CAST[order(rowMeans(Data[genes.CAST,]), decreasing = TRUE)],] + 1), cluster_rows = FALSE, fontsize_row = 2, fontsize_col = 4, cellwidth = 3, cellheight = 3, colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100)[1:round(log10(max(Data[genes.CAST,]))/6*100, 0)], border_color = NA, cluster_cols = FALSE, gaps_col = c(30,60,90))

### Look at frequencies - downsampled genes to 70
genes.shared <- shared.a[sample(1:length(shared.a), 70)]
genes.B6 <- B6.spec.a[sample(1:length(B6.spec.a), 70)]
genes.CAST <- CAST.spec.a[sample(1:length(CAST.spec.a), 70)]

shared.cells <- Data[genes.shared, which(grepl("active", colnames(Data)))]
B6.cells <- exp.data[genes.B6, which(grepl("active", colnames(Data)))]
CAST.cells <- exp.data[genes.CAST, which(grepl("active", colnames(Data)))]


pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4D_shared.pdf"),  width = 7, height = 7, useDingbats = TRUE)
plot(hist(apply(shared.cells[sample(1:nrow(shared.cells), 70),], 1, function(n){length(which(n > 0))/ncol(shared.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(shared.cells, 1, function(n){length(which(n > 0))/ncol(shared.cells)})))
dev.off()

pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4D_B6.pdf"),  width = 7, height = 7, useDingbats = TRUE)
plot(hist(apply(B6.cells[sample(1:nrow(B6.cells), 70),], 1, function(n){length(which(n > 0))/ncol(B6.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(B6.cells, 1, function(n){length(which(n > 0))/ncol(B6.cells)})))
dev.off()

pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4D_CAST.pdf"),  width = 7, height = 7, useDingbats = TRUE)
plot(hist(apply(CAST.cells[sample(1:nrow(CAST.cells), 70),], 1, function(n){length(which(n > 0))/ncol(CAST.cells)}), breaks = 20), xlim = c(0,1),ylim = c(0,30))
abline(v = median(apply(CAST.cells, 1, function(n){length(which(n > 0))/ncol(CAST.cells)})))
dev.off()


### D) GO categories
write.table(as.data.frame(B6_CAST_active$GeneNames[B6_CAST_active$ResultDiffExp != "ExcludedByUser"]), "/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Fig4D_background.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(as.character(genenames[match(shared.a, genenames[,2]),1])), "/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Fig4D_shared.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(as.character(genenames[match(B6.spec.a, genenames[,2]),1])), "/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Fig4D_B6.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(as.character(genenames[match(CAST.spec.a, genenames[,2]),1])), "/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Fig4D_CAST.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

results <- read.table("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/GO_shared_genes.txt", header = TRUE, sep = '\t')
results <- results[results$Bonferroni < 0.8,]
pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fi4D_shared.pdf"),  width = 10, height = 10, useDingbats = TRUE)
par(mar = c(16,4,1,4))
barplot(height = -log10(results$Bonferroni), col = "white", names.arg = results$Term, las = 2, cex.names = 0.5, ylim = c(0,5))
abline(a = 1, b = 0, col = "red", lwd = 2)
dev.off()

results.B6 <- read.table("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/GO_B6_genes.txt", header = TRUE, sep = '\t')
pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4D_B6.pdf"),  width = 10, height = 10, useDingbats = TRUE)
par(mar = c(16,4,1,4))
barplot(height = -log10(results.B6$Bonferroni), col = "white", names.arg = results.B6$Term, las = 2, cex.names = 0.5, ylim = c(0,5))
abline(a = 1, b = 0, col = "red", lwd = 2)
dev.off()

results.CAST <- read.table("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/GO_CAST_genes.txt", header = TRUE, sep = '\t')
pdf(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4D_CAST.pdf"),  width = 10, height = 10, useDingbats = TRUE)
par(mar = c(16,4,1,4))
barplot(height = -log10(results.CAST$Bonferroni), col = "white", names.arg = results.CAST$Term, las = 2, cex.names = 0.5, ylim = c(0,5))
abline(a = 1, b = 0, col = "red", lwd = 2)
dev.off()
