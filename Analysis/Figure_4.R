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

### Figure 4
# Core activation process is robust during ageing
# Read in shared activation genes - supplementary table 3
library(xlsx)
shared.a <- read.xlsx("Data/S3.xlsx", sheetIndex = 1)
shared.a  <- as.character(shared.a $Gene.ID)

# Read in the BASiCS test objects tested on logFC0: B6 active young vs old and CAST active young vs old
### Figure 4A
# Evaluate changes in B6
load("/Users/nils/Downloads//DV_DE_all_B6CAST_logFC0.RData")
DV.B6.old <- TestLogFC0.B6
# Select the rows with shared activation genes
DV.B6.old <- DV.B6.old[match(shared.a, as.character(DV.B6.old$GeneNames)),]

# Order based on mean expression
DV.B6.old <- DV.B6.old[order(DV.B6.old$ExpOverall, decreasing = TRUE),]

# Plot Log2FC in mean values
plot(DV.B6.old$ExpLogFC, pch=16, col="grey", ylim = c(-6, 6))
abline(h = 0, lwd = 3, col = "red")

# Calculated percentages
length(which(DV.B6.old$ExpLogFC > 0))/length(DV.B6.old$ExpLogFC)
length(which(DV.B6.old$ExpLogFC < 0))/length(DV.B6.old$ExpLogFC)

# Plot Log2FC in over-dispersion values - for genes that do not change in mean expression
plot(DV.B6.old$OverDispLogFC[DV.B6.old$ResultDiffExp == "NoDiff"], pch=16, ylim = c(-2, 2), col = "grey")
abline(h = 0, lwd = 3, col = "red")

# Calculated percentages
length(which(DV.B6.old$OverDispLogFC[DV.B6.old$ResultDiffExp == "NoDiff"] > 0))/length(DV.B6.old$OverDispLogFC[DV.B6.old$ResultDiffExp == "NoDiff"])
length(which(DV.B6.old$OverDispLogFC[DV.B6.old$ResultDiffExp == "NoDiff"] < 0))/length(DV.B6.old$OverDispLogFC[DV.B6.old$ResultDiffExp == "NoDiff"])


### Figure 4B
# Evaluate changes in CAST
DV.CAST.old <- TestLogFC0.CAST
# Select the rows with shared activation genes
DV.CAST.old <- DV.CAST.old[match(shared.a, as.character(DV.CAST.old$GeneNames)),]

# Order based on mean expression
DV.CAST.old <- DV.CAST.old[order(DV.CAST.old$ExpOverall, decreasing = TRUE),]

# Plot Log2FC in mean values
plot(DV.CAST.old$ExpLogFC, pch=16, col="grey", ylim = c(-6, 6))
abline(h = 0, lwd = 3, col = "red")

# Calculated percentages
length(which(DV.CAST.old$ExpLogFC > 0))/length(DV.CAST.old$ExpLogFC)
length(which(DV.CAST.old$ExpLogFC < 0))/length(DV.CAST.old$ExpLogFC)

# Plot Log2FC in over-dispersion values - for genes that do not change in mean expression
plot(DV.CAST.old$OverDispLogFC[DV.CAST.old$ResultDiffExp == "NoDiff"], pch=16, ylim = c(-2, 2), col = "grey")
abline(h = 0, lwd = 3, col = "red")

# Calculated percentages
length(which(DV.CAST.old$OverDispLogFC[DV.CAST.old$ResultDiffExp == "NoDiff"] > 0))/length(DV.CAST.old$OverDispLogFC[DV.CAST.old$ResultDiffExp == "NoDiff"])
length(which(DV.CAST.old$OverDispLogFC[DV.CAST.old$ResultDiffExp == "NoDiff"] < 0))/length(DV.CAST.old$OverDispLogFC[DV.CAST.old$ResultDiffExp == "NoDiff"])

### Figure 4C
### Look at expression ratio
B6.young <- exp.data[,which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]
B6.old <- exp.data[,which(grepl("SS63_active", colnames(exp.data)) | grepl("SS64_active", colnames(exp.data)))]
CAST.young <- exp.data[,which(grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]
CAST.old <- exp.data[,which(grepl("SS53_active", colnames(exp.data)) | grepl("SS54_active", colnames(exp.data)))]

# Plot density of expression ratios in the 4 dfferent conditions
plot(density(apply(B6.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.young)})), ylim = c(0,3), xlim = c(0,1))
abline(v = median(apply(B6.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.young)})))

plot(density(apply(B6.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.old)})), ylim = c(0,3), xlim = c(0,1))
abline(v = median(apply(B6.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.old)})))

plot(density(apply(CAST.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.young)})), ylim = c(0,3), xlim = c(0,1))
abline(v = median(apply(CAST.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.young)})))

plot(density(apply(CAST.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.old)})), ylim = c(0,3), xlim = c(0,1))
abline(v = median(apply(CAST.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.old)})))

# Plot the ratios as scatter plots: old vs young - colour by multiple testing corrected p-values
# B6 

### Do statistical testing
stats <- list()
for(i in 1:length(shared.a)){
  bin.t <- binom.test(x = length(which(B6.young[shared.a[i],] > 0)), n = ncol(B6.young), p = length(which(B6.old[shared.a[i],] > 0))/ncol(B6.old))
  stats[[i]] <- bin.t
}

p.values <- sapply(stats, function(n){print(n$p.value)})
p.values[p.values == TRUE] <- 1
q.values <- p.adjust(p.values, method = "bonferroni")

svg(file = paste("/Users/eling01/Google Drive/scMouse_Immun/Manuscript/Figures/Fig4C_B6_scatter.svg"),  width = 7, height = 7)
plot(apply(B6.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.young)}), apply(B6.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(B6.old)}), pch = 16, col = ifelse(q.values < 0.1, "red", "black"), xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red")
dev.off()

# Now for CAST
stats <- list()
for(i in 1:length(shared.a)){
  bin.t <- binom.test(x = length(which(CAST.young[shared.a[i],] > 0)), n = ncol(CAST.young), p = length(which(CAST.old[shared.a[i],] > 0))/ncol(CAST.old))
  stats[[i]] <- bin.t
}

p.values <- sapply(stats, function(n){print(n$p.value)})
p.values[p.values == TRUE] <- 1
p.values[p.values == FALSE] <- 0
q.values <- p.adjust(p.values, method = "bonferroni")

plot(apply(CAST.young[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.young)}), apply(CAST.old[shared.a,], 1, function(n){length(which(n > 0))/ncol(CAST.old)}), pch = 16, col = ifelse(q.values < 0.1, "red", "black"), xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red")




