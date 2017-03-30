#############################################################
#### Code to generate Figure 1 ##############################
#############################################################
library(Rtsne)
library(pheatmap)
library(ggplot2)
library(BASiCS)
setwd("/Users/eling01/GitHub/ImmuneAging2017/")

# Load gene names
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

### Figure 2
# Fig 2A
exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")

# Collect naive and activated cells from young B6
Data <- exp.data[,c(which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data))),
                       which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data))))]

# Calculate tSNE
tsne <- Rtsne(t(log10(Data + 1)), perplexity = 10)
tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], activation = c(rep("-", 93), rep("+", 53)))

# Plot tSNE
ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
  geom_point(size = 8, aes(colour = activation),shape = "-") +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_minimal() +
  ylab("tSNE 1") +
  xlab("tSNE 2") 


#### Fig 2B
# Scatterplot of DE genes
# Load a BASiCS Test dataframe tested on logFC=2 (Test_logFC2) and logFC=0 (Test_logFC0) between naive and activated cells
for.plot <- Test_logFC2$Table
for.plot.DV <- Test_logFC0$Table

for.plot$cols <- ifelse(for.plot$ResultDiffExp == "SS51_active_SS52_active+", "Active", ifelse(for.plot$ResultDiffExp == "SS51_naive_SS52_naive+", "Naive", ifelse(for.plot.DV$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("dark red", "black", "dark blue", "grey")) + 
  theme_minimal(base_size = 20)

# Fig 2C
# Plot the variability measures in naive and activated cells for genes taht don't change in mean expression

var.gene <- data.frame(exp = log10(c(for.plot.DV$OverDispRef[which(for.plot.DV$ResultDiffExp == "NoDiff")], for.plot.DV$OverDispTest[which(for.plot.DV$ResultDiffExp == "NoDiff")])),
                       cond = factor(c(rep("Naive", length(for.plot.DV$OverDispRef[which(for.plot.DV$ResultDiffExp == "NoDiff")])),
                                       rep("Active", length(for.plot.DV$OverDispTest[which(for.plot.DV$ResultDiffExp == "NoDiff")]))), levels = c("Naive", "Active")))

ggplot(var.gene, aes(cond, exp)) + 
  geom_boxplot(aes(fill = cond), outlier.shape = NA) +
  theme_minimal(base_size = 30) +
  ylim(0,2.5) +
  xlab("State") + ylab("log10[Over-dispersion]") + 
  scale_fill_manual(values = c("red", "blue"))

# Test for difference
wilcox.test(log10(for.plot.DV$OverDispRef[which(for.plot.DV$ResultDiffExp == "NoDiff")]),
            log10(for.plot.DV$ExpTest[which(for.plot.DV$ResultDiffExp == "NoDiff")]))

# Fig 2C
# Calculate the frequency of expression of up- and down-regulated genes in naive and active

active.cells <- exp.data[as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS51_active_SS52_active+"]),
                         which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]
naive.cells <- exp.data[as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS51_naive_SS52_naive+"]),
                        which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]


### Histogram of frequency - make sure the genes are downsampled otherwise the y-scale is off
plot(hist(apply(naive.cells[sample(1:nrow(naive.cells), 600),], 1, function(n){length(which(n > 0))/ncol(naive.cells)}), breaks = 50), xlim = c(0,1),ylim = c(0,70))
abline(v = median(apply(naive.cells, 1, function(n){length(which(n > 0))/ncol(naive.cells)})))

plot(hist(apply(active.cells[sample(1:nrow(active.cells), 600),], 1, function(n){length(which(n > 0))/ncol(active.cells)}), breaks = 50), xlim = c(0,1), ylim = c(0,70))
abline(v = median(apply(active.cells, 1, function(n){length(which(n > 0))/ncol(active.cells)})))


# Fig 2E-G
# Downsample to 50 cells for visualisation purposes
cells.naive <- colnames(exp.data)[which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]
cells.active <- colnames(exp.data)[which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]

gene <-"Il2ra"
gene.id <-rownames(genenames)[match(gene, as.character(genenames[,2]))]

barplot(as.numeric(log10(exp.data[gene.id,cells.naive[sample(1:length(cells.naive), 50)]] + 1)), ylim = c(0,6))
abline(h = log10(for.plot$ExpRef[match(gene.id, for.plot$GeneNames)]))

barplot(as.numeric(log10(exp.data[gene.id,cells.active[sample(1:length(cells.active), 50)]] + 1)), ylim = c(0,6))
abline(h = log10(for.plot$ExpTest[match(gene.id, for.plot$GeneNames)]))

