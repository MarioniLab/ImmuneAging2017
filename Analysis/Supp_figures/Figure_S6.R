################################################################################
#### Code to generate Figure S6 ################################################
################################################################################
library(pheatmap)
library(ggplot2)

exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

########## Dimensionlity reduction
### Fig S6A
Data <- exp.data[,which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)) |
                                 grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)) |
                                 grepl("SS63_naive", colnames(exp.data)) | grepl("SS64_naive", colnames(exp.data)) | 
                                 grepl("SS53_naive", colnames(exp.data)) | grepl("SS54_naive", colnames(exp.data)))]

species <- vector(length = ncol(Data))
species[which(grepl("SS51_naive", colnames(Data)) | grepl("SS52_naive", colnames(Data)))] <- "B6_young" 
species[which(grepl("SS19_naive", colnames(Data)) | grepl("SS25_naive", colnames(Data)))] <- "CAST_young" 
species[which(grepl("SS63_naive", colnames(Data)) | grepl("SS64_naive", colnames(Data)))] <- "B6_old" 
species[which(grepl("SS53_naive", colnames(Data)) | grepl("SS54_naive", colnames(Data)))] <- "CAST_old" 

pca <- prcomp(t(log10(Data + 1)), center = TRUE)

pca.df <- data.frame(pca.1 = pca$x[,1], pca.2 = pca$x[,2], species = species)

ggplot(data = pca.df, aes(pca.1, pca.2)) + 
  geom_point(data = subset(pca.df, species == "B6_young" | species == "CAST_young"), size = 7, col = "black") +
  geom_point(data = subset(pca.df, species == "B6_young" | species == "CAST_young"), aes(colour = species), size = 6) +
  geom_point(data = subset(pca.df, species == "B6_old" | species == "CAST_old"), size = 9, col = "black", pch = 17) +
  geom_point(data = subset(pca.df, species == "B6_old" | species == "CAST_old"), aes(colour = species), size = 6, pch = 17) +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") 

### Fig S6B
Data <- exp.data[,which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)) |
                          grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)) |
                          grepl("SS63_active", colnames(exp.data)) | grepl("SS64_active", colnames(exp.data)) | 
                          grepl("SS53_active", colnames(exp.data)) | grepl("SS54_active", colnames(exp.data)))]

species <- vector(length = ncol(Data))
species[which(grepl("SS51_active", colnames(Data)) | grepl("SS52_active", colnames(Data)))] <- "B6_young" 
species[which(grepl("SS19_active", colnames(Data)) | grepl("SS25_active", colnames(Data)))] <- "CAST_young" 
species[which(grepl("SS63_active", colnames(Data)) | grepl("SS64_active", colnames(Data)))] <- "B6_old" 
species[which(grepl("SS53_active", colnames(Data)) | grepl("SS54_active", colnames(Data)))] <- "CAST_old" 

pca <- prcomp(t(log10(Data + 1)), center = TRUE)

pca.df <- data.frame(pca.1 = pca$x[,1], pca.2 = pca$x[,2], species = species)

ggplot(data = pca.df, aes(pca.1, pca.2)) + 
  geom_point(data = subset(pca.df, species == "B6_young" | species == "CAST_young"), size = 7, col = "black") +
  geom_point(data = subset(pca.df, species == "B6_young" | species == "CAST_young"), aes(colour = species), size = 6) +
  geom_point(data = subset(pca.df, species == "B6_old" | species == "CAST_old"), size = 9, col = "black", pch = 17) +
  geom_point(data = subset(pca.df, species == "B6_old" | species == "CAST_old"), aes(colour = species), size = 6, pch = 17) +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") 

############ Differential expression
# For each comparison, a BASiCS Test (Test_logFC2) object is needed. Optionally, testing on logFC=0 can be incorporated (Test_logFC0)
### Fig S6C
# B6 naive
for.plot <- Test_logFC2$Table
for.plot.DV <- Test_logFC0$Table

for.plot$cols <- ifelse(for.plot$ResultDiffExp == "SS63_naive_SS64_naive+", "Old", ifelse(for.plot$ResultDiffExp == "SS51_naive_SS52_naive+", "Young", ifelse(for.plot.DV$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("black", "dark blue", "grey", "dark red")) + 
  theme_minimal(base_size = 20)

# Fraction of cells expressing DE genes
young.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS51_naive_SS52_naive+"])
young.cells <- exp.data[,which(grepl("SS51_naive", colnames(exp.data)) | grepl("SS52_naive", colnames(exp.data)))]
old.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS63_naive_SS64_naive+"])
old.cells <- exp.data[,which(grepl("SS63_naive", colnames(exp.data)) | grepl("SS64_naive", colnames(exp.data)))]
plot(density(apply(young.cells[young.genes[sample(1:length(young.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(young.cells)})), ylim = c(0,9), xlim = c(0,1))
plot(density(apply(old.cells[old.genes[sample(1:length(old.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(old.cells)})), ylim = c(0,9), xlim = c(0,1))


# CAST naive
for.plot <- Test_logFC2$Table
for.plot.DV <- Test_logFC0$Table

for.plot$cols <- ifelse(for.plot$ResultDiffExp == "SS53_naive_SS54_naive+", "Old", ifelse(for.plot$ResultDiffExp == "SS19_naive_SS25_naive+", "Young", ifelse(for.plot.DV$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("black", "dark blue", "grey", "dark red")) + 
  theme_minimal(base_size = 20)

# Fraction of cells expressing DE genes
young.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS19_naive_SS25_naive+"])
young.cells <- exp.data[,which(grepl("SS19_naive", colnames(exp.data)) | grepl("SS25_naive", colnames(exp.data)))]
old.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS53_naive_SS54_naive+"])
old.cells <- exp.data[,which(grepl("SS53_naive", colnames(exp.data)) | grepl("SS54_naive", colnames(exp.data)))]
plot(density(apply(young.cells[young.genes[sample(1:length(young.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(young.cells)})), ylim = c(0,9), xlim = c(0,1))
plot(density(apply(old.cells[old.genes[sample(1:length(old.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(old.cells)})), ylim = c(0,9), xlim = c(0,1))

### Fig S6D
# B6 active
for.plot <- Test_logFC2$Table
for.plot.DV <- Test_logFC0$Table

for.plot$cols <- ifelse(for.plot$ResultDiffExp == "SS63_active_SS64_active+", "Old", ifelse(for.plot$ResultDiffExp == "SS51_active_SS52_active+", "Young", ifelse(for.plot.DV$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("black", "dark blue", "grey", "dark red")) + 
  theme_minimal(base_size = 20)

# Fraction of cells expressing DE genes
young.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS51_active_SS52_active+"])
young.cells <- exp.data[,which(grepl("SS51_active", colnames(exp.data)) | grepl("SS52_active", colnames(exp.data)))]
old.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS63_active_SS64_active+"])
old.cells <- exp.data[,which(grepl("SS63_active", colnames(exp.data)) | grepl("SS64_active", colnames(exp.data)))]
plot(density(apply(young.cells[young.genes[sample(1:length(young.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(young.cells)})), ylim = c(0,9), xlim = c(0,1))
plot(density(apply(old.cells[old.genes[sample(1:length(old.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(old.cells)})), ylim = c(0,9), xlim = c(0,1))


# CAST active
for.plot <- Test_logFC2$Table
for.plot.DV <- Test_logFC0$Table

for.plot$cols <- ifelse(for.plot$ResultDiffExp == "SS53_active_SS54_active+", "Old", ifelse(for.plot$ResultDiffExp == "SS19_active_SS25_active+", "Young", ifelse(for.plot.DV$ResultDiffExp == "NoDiff", "logFC0", "NoDiff")))

ggplot(for.plot, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols)) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("black", "dark blue", "grey", "dark red")) + 
  theme_minimal(base_size = 20)

# Fraction of cells expressing DE genes
young.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS19_active_SS25_active+"])
young.cells <- exp.data[,which(grepl("SS19_active", colnames(exp.data)) | grepl("SS25_active", colnames(exp.data)))]
old.genes <- as.character(for.plot$GeneNames[for.plot$ResultDiffExp == "SS53_active_SS54_active+"])
old.cells <- exp.data[,which(grepl("SS53_active", colnames(exp.data)) | grepl("SS54_active", colnames(exp.data)))]
plot(density(apply(young.cells[young.genes[sample(1:length(young.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(young.cells)})), ylim = c(0,9), xlim = c(0,1))
plot(density(apply(old.cells[old.genes[sample(1:length(old.genes), 100)],], 1, function(n){length(which(n > 0))/ncol(old.cells)})), ylim = c(0,9), xlim = c(0,1))


# Overlap in aging related gene expression differences
# Jaccar index
jac.index <- function(a,b){
  length(which(!is.na(match(a, b))))/(length(a) + length(b) - length(which(!is.na(match(a, b)))))
}

### Fig S6E - DE comparison between young and old in naive cells
DE.B6 <- Test_logFC2_B6_naive$Table
DE.B6.young <- DE.B6[DE.B6$ResultDiffExp == "SS51_naive_SS52_naive+",]
DE.B6.old <- DE.B6[DE.B6$ResultDiffExp == "SS63_naive_SS64_naive+",]
DE.CAST <- Test_logFC2_CAST_naive$Table
DE.CAST.young <- DE.CAST[DE.CAST$ResultDiffExp == "SS19_naive_SS25_naive+",]
DE.CAST.old <- DE.CAST[DE.CAST$ResultDiffExp == "SS53_naive_SS54_naive+",]

jac.index(DE.B6.young$GeneNames, DE.CAST.young$GeneNames)
jac.index(DE.B6.old$GeneNames, DE.CAST.old$GeneNames)

### Fig S6F - DE comparison between young and old in active cells
DE.B6 <- Test_logFC2_B6_active$Table
DE.B6.young <- DE.B6[DE.B6$ResultDiffExp == "SS51_naive_SS52_naive+",]
DE.B6.old <- DE.B6[DE.B6$ResultDiffExp == "SS63_naive_SS64_naive+",]
DE.CAST <- Test_logFC2_CAST_active$Table
DE.CAST.young <- DE.CAST[DE.CAST$ResultDiffExp == "SS19_naive_SS25_naive+",]
DE.CAST.old <- DE.CAST[DE.CAST$ResultDiffExp == "SS53_naive_SS54_naive+",]

jac.index(DE.B6.young$GeneNames, DE.CAST.young$GeneNames)
jac.index(DE.B6.old$GeneNames, DE.CAST.old$GeneNames)
