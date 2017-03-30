################################################################################
#### Code to generate Figure S3 ################################################
################################################################################
library(pheatmap)
library(ggplot2)

exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

# Load BASiCS Test file comparing B6 young naive and CAST young naive
load("/path/to/Test_logFC2_file/")
B6_CAST <- Test_logFC2$Table

### Fig S3A
B6_CAST$cols <- ifelse(B6_CAST$ResultDiffExp == "SS19_naive_SS25_naive+", "CAST", ifelse(B6_CAST$ResultDiffExp == "SS51_naive_SS52_naive+", "Mus", "NoDiff"))

ggplot(B6_CAST, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols), size = 2) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("dark red", "dark blue", "grey")) + 
  theme_minimal(base_size = 20)

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

### Fig S3B
DE.CAST.map$cols <- ifelse(DE.CAST.map$ResultDiffExp == "CAST_SS19_naive_CAST_SS25_naive+", "CAST", ifelse(DE.CAST.map$ResultDiffExp == "Mus_SS19_naive_Mus_SS25_naive+", "Mus", "NoDiff"))

ggplot(DE.CAST.map, aes(log10(ExpRef + 1), log10(ExpTest + 1))) +
  geom_point(aes(col = cols), size = 2) + 
  ylim(0,6) + xlim(0,6) +
  scale_colour_manual(values = c("dark red", "dark blue", "grey")) + 
  theme_minimal(base_size = 20)

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


### Fig S3C
# Plot species specifc heatmaps
pheatmap(log10(Data[which(!is.na(match(rownames(exp.data), B6.spec))),] + 1), fontsize = 2, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(93), cellwidth = 2, cellheight = 2)
pheatmap(log10(Data[which(!is.na(match(rownames(exp.data), CAST.spec))),] + 1), fontsize = 2, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, gaps_col = c(93), cellwidth = 2, cellheight = 2)


### Fig S3D
# See Fig 3C

### Fig S3E
# Global variability and species-specific variability
boxplot(list(B6.all = log10(B6_CAST$OverDispRef[B6_CAST$ResultDiffExp != "ExcludedByUser"]),
             B6.spec = log10(B6_CAST$OverDispRef[match(rownames(genenames)[match(B6.spec, as.character(genenames[,2]))], B6_CAST$GeneNames)]),
             CAST.all = log10(B6_CAST$OverDispTest[B6_CAST$ResultDiffExp != "ExcludedByUser"]),
             CAST.spec = log10(B6_CAST$OverDispTest[match(rownames(genenames)[match(CAST.spec, as.character(genenames[,2]))], B6_CAST$GeneNames)])),
        pch = 16)

### Test for significance
wilcox.test(log10(B6_CAST$OverDispRef[B6_CAST$ResultDiffExp != "ExcludedByUser"]), log10(B6_CAST$OverDispRef[match(rownames(genenames)[match(B6.spec, as.character(genenames[,2]))], B6_CAST$GeneNames)]))
wilcox.test(log10(B6_CAST$OverDispTest[B6_CAST$ResultDiffExp != "ExcludedByUser"]), log10(B6_CAST$OverDispTest[match(rownames(genenames)[match(CAST.spec, as.character(genenames[,2]))], B6_CAST$GeneNames)]))
