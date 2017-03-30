################################################################################
#### Code to generate Figure S7 ################################################
################################################################################
library(pheatmap)
library(ggplot2)
library(xlsx)

exp.data <- read.table("Data/normalized_data.txt", header = TRUE, sep = "\t")
genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

### Fig S7A
Data <- exp.data[,c(which(grepl("SS51_active", colnames(exp.data))), which(grepl("SS52_active", colnames(exp.data))),
            which(grepl("SS63_active", colnames(exp.data))), which(grepl("SS64_active", colnames(exp.data))),
            which(grepl("SS19_active", colnames(exp.data))), which(grepl("SS25_active", colnames(exp.data))),
            which(grepl("SS53_active", colnames(exp.data))), which(grepl("SS54_active", colnames(exp.data))))]

# Load shared activation genes
shared <- read.xlsx("Data/S3.xlsx", sheetIndex = 1)
shared <- as.character(shared$Gene.ID)

# order by mean expression
shared <- shared[order(rowMeans(Data[shared,]), decreasing = TRUE)]

pheatmap(log10(Data[shared,] + 1), fontsize = 4, cluster_cols = FALSE, cluster_rows = FALSE, gaps_col = c(53, 127, 162), col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100))

### Fig S7B
# Remove cells from the activated state that don't show correct activation markers
# Only on B6 cells

input <- read.table("Data/raw_data.txt", header = TRUE, sep = "\t")
active <- input[,which(grepl("SS51_active", colnames(input)) | grepl("SS52_active", colnames(input)) | grepl("SS63_active", colnames(input)) | grepl("SS64_active", colnames(input)))]

# Cd69 = ENSMUSG00000030156
plot(as.numeric(log10(active["ENSMUSG00000030156",order(log10(active["ENSMUSG00000030156",] + 1), decreasing = FALSE)] + 1)), pch = 16)
active <- active[,-which(log10(active["ENSMUSG00000030156",] + 1) < 2.5)]

# Sell/Cd62l = ENSMUSG00000026581
plot(as.numeric(log10(active["ENSMUSG00000026581",order(log10(active["ENSMUSG00000026581",] + 1), decreasing = FALSE)] + 1)), pch = 16)
active <- active[,-which(log10(active["ENSMUSG00000026581",] + 1) > 1)]

# Trac = ENSMUSG00000076928
plot(as.numeric(log10(active["ENSMUSG00000076928",order(log10(active["ENSMUSG00000076928",] + 1), decreasing = FALSE)] + 1)), pch = 16)

# Il2ra = ENSMUSG00000026770
plot(as.numeric(log10(active["ENSMUSG00000026770",order(log10(active["ENSMUSG00000026770",] + 1), decreasing = FALSE)] + 1)), pch = 16)
active <- active[,-which(log10(active["ENSMUSG00000026770",] + 1) < 2)]

# Ifng = ENSMUSG00000055170
plot(as.numeric(log10(active["ENSMUSG00000055170",order(log10(active["ENSMUSG00000055170",] + 1), decreasing = FALSE)] + 1)), pch = 16)
active <- active[,-which(log10(active["ENSMUSG00000055170",] + 1) > 0)]

# Heatmpas before filtering
# old active
pheatmap(log10(input[c("ENSMUSG00000026581", "ENSMUSG00000055170", "ENSMUSG00000002897", "ENSMUSG00000005087", "ENSMUSG00000026012", "ENSMUSG00000032093", "ENSMUSG00000076928", "ENSMUSG00000026770", "ENSMUSG00000030156"),grepl("SS63_active", colnames(input)) | grepl("SS64_active", colnames(input))] + 1), cluster_rows = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100))
# young active
pheatmap(log10(input[c("ENSMUSG00000026581", "ENSMUSG00000055170", "ENSMUSG00000002897", "ENSMUSG00000005087", "ENSMUSG00000026012", "ENSMUSG00000032093", "ENSMUSG00000076928", "ENSMUSG00000026770", "ENSMUSG00000030156"),grepl("SS51_active", colnames(input)) | grepl("SS52_active", colnames(input))] + 1), cluster_rows = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100))

# After filtering
# old active
pheatmap(log10(active[c("ENSMUSG00000026581", "ENSMUSG00000055170", "ENSMUSG00000002897", "ENSMUSG00000005087", "ENSMUSG00000026012", "ENSMUSG00000032093", "ENSMUSG00000076928", "ENSMUSG00000026770", "ENSMUSG00000030156"),grepl("SS63_active", colnames(active)) | grepl("SS64_active", colnames(active))] + 1), cluster_rows = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100))
# young active
pheatmap(log10(active[c("ENSMUSG00000026581", "ENSMUSG00000055170", "ENSMUSG00000002897", "ENSMUSG00000005087", "ENSMUSG00000026012", "ENSMUSG00000032093", "ENSMUSG00000076928", "ENSMUSG00000026770", "ENSMUSG00000030156"),grepl("SS51_active", colnames(active)) | grepl("SS52_active", colnames(active))] + 1), cluster_rows = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100))

### Fig S7C
# Run the BASiCS model on activated cells after filtering and compare between filtered young activated B6 and filtered old activated B6
# Analysis can be seen in Fig 4A

### Fig S7D
# Run the BASiCS model on activated cells from SS57+SS58+SS59 and compare the variability measures to activated cells from young B6 (SS51+SS52)
# Analysis can be seen in Fig 4A

### Fig S7E
# Downsample cells to a total of 30 per condition, run the BASiCS model on the downsampled populations
# Perform analysis as in Fig 4A