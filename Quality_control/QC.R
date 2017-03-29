################################################################################
#### Code for computational filtering and quality control ######################
################################################################################
library(RColorBrewer)
library(Rtsne)
library(pheatmap)
library(ggplot2)
library(scran)
library(plot3Drgl)
setwd("/Users/eling01/GitHub/ImmuneAging2017/")

genenames <- read.table("Data/Genenames.txt", sep = "\t", header = TRUE)
rownames(genenames) <- genenames[,1]

### For quality control, also mitochondrial genes need to be read in
mit_genes <- read.table("/path/to/MT_genes.txt")

### Read in the HTSeq output 
input.CAST <- read.table("/path/to/output_CAST.txt", sep="\t")
input.B6 <- read.table("/path/to/output_B6.txt", sep="\t")

### Merge datsets based on gene orthologs
CAST.genes <- unlist(lapply(rownames(input.CAST)[which(grepl("ENSMU", rownames(input.CAST)))], function(x){unlist(strsplit(x, "\\."))[1]}))
B6.genes <- rownames(input.B6)
ERCC <- rownames(input.B6)[which(grepl("ERCC", rownames(input.B6)))]

common.genes <- intersect(B6.genes, CAST.genes)

### Adjust rownames for the CAST dataset if genenames are not consistent with B6
rownames(input.CAST)[which(grepl("ENSMU", rownames(input.CAST)))] <- unlist(lapply(CAST.genes, function(n){unlist(strsplit(n, "\\."))[1]}))

### Reassemble matrices
input.CAST <- rbind(input.CAST[common.genes,], input.CAST[ERCC,], tail(input.CAST, n=5))
input.B6 <- rbind(input.B6[common.genes,], input.B6[ERCC,], tail(input.B6, n=5))

### Merge datasets
input <- cbind(input.CAST, input.B6)

### Read in the decoding files for C1 experiments - one file per batch is needed
# Here I only demonstrate the code using two decoding files

Decoding.active <- read.table("/path/to/Decoding_active.txt", sep="\t") 
Decoding.naive <- read.table("/path/to/Decoding_naive.txt", sep="\t") 

### Order by position on plate to match with visual inspection
ordering <- order(Decoding.active[,2])

### Create data matrix to store quality features
Master.QC <- as.data.frame(matrix(data=NA, nrow=192, ncol=3))
colnames(Master.QC) <- c("PositionOnChip", "PositionOnPlate", "Microscopy")
rownames(Master.QC) <- c(paste("SS51_naive", ordering, sep="_"), paste("SS52_active", ordering, sep="_"))

Master.QC[1:96,1:3] <- as.matrix(Decoding.naive[ordering,])
Master.QC[(1*96+1):(2*96),] <- as.matrix(Decoding.active[ordering,])

### Expand the dataframe by features to examine
Master.QC$Genomic <- NA
Master.QC$ERCC <- NA
Master.QC$NoFeature <- NA
Master.QC$Ambigous <- NA
Master.QC$LowQual <- NA
Master.QC$NotAligned <- NA
Master.QC$NotUnique <- NA
Master.QC$TotalReads <- NA
Master.QC$Mitochondrial <- NA
Master.QC$GenesDetected <- NA
Master.QC$Libraries <- colnames(input)

### Fill in the dataframe with quality features - for demostration purposes only two batches are used
Master.QC[1:192,4:13] <- matrix(data = as.numeric(c(colSums(input[1:31089,]), colSums(input[31090:31181,]), input[31182,], input[31183,], input[31184,], input[31185,], input[31186,], colSums(input), colSums(input[as.character(mit_genes[,1]),]), unlist(apply(input[1:31089,], 2, function(n){length(which(n > 0))})))), ncol=10, nrow=192, byrow = FALSE)

### Remove libraries for which the mapping did not work - only 1 or 2 in the full dataset
Master.QC <- Master.QC[-which(is.na(Master.QC$Genomic)),]

### Collect the batch information
chips <- unlist(lapply(rownames(Master.QC), function(n){paste(unlist(strsplit(n, split = "_"))[1:2], collapse = "_")}))

### Fig S1B
plot((Master.QC$Genomic/Master.QC$TotalReads), pch=16, col= ifelse(Master.QC$Microscopy == "-", "red", ifelse(Master.QC$Microscopy == "single", "black", ifelse(Master.QC$Microscopy == "2 cells", "yellow", "green"))), xlab="%Genomic", ylab="%ERCC")

### Collect only single cells
Singles <- Master.QC[Master.QC$Microscopy == "single",]

### Fig S1C
plot(Singles$Genomic/Singles$TotalReads, Singles$ERCC/Singles$TotalReads, col= ifelse(Singles$ERCC/Singles$TotalReads < 0.5 & Singles$Genomic/Singles$TotalReads > 0.2, "black", "red"), xlab="%Genomic", ylab="%ERCC", pch=16)

### Collect only cells with high exonic content
Singles.ERCC <- Singles[which(Singles$ERCC/Singles$TotalReads < 0.5 & Singles$Genomic/Singles$TotalReads > 0.2),]

### Fig S1D
plot(log10(Singles.ERCC$TotalReads), pch=16, xlab="Cells", ylab="log[Total reads]", col= ifelse(log10(Singles.ERCC$TotalReads) > 6, "black", "red"), cex=0.7)
abline(h = mean(log10(Singles.ERCC$TotalReads)), col = "red")

### Collect cells with more than 1mio reads 
Singles.Total <- Singles.ERCC[which(log10(Singles.ERCC$TotalReads) > 6),]

### Fig S1E
plot(Singles.Total$GenesDetected, pch=16, col=ifelse(Singles.Total$GenesDetected > 3000 | Singles.Total$GenesDetected < 1250, "red", "black"), xlab="Cells", ylab="# Genes Detected")

### Collect cells with more than 1250 and less than 3000 genes detected
Singles.final <- Singles.Total[which(Singles.Total$GenesDetected < 3000 & Singles.Total$GenesDetected > 1250),]

### Quality control for mitochondirial genes can only be done based on reads mapped to the B6 genome
B6.mus <- read.table("/path/to/B6_B6mappped.txt", header = TRUE, sep = "\t")
CAST.mus <- read.table("/path/to/CAST_B6mappped.txt", header = TRUE, sep = "\t")

B6.mus <- B6.mus[,which(!is.na(match(colnames(B6.mus), Singles.final$Libraries)))]
CAST.mus <- CAST.mus[,which(!is.na(match(colnames(CAST.mus), Singles.final$Libraries)))]

input.mus <- cbind(B6.mus, CAST.mus)

Mitos <- data.frame(row.names = colnames(input.mus))
Mitos$mitos <- c(colSums(B6.mus[as.character(mit_genes[,1]),])/colSums(B6.mus[1:45513,]),
                 colSums(CAST.mus[as.character(mit_genes[,1]),])/colSums(CAST.mus[1:45513,]))

### Fig S1F
plot(Mitos$mitos, pch=16, col = ifelse(Mitos$mitos > 0.1 | Mitos$mitos < 0.005, "red", "black"))

### Collect cells with more than 0.5% and less than 10% mitochondrial reads
Singles.Mito <- Singles.final[which(!is.na(match(Singles.final$Libraries, rownames(Mitos)[Mitos$mitos < 0.1 & Mitos$mitos > 0.005]))),]

### Fig S1G
# Check for batch effects
input.1 <- input[-which(grepl("ERCC-", rownames(input))),]
input.1 <- input.1[,which(grepl("SS51", colnames(input)) | grepl("SS52", colnames(input)))]
pca <- prcomp(log10(t(input.1) + 1))

chips <- sapply(colnames(input.1), function(n){paste(unlist(strsplit(n, "_"))[1:2], collapse = "_")})

df <- data.frame(row.names = colnames(input.1), PC1=pca$x[,1], PC2=pca$x[,2], activation = ifelse(grepl("active", colnames(input.1)), "active", "naive"), batch = sapply(colnames(input.1), function(n){paste(unlist(strsplit(n, "_"))[1:2], collapse = "_")}))

ggplot(df, aes(x=PC1, y=PC2)) + geom_point(aes(shape=activation, colour=batch), size=6) + scale_shape_manual(values=c(124,95))

### Fig S1H
# Check for cell contaminations
# Immune genes - again use the reads mapped to B6 genome since not all immune genes are annotated for CAST
input.mus <- input.mus[1:45513,as.character(Singles.Mito$Libraries)]
input.mus <- input.mus[match(unique(as.character(genenames[as.character(rownames(input.mus)),2])), as.character(genenames[as.character(rownames(input.mus)),2])),]
rownames(input.mus) <- as.character(genenames[as.character(rownames(input.mus)),2])
genelist <- c("Cd8a", "Cd4", "Cd3e", "Cd19", "H2-Aa", "Cd28", "Il2ra", "Cd44", "Trac", "Sell", "Il7r", "Cd69", "Klrg1")
order.cells <- colnames(input.mus)[input.mus["H2-Aa",] > 20]
cd8 <- colnames(input.mus)[input.mus["Cd8a",] > 20]
order.cells <- append(order.cells, cd8[which(is.na(match(cd8, order.cells)))])
il2ra <- colnames(input.mus)[input.mus["Il2ra",] > 0]
order.cells <- append(order.cells, il2ra[which(is.na(match(il2ra, order.cells)))])
order.cells <- append(order.cells, colnames(input.mus)[which(is.na(match(colnames(input.mus), order.cells)))])

cells <- rownames(Singles.Mito)[match(order.cells, Singles.Mito$Libraries)]

# Plot the heatmap with annotating activation state, age and species
pheatmap(log10(input.mus[genelist,order.cells] + 1), fontsize_col = 6, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), border_color = NA, annotation_col = data.frame(row.names = order.cells, strain = ifelse(grepl("SS19", cells) | grepl("SS25", cells) | grepl("SS53", cells) | grepl("SS54", cells), "CAST", "B6"), age = ifelse(grepl("SS19", cells) | grepl("SS25", cells) | grepl("SS51", cells) | grepl("SS52", cells), "young", "old"), activation = ifelse(grepl("active", cells), "active", "naive")))

# Pre-normalisation to match library sizes
rpm <- (input.mus/colSums(input.mus))*1000000

# Remove CB8 T cells and B cells
CD4 <- colnames(input.mus)[which(rpm["Cd8a",] < 20 & rpm["H2-Aa",] < 20 & rpm["Cd19",] < 20)]

Singles.Immu <- Singles.Mito[which(!is.na(match(Singles.Mito$Libraries, CD4))),]

### Fig S1I
# Remove activated cells that were not fully activated (eg. cells from activated chips that cluster with naive cells) 
input <- input.all[,as.character(Singles.Immu$Libraries)]
tsne <- Rtsne(t(log10(input + 1)), perplexity = 10)

cells <- rownames(Singles.Immu)
df <- data.frame(row.names = colnames(input), tSNE1 = tsne$Y[,1], tSNE2 = tsne$Y[,2],
                 strain = ifelse(grepl("SS19", cells) | grepl("SS25", cells) | grepl("SS53", cells) | grepl("SS54", cells), "CAST", "B6"), 
                 age = ifelse(grepl("SS19", cells) | grepl("SS25", cells) | grepl("SS51", cells) | grepl("SS52", cells), "young", "old"), 
                 activation = ifelse(grepl("active", cells), "active", "naive"))

library(ggplot2)

ggplot(df, aes(x=tSNE1, y=tSNE2)) + geom_point(aes(colour = strain, shape = activation), size = 6) + scale_shape_manual(values=c(124, 95))

# manually select cells that need to be removed
x <- rownames(Singles.Immu)[which(tsne$Y[,2] > 20)]
x[which(grepl("active", x))]
x <- rownames(Singles.Immu)[which(tsne$Y[,1] > 18)]
x[which(grepl("active", x))]

### inactive
inactive <- c("SS19_active_5",  "SS19_active_80", "SS25_active_5",  "SS25_active_89", "SS25_active_20", "SS25_active_68",
              "SS51_active_78", "SS53_active_77", "SS54_active_77", "SS54_active_73", "SS54_active_71", "SS63_active_15",
              "SS64_active_40", "SS64_active_42")
Singles.active <- Singles.Immu[-match(inactive, rownames(Singles.Immu)),]

################################################################################
#### From this point on filtering on the cells is done #########################
#### Now the characterisation of cells starts ##################################
################################################################################

### Fig S2A
# Estimate cell cycle state for each cell
exp.data <- read.table("Data/normalized_data.txt")

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(as.matrix(exp.data), mm.pairs)
head(assigned$scores)

phase <- rep("S", ncol(exp.data))
names(phase) <- colnames(exp.data)
phase[assigned$scores$G1 > 0.5] <- "G1"
phase[assigned$scores$G2M > 0.5] <- "G2M"
phase[assigned$scores$G1 > 0.5 & assigned$scores$G2M > 0.5] <- "unknown"
table(phase)
phase[phase == "S"]
phase[phase == "G1"]
phase[phase == "G2M"]

# Plot the results
scatter3D(assigned$normalized.scores$G1, assigned$normalized.scores$S, assigned$normalized.scores$G2M,
          colvar = match(phase, unique(phase)), pch = 16, xlab = "G1", ylab = "S", zlab = "G2M")

### Characterise possible differentiation processes after activation
### Fig S2H
markers.diff <- c("Tbx21", "Ifng", "Gata3", "Il4", "Rorc", "Il17a", "Foxp3", "Il10", "Bcl6", "Il21")

# Check gene expression for raw counts before gene filtering
input.raw <- read.table("Data/raw_data.txt", header = TRUE, sep = "\t")
genenames <- read.table("Data/Genenames.txt", header = TRUE, sep = "\t")
rownames(genenames) <- as.character(genenames[,1])

dat.sc.young <- input.raw[,c(which(grepl("SS51_naive", colnames(input.raw))), which(grepl("SS52_naive", colnames(input.raw))), which(grepl("SS51_active", colnames(input.raw))), which(grepl("SS52_active", colnames(input.raw))))]
dat.sc.old <- input.raw[,c(which(grepl("SS63_naive", colnames(input.raw))), which(grepl("SS64_naive", colnames(input.raw))), which(grepl("SS63_active", colnames(input.raw))), which(grepl("SS64_active", colnames(input.raw))))]

genes <- rownames(genenames)[match(markers.diff, genenames[,2])]
pheatmap(log10(cbind(dat.sc.young[genes,],dat.sc.old[genes,]) + 1), labels_row = genenames[genes,2], cellwidth = 2, cellheight = 8, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), gaps_col = c(93, 146, 268), gaps_row = c(2,4,6,8,10))

### Fig S2I - exhaustion markers
markers.exhaustion <- c("Cd5", "Pdcd1", "Lag3", "Havcr2", "Ctla4")

dat.sc.naive <- input.raw[,c(which(grepl("SS51_naive", colnames(input.raw))), which(grepl("SS52_naive", colnames(input.raw))), which(grepl("SS63_naive", colnames(input.raw))), which(grepl("SS64_naive", colnames(input.raw))))]
dat.sc.active <- input.raw[,c(which(grepl("SS51_active", colnames(input.raw))), which(grepl("SS52_active", colnames(input.raw))), which(grepl("SS63_active", colnames(input.raw))), which(grepl("SS64_active", colnames(input.raw))))]

genes <- rownames(genenames)[match(markers.exhaustion, genenames[,2])]
pheatmap(log10(cbind(dat.sc.naive[genes,], dat.sc.active[genes,]) + 1), labels_row = genenames[genes,2], cellwidth = 2, cellheight = 8, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE, col = colorRampPalette(colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(100), gaps_col = c(93, 215, 269), gaps_row = c(1))

### Fig S2J + K
# Just replace the name of the gene to plot the boxplots
gene <- "Ifng"
gene <- rownames(genenames)[match(gene, genenames[,2])]
par(mar=c(10,4,4,4))
boxplot(list(B6_young_naive = as.numeric(log10(input.raw[gene,which(grepl("SS51_naive", colnames(input.raw)) | grepl("SS52_naive", colnames(input.raw)))] + 1)),
             CAST_young_naive = as.numeric(log10(input.raw[gene,which(grepl("SS19_naive", colnames(input.raw)) | grepl("SS25_naive", colnames(input.raw)))] + 1)),
             B6_young_active = as.numeric(log10(input.raw[gene,which(grepl("SS51_active", colnames(input.raw)) | grepl("SS52_active", colnames(input.raw)))] + 1)),
             CAST_young_active = as.numeric(log10(input.raw[gene,which(grepl("SS19_active", colnames(input.raw)) | grepl("SS25_active", colnames(input.raw)))] + 1)),
             B6_old_naive = as.numeric(log10(input.raw[gene,which(grepl("SS63_naive", colnames(input.raw)) | grepl("SS64_naive", colnames(input.raw)))] + 1)),
             CAST_old_naive = as.numeric(log10(input.raw[gene,which(grepl("SS53_naive", colnames(input.raw)) | grepl("SS54_naive", colnames(input.raw)))] + 1)),
             B6_old_active = as.numeric(log10(input.raw[gene,which(grepl("SS63_active", colnames(input.raw)) | grepl("SS64_active", colnames(input.raw)))] + 1)),
             CAST_old_active = as.numeric(log10(input.raw[gene,which(grepl("SS53_active", colnames(input.raw)) | grepl("SS54_active", colnames(input.raw)))] + 1))), pch=16, las=2, ylim=c(0,4))

