#####################################################################
#### DE and DV testing based on previously simulated MCMC Chains  ###
#####################################################################

# Set current working directory 
setwd("/Users/eling01/GitHub/ImmuneAging2017/")
library(BASiCS)

### For testing 2 MCMC objects are needed - both should contain the same genes
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

# Read in ERCC concentrations
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

# Calculate total number of ERCC molecules per well
ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/50000 # dilution factor
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

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
genes <- rownames(CountsRef)[which(grepl("ENS", rownames(CountsRef)))]
sel <- rep(TRUE, length(genes))
sel[match(excl.genes, genes)] <- FALSE

# Testing for Differential expression 
Test_logFC2 <- BASiCS_D_TestDE(Data_DV, MCMC_DV, GeneNames = genes,
                               EpsilonM = 2, EpsilonD = 0.4, GenesSelect = sel, 
                               EFDR_M = 0.05, EFDR_D = 0.05,
                               OrderVariable = "GeneNames", GroupLabelRef = "Naive", 
                               GroupLabelTest = "Active")

# Testing for Differential variability
Test_logFC0 <- BASiCS_D_TestDE(Data_DV, MCMC_DV, GeneNames = genes,
                               EpsilonM = 0, EpsilonD = 0.4, GenesSelect = sel, 
                               EFDR_M = 0.05, EFDR_D = 0.05,
                               OrderVariable = "GeneNames",  GroupLabelRef = "Naive", 
                               GroupLabelTest = "Active")
