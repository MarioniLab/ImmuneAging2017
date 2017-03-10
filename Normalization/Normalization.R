#######################################################
#### Normalization of scRNAseq data with BASiCS #######
#######################################################

# Set current working directory 
setwd("/Users/eling01/GitHub/ImmuneAging2017/")
library(BASiCS)

# Read in data
input <- read.table("Data/raw_data.txt", sep = "\t", header = TRUE)

# Split data into biological genes and ERCC spike ins
bio.g <- input[which(grepl("ENSMUS", rownames(input))),]
ERCC <- input[which(grepl("ERCC", rownames(input))),]

# Read in ERCC concentrations
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

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
plot(MCMC_Output, Param = "delta", Gene = 1)
plot(MCMC_Output, Param = "mu", Gene = 1)
plot(MCMC_Output, Param = "s", Cell = 1)

#### Retrieving normalized counts ###############################
DenoisedCounts = BASiCS_DenoisedCounts(Data = Data, Chain = MCMC_Output)
