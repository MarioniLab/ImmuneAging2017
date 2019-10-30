setwd("~/Dropbox (Personal)/BASiCSplus/Comparisons/Science/Supplements/")
website <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4888/"
file <- "E-MTAB-4888.processed.1.zip"
destfile <- "raw_data.txt.zip"
download.file(paste(website, file, sep = ""), 
              destfile = destfile)
  
# Unzip file
unzip("raw_data.txt.zip")
file.remove("raw_data.txt.zip")

# Read in raw data
CD4.raw <- read.table("raw_data.txt", header = TRUE, sep = "\t")
bio.g <- CD4.raw[grepl("ENSM", rownames(CD4.raw)),]
ERCC <- CD4.raw[grepl("ERCC", rownames(CD4.raw)),]

# Remove lowly expressed genes
# Pre-scale reads to adjust for differences in library size
rpm <- (bio.g/colSums(bio.g))*1000000

# Filtering on genes
cell.count <- apply(rpm, 1, function(n){length(which(n > 20))})
bio.g.1 <- bio.g[names(which(cell.count > 2)),]

# Download raw counts file
file <- "E-MTAB-4888.additional.1.zip"
destfile <- "metadata.txt.zip"
download.file(paste(website, file, sep = ""), 
              destfile = destfile)
  
# Unzip file
unzip("metadata.txt.zip")
file.remove("metadata.txt.zip")

# Read in metadata file
CD4.metadata <- read.table("metadata_file.txt", header = TRUE, sep = "\t")

# Save library identifier as rownames
rownames(CD4.metadata) <- CD4.metadata$X

# Order metadata based on counts
CD4.metadata <- CD4.metadata[colnames(bio.g.1),]

# Read in the spike-in concentrations
website <- "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/"
file <- "cms_095046.txt"
ERCC.conc <- read.table(url(paste(website, file, sep = "")), 
                        sep = "\t", header = TRUE)

ERCC.mmul <- ERCC.conc$concentration.in.Mix.1..attomoles.ul. * (10^(-18))

ERCC.countmul <- ERCC.mmul*(6.02214076*(10^23))

ERCC.count <- ERCC.countmul / 50000

ERCC.count.final <- ERCC.count * 0.009

ERCC.count <- data.frame(row.names = ERCC.conc$ERCC.ID,
                         Names = ERCC.conc$ERCC.ID,
                         count = ERCC.count.final)

# Generate the SingleCellExperiment object
library(BASiCS)

# Naive old B6 
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "B6 old 1" |
                      CD4.metadata$Individuals == "B6 old 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "B6 old 1" |
                      CD4.metadata$Individuals == "B6 old 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data.B6.naive.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                             Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                      rep(TRUE, nrow(cur_ERCC))), 
                             SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.naive.old, "Data/Data_B6.naive.old.rds")

# Naive old CAST 
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus castaneus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Unstimulated" &
                       (CD4.metadata$Individuals == "CAST old 1" |
                          CD4.metadata$Individuals == "CAST old 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus castaneus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "CAST old 1" |
                      CD4.metadata$Individuals == "CAST old 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.CAST.naive.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                    Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                             rep(TRUE, nrow(cur_ERCC))), 
                                    SpikeInfo = cur_SpikeInput)
saveRDS(Data.CAST.naive.old, "Data/Data_CAST.naive.old.rds")


# Remove contaminating cell types
# B6 young active
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 young 1" |
                          CD4.metadata$Individuals == "B6 young 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000030156",] + 1) < 2.5)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000026581",] + 1) > 1)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000026770",] + 1) < 2)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000055170",] + 1) > 0)]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.woContamination <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                      Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                               rep(TRUE, nrow(cur_ERCC))), 
                                      SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.woContamination, "Data/Data_B6.active.woContamination.rds")

# B6 old active
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 old 1" |
                          CD4.metadata$Individuals == "B6 old 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000030156",] + 1) < 2.5)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000026581",] + 1) > 1)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000026770",] + 1) < 2)]
cur_bio.g <- cur_bio.g[,-which(log10(cur_bio.g["ENSMUSG00000055170",] + 1) > 0)]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.old.woContamination <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                                 Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                          rep(TRUE, nrow(cur_ERCC))), 
                                                 SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.old.woContamination, "Data/Data_B6.active,old.woContamination.rds")

# Additional replicates
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Activated" &
                       (CD4.metadata$Individuals == "B6 old 8" |
                          CD4.metadata$Individuals == "B6 old 9" |
                          CD4.metadata$Individuals == "B6 old 10") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Activated" &
                   (CD4.metadata$Individuals == "B6 old 8" |
                      CD4.metadata$Individuals == "B6 old 9" |
                      CD4.metadata$Individuals == "B6 old 10") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.old.replicates <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                      Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                               rep(TRUE, nrow(cur_ERCC))), 
                                      SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.old.replicates, "Data/Data_B6.active.old.replicates.rds")

# Downsampling
# B6 active young
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 young 1" |
                          CD4.metadata$Individuals == "B6 young 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

set.seed(12345)
cur_bio.g <- cur_bio.g[,sample(1:ncol(cur_bio.g), 30)]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.young.downsampled <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                                Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                         rep(TRUE, nrow(cur_ERCC))), 
                                                SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.young.downsampled, "Data/Data_B6.active.young.downsampled.rds")

# B6 active old
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 old 1" |
                          CD4.metadata$Individuals == "B6 old 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

set.seed(12345)
cur_bio.g <- cur_bio.g[,sample(1:ncol(cur_bio.g), 30)]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.old.downsampled <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                                   Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                            rep(TRUE, nrow(cur_ERCC))), 
                                                   SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.old.downsampled, "Data/Data_B6.active.old.downsampled.rds")

# FACS sorted Naive - active young
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Activated" &
                       (CD4.metadata$Individuals == "B6 young 3" |
                          CD4.metadata$Individuals == "B6 young 4") &
                       CD4.metadata$Celltype == "FACS-purified Naive"]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.FACSnaive.active <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                                 Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                          rep(TRUE, nrow(cur_ERCC))), 
                                                 SpikeInfo = cur_SpikeInput)
saveRDS(Data.FACSnaive.active, "Data/Data_FACSnaive.active.rds")

# FACS sorted Naive - active old
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Activated" &
                       (CD4.metadata$Individuals == "B6 old 3" |
                          CD4.metadata$Individuals == "B6 old 4") &
                       CD4.metadata$Celltype == "FACS-purified Naive"]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.FACSnaive.active.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                        Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                 rep(TRUE, nrow(cur_ERCC))), 
                                        SpikeInfo = cur_SpikeInput)
saveRDS(Data.FACSnaive.active.old, "Data/Data_FACSnaive.active.old.rds")

# FACS sorted Effector Memory - active young
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Activated" &
                       (CD4.metadata$Individuals == "B6 young 5" |
                          CD4.metadata$Individuals == "B6 young 6" |
                          CD4.metadata$Individuals == "B6 young 7") &
                       CD4.metadata$Celltype == "FACS-purified Effector Memory"]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.FACSem.active <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                        Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                 rep(TRUE, nrow(cur_ERCC))), 
                                        SpikeInfo = cur_SpikeInput)
saveRDS(Data.FACSem.active, "Data/Data_FACSem.active.rds")

# FACS sorted Naive - active old
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Activated" &
                       (CD4.metadata$Individuals == "B6 old 5" |
                          CD4.metadata$Individuals == "B6 old 6" |
                          CD4.metadata$Individuals == "B6 old 7") &
                       CD4.metadata$Celltype == "FACS-purified Effector Memory"]

cur_ERCC <- ERCC[,colnames(cur_bio.g)] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.FACSem.active.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                            Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                                     rep(TRUE, nrow(cur_ERCC))), 
                                            SpikeInfo = cur_SpikeInput)
saveRDS(Data.FACSem.active.old, "Data/Data_FACSem.active.old.rds")


# Run MCMC
library(BASiCS)
setwd("~/Dropbox (Personal)/BASiCSplus/Comparisons/Science/Supplements/")
chains.path <- "~/Dropbox (Personal)/BASiCSplus/Comparisons/Science/Supplements/MCMCs/"
cur_dat <- readRDS("Data/Data_B6.active.old.downsampled.rds")
MCMC.B6.active.old.downsampled <- BASiCS_MCMC(Data = cur_dat, 
                          N = 40000, 
                          Thin = 20, 
                          Burn = 20000, 
                          Regression = FALSE, 
                          WithSpikes = TRUE, 
                          PriorDelta = "log-normal",
                          StoreChains = TRUE,
                          StoreDir = chains.path,
                          RunName = "B6.active.old.downsampled")

cur_dat <- readRDS("Data/Data_B6.active.old.replicates.rds")
MCMC.B6.active.old.replicates <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "B6.active.old.replicates")

cur_dat <- readRDS("Data/Data_B6.active.old.woContamination.rds")
MCMC.B6.active.old.woContamination <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "B6.active.old.woContamination")

cur_dat <- readRDS("Data/Data_B6.active.woContamination.rds")
MCMC.B6.active.woContamination <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "B6.active.woContamination")

cur_dat <- readRDS("Data/Data_B6.active.young.downsampled.rds")
MCMC.B6.active.young.downsampled <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "B6.active.young.downsampled")

cur_dat <- readRDS("Data/Data_B6.naive.old.rds")
MCMC.B6.naive.old <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "B6.naive.old")

cur_dat <- readRDS("Data/Data_CAST.naive.old.rds")
MCMC.CAST.naive.old <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "CAST.naive.old")

cur_dat <- readRDS("Data/Data_FACSem.active.old.rds")
MCMC.FACSem.active.old <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "FACSem.active.old")

cur_dat <- readRDS("Data/Data_FACSem.active.rds")
MCMC.FACSem.active <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "FACSem.active")

cur_dat <- readRDS("Data/Data_FACSnaive.active.old.rds")
MCMC.FACSnaive.active.old <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "FACSnaive.active.old")

cur_dat <- readRDS("Data/Data_FACSnaive.active.rds")
MCMC.FACSnaive.active <- BASiCS_MCMC(Data = cur_dat, 
                                              N = 40000, 
                                              Thin = 20, 
                                              Burn = 20000, 
                                              Regression = FALSE, 
                                              WithSpikes = TRUE, 
                                              PriorDelta = "log-normal",
                                              StoreChains = TRUE,
                                              StoreDir = chains.path,
                                              RunName = "FACSnaive.active")

