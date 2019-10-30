setwd("Dropbox (Personal)//BASiCSplus/Comparisons/Science/")
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
# Naive young B6 - copied to md
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Young" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "B6 young 1" |
                      CD4.metadata$Individuals == "B6 young 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Young" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "B6 young 1" |
                      CD4.metadata$Individuals == "B6 young 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data.B6.naive <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                             Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                      rep(TRUE, nrow(cur_ERCC))), 
                             SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.naive, "Data/Data_B6.naive.rds")

# Active young B6 - copied to md
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 young 1" |
                          CD4.metadata$Individuals == "B6 young 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Young" &
                   CD4.metadata$Stimulus == "Active" &
                   (CD4.metadata$Individuals == "B6 young 1" |
                      CD4.metadata$Individuals == "B6 young 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                             Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                      rep(TRUE, nrow(cur_ERCC))), 
                             SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active, "Data/Data_B6.active.rds")

# Active old B6 - copied to md
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus domesticus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "B6 old 1" |
                          CD4.metadata$Individuals == "B6 old 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus domesticus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Active" &
                   (CD4.metadata$Individuals == "B6 old 1" |
                      CD4.metadata$Individuals == "B6 old 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.B6.active.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                              Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                       rep(TRUE, nrow(cur_ERCC))), 
                              SpikeInfo = cur_SpikeInput)
saveRDS(Data.B6.active.old, "Data/Data_B6.active.old.rds")

# Naive young CAST - copied to md
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus castaneus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Unstimulated" &
                       (CD4.metadata$Individuals == "CAST young 1" |
                          CD4.metadata$Individuals == "CAST young 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus castaneus" &
                   CD4.metadata$Age == "Young" &
                   CD4.metadata$Stimulus == "Unstimulated" &
                   (CD4.metadata$Individuals == "CAST young 1" |
                      CD4.metadata$Individuals == "CAST young 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.CAST.naive <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                         rep(TRUE, nrow(cur_ERCC))), 
                                SpikeInfo = cur_SpikeInput)
saveRDS(Data.CAST.naive, "Data/Data_CAST.naive.rds")

# Active young CAST - copied to md
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus castaneus" &
                       CD4.metadata$Age == "Young" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "CAST young 1" |
                          CD4.metadata$Individuals == "CAST young 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus castaneus" &
                   CD4.metadata$Age == "Young" &
                   CD4.metadata$Stimulus == "Active" &
                   (CD4.metadata$Individuals == "CAST young 1" |
                      CD4.metadata$Individuals == "CAST young 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.CAST.active <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                 Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                          rep(TRUE, nrow(cur_ERCC))), 
                                 SpikeInfo = cur_SpikeInput)
saveRDS(Data.CAST.active, "Data/Data_CAST.active.rds")

# Active old CAST
cur_bio.g <- bio.g.1[,CD4.metadata$Strain == "Mus musculus castaneus" &
                       CD4.metadata$Age == "Old" &
                       CD4.metadata$Stimulus == "Active" &
                       (CD4.metadata$Individuals == "CAST old 1" |
                          CD4.metadata$Individuals == "CAST old 2") &
                       CD4.metadata$Celltype == "MACS-purified Naive"]

cur_ERCC <- ERCC[,CD4.metadata$Strain == "Mus musculus castaneus" &
                   CD4.metadata$Age == "Old" &
                   CD4.metadata$Stimulus == "Active" &
                   (CD4.metadata$Individuals == "CAST old 1" |
                      CD4.metadata$Individuals == "CAST old 2") &
                   CD4.metadata$Celltype == "MACS-purified Naive"] 
cur_ERCC <- cur_ERCC[rowSums(cur_ERCC) > 0,]

SpikeInput <- ERCC.count[rownames(cur_ERCC),2]
cur_SpikeInput <- data.frame("Name" = rownames(cur_ERCC),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)

Data.CAST.active.old <- newBASiCS_Data(Counts = rbind(cur_bio.g, cur_ERCC), 
                                     Tech = c(rep(FALSE, nrow(cur_bio.g)), 
                                              rep(TRUE, nrow(cur_ERCC))), 
                                     SpikeInfo = cur_SpikeInput)
saveRDS(Data.CAST.active.old, "Data/Data_CAST.active.old.rds")

chains.path <- "~/Dropbox/BASiCSplus/Comparisons/Science/MCMCs"

# Run MCMC
MCMC.B6.naive <- BASiCS_MCMC(Data = Data.B6.naive, 
                             N = 40000, 
                             Thin = 20, 
                             Burn = 20000, 
                             Regression = FALSE, 
                             WithSpikes = TRUE, 
                             PriorDelta = "log-normal",
                             StoreChains = TRUE,
                             StoreDir = chains.path,
                             RunName = "B6.naive")

MCMC.B6.active <- BASiCS_MCMC(Data = Data.B6.active, 
                              N = 40000, 
                              Thin = 20, 
                              Burn = 20000, 
                              Regression = FALSE, 
                              WithSpikes = TRUE, 
                              PriorDelta = "log-normal",
                              StoreChains = TRUE,
                              StoreDir = chains.path,
                              RunName = "B6.active")

MCMC.B6.active.old <- BASiCS_MCMC(Data = Data.B6.active.old, 
                                  N = 40000, 
                                  Thin = 20, 
                                  Burn = 20000, 
                                  Regression = FALSE, 
                                  WithSpikes = TRUE, 
                                  PriorDelta = "log-normal",
                                  StoreChains = TRUE,
                                  StoreDir = chains.path,
                                  RunName = "B6.active.old")

MCMC.CAST.naive <- BASiCS_MCMC(Data = Data.CAST.naive, 
                               N = 40000, 
                               Thin = 20, 
                               Burn = 20000, 
                               Regression = FALSE, 
                               WithSpikes = TRUE, 
                               PriorDelta = "log-normal",
                               StoreChains = TRUE,
                               StoreDir = chains.path,
                               RunName = "CAST.naive")

MCMC.CAST.active <- BASiCS_MCMC(Data = Data.CAST.active, 
                                N = 40000, 
                                Thin = 20, 
                                Burn = 20000, 
                                Regression = FALSE, 
                                WithSpikes = TRUE, 
                                PriorDelta = "log-normal",
                                StoreChains = TRUE,
                                StoreDir = chains.path,
                                RunName = "CAST.active")

MCMC.CAST.active.old <- BASiCS_MCMC(Data = Data.CAST.active.old, 
                                    N = 40000, 
                                    Thin = 20, 
                                    Burn = 20000, 
                                    Regression = FALSE, 
                                    WithSpikes = TRUE, 
                                    PriorDelta = "log-normal",
                                    StoreChains = TRUE,
                                    StoreDir = chains.path,
                                    RunName = "CAST.active.old")
