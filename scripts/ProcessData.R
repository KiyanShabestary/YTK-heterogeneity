# Prepare Jackson dataset for visualization

require(scater) 
require(scran) 
require(SingleCellExperiment)
require(SummarizedExperiment)
require(BiocParallel)
require(dplyr)
require(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory to current location


# Custom normalization (Jackson2020) ----------------------------------------------------

# Load dataset as sce type
sce.conditions <- readRDS(file = file.path("/Users/kshab/Documents/ICL/MacroHetFinder/Rdata","Jackson2020_sce.Rds"))

# Could also add Jackson script from tsv

# Normalize dataset
sce.conditions <- lapply(sce.conditions, logNormCounts)

# Save dataset for shiny
#saveRDS(sce.conditions, file = file.path("/Users/kshab/Documents/shiny/ShinyScRNASeq-v2/data","Jackson2020_lognorm.Rds"))


# Select subset (log counts) -----------------------------------------------------------
# Only select YPD, MinimalGlucose, YPEtOH, CStarve conditions
# Only select WT genotype

# Filter condition
#sce.conditions.subset <- sce.conditions[c('YPD','MinimalGlucose','MinimalEthanol','CStarve')]

# Extract countings
YPD.countings = sce.conditions[["YPD"]]@assays@data@listData[["logcounts"]]
MinimalGlucose.countings = sce.conditions[["MinimalGlucose"]]@assays@data@listData[["logcounts"]]
YPEtOH.countings = sce.conditions[["YPEtOH"]]@assays@data@listData[["logcounts"]]
CStarve.countings = sce.conditions[["CStarve"]]@assays@data@listData[["logcounts"]]

# Extract labels
YPD.labels = sce.conditions[["YPD"]]@colData@listData[["Genotype_Group"]]
MinimalGlucose.labels = sce.conditions[["MinimalGlucose"]]@colData@listData[["Genotype_Group"]]
YPEtOH.labels = sce.conditions[["YPEtOH"]]@colData@listData[["Genotype_Group"]]
CStarve.labels = sce.conditions[["CStarve"]]@colData@listData[["Genotype_Group"]]

# Get index for WT phenotypes
YPD.idx = which(YPD.labels %in% "WT(ho)")
MinimalGlucose.idx = which(MinimalGlucose.labels %in% "WT(ho)")
YPEtOH.idx = which(YPEtOH.labels %in% "WT(ho)")
CStarve.idx = which(CStarve.labels %in% "WT(ho)")

# Keep countings associated to WT phenotype
YPD.countings.WT = YPD.countings[,YPD.idx]
MinimalGlucose.countings.WT = MinimalGlucose.countings[,MinimalGlucose.idx]
YPEtOH.countings.WT = YPEtOH.countings[,YPEtOH.idx]
CStarve.countings.WT = CStarve.countings[,CStarve.idx]

# Invert and add common gene names
geneNames <- read.table(file = file.path("/Users/kshab/Documents/shiny/ShinyScRNASeq-v2/data","yeast_gene_names.tsv"), sep = '\t', header = TRUE)

YPD.countings.WT.geneNames <- YPD.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

MinimalGlucose.countings.WT.geneNames <- MinimalGlucose.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

YPEtOH.countings.WT.geneNames <- YPEtOH.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

CStarve.countings.WT.geneNames <- CStarve.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

# Get rid of NA and blank
YPD.countings.WT.geneNames$Name <- ifelse(is.na(YPD.countings.WT.geneNames$Name), YPD.countings.WT.geneNames$SystematicName, YPD.countings.WT.geneNames$Name)
MinimalGlucose.countings.WT.geneNames$Name <- ifelse(is.na(MinimalGlucose.countings.WT.geneNames$Name), MinimalGlucose.countings.WT.geneNames$SystematicName, MinimalGlucose.countings.WT.geneNames$Name)
YPEtOH.countings.WT.geneNames$Name <- ifelse(is.na(YPEtOH.countings.WT.geneNames$Name), YPEtOH.countings.WT.geneNames$SystematicName, YPEtOH.countings.WT.geneNames$Name)
CStarve.countings.WT.geneNames$Name <- ifelse(is.na(CStarve.countings.WT.geneNames$Name), CStarve.countings.WT.geneNames$SystematicName, CStarve.countings.WT.geneNames$Name)

YPD.countings.WT.geneNames$Name <- ifelse(YPD.countings.WT.geneNames$Name == '', YPD.countings.WT.geneNames$SystematicName, YPD.countings.WT.geneNames$Name)
MinimalGlucose.countings.WT.geneNames$Name <- ifelse(MinimalGlucose.countings.WT.geneNames$Name == '', MinimalGlucose.countings.WT.geneNames$SystematicName, MinimalGlucose.countings.WT.geneNames$Name)
YPEtOH.countings.WT.geneNames$Name <- ifelse(YPEtOH.countings.WT.geneNames$Name == '', YPEtOH.countings.WT.geneNames$SystematicName, YPEtOH.countings.WT.geneNames$Name)
CStarve.countings.WT.geneNames$Name <- ifelse(CStarve.countings.WT.geneNames$Name == '', CStarve.countings.WT.geneNames$SystematicName, CStarve.countings.WT.geneNames$Name)

# Assemble in a list
countings.WT = list(YPD.countings.WT.geneNames, MinimalGlucose.countings.WT.geneNames, YPEtOH.countings.WT.geneNames, CStarve.countings.WT.geneNames)
names(countings.WT) = c("YPD","MinimalGlucose","YPEtOH","CStarve")

#exported to shiny
#saveRDS(countings.WT, file = file.path("/Users/kshab/Documents/shiny/ShinyScRNASeq-v2/data","Jackson2020_lognorm_WT.Rds"))

# Select subset (log counts) -----------------------------------------------------------

# Extract counts
YPD.countings = sce.conditions[["YPD"]]@assays@data@listData[["counts"]]
MinimalGlucose.countings = sce.conditions[["MinimalGlucose"]]@assays@data@listData[["counts"]]
YPEtOH.countings = sce.conditions[["YPEtOH"]]@assays@data@listData[["counts"]]
CStarve.countings = sce.conditions[["CStarve"]]@assays@data@listData[["counts"]]

# Extract labels
YPD.labels = sce.conditions[["YPD"]]@colData@listData[["Genotype_Group"]]
MinimalGlucose.labels = sce.conditions[["MinimalGlucose"]]@colData@listData[["Genotype_Group"]]
YPEtOH.labels = sce.conditions[["YPEtOH"]]@colData@listData[["Genotype_Group"]]
CStarve.labels = sce.conditions[["CStarve"]]@colData@listData[["Genotype_Group"]]

# Get index for WT phenotypes
YPD.idx = which(YPD.labels %in% "WT(ho)")
MinimalGlucose.idx = which(MinimalGlucose.labels %in% "WT(ho)")
YPEtOH.idx = which(YPEtOH.labels %in% "WT(ho)")
CStarve.idx = which(CStarve.labels %in% "WT(ho)")

# Keep countings associated to WT phenotype
YPD.countings.WT = YPD.countings[,YPD.idx]
MinimalGlucose.countings.WT = MinimalGlucose.countings[,MinimalGlucose.idx]
YPEtOH.countings.WT = YPEtOH.countings[,YPEtOH.idx]
CStarve.countings.WT = CStarve.countings[,CStarve.idx]

# Invert and add common gene names
geneNames <- read.table(file = file.path("/Users/kshab/Documents/shiny/ShinyScRNASeq-v2/data","yeast_gene_names.tsv"), sep = '\t', header = TRUE)

YPD.countings.WT.geneNames <- YPD.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

MinimalGlucose.countings.WT.geneNames <- MinimalGlucose.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

YPEtOH.countings.WT.geneNames <- YPEtOH.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

CStarve.countings.WT.geneNames <- CStarve.countings.WT %>% as.data.frame() %>%
  rownames_to_column('SystematicName') %>% merge(geneNames, all.x = TRUE, all.y=TRUE) %>%
  relocate(Name,.after = SystematicName)

# Get rid of NA and blank
YPD.countings.WT.geneNames$Name <- ifelse(is.na(YPD.countings.WT.geneNames$Name), YPD.countings.WT.geneNames$SystematicName, YPD.countings.WT.geneNames$Name)
MinimalGlucose.countings.WT.geneNames$Name <- ifelse(is.na(MinimalGlucose.countings.WT.geneNames$Name), MinimalGlucose.countings.WT.geneNames$SystematicName, MinimalGlucose.countings.WT.geneNames$Name)
YPEtOH.countings.WT.geneNames$Name <- ifelse(is.na(YPEtOH.countings.WT.geneNames$Name), YPEtOH.countings.WT.geneNames$SystematicName, YPEtOH.countings.WT.geneNames$Name)
CStarve.countings.WT.geneNames$Name <- ifelse(is.na(CStarve.countings.WT.geneNames$Name), CStarve.countings.WT.geneNames$SystematicName, CStarve.countings.WT.geneNames$Name)

YPD.countings.WT.geneNames$Name <- ifelse(YPD.countings.WT.geneNames$Name == '', YPD.countings.WT.geneNames$SystematicName, YPD.countings.WT.geneNames$Name)
MinimalGlucose.countings.WT.geneNames$Name <- ifelse(MinimalGlucose.countings.WT.geneNames$Name == '', MinimalGlucose.countings.WT.geneNames$SystematicName, MinimalGlucose.countings.WT.geneNames$Name)
YPEtOH.countings.WT.geneNames$Name <- ifelse(YPEtOH.countings.WT.geneNames$Name == '', YPEtOH.countings.WT.geneNames$SystematicName, YPEtOH.countings.WT.geneNames$Name)
CStarve.countings.WT.geneNames$Name <- ifelse(CStarve.countings.WT.geneNames$Name == '', CStarve.countings.WT.geneNames$SystematicName, CStarve.countings.WT.geneNames$Name)

# Assemble in a list
countings.WT = list(YPD.countings.WT.geneNames, MinimalGlucose.countings.WT.geneNames, YPEtOH.countings.WT.geneNames, CStarve.countings.WT.geneNames)
names(countings.WT) = c("YPD","MinimalGlucose","YPEtOH","CStarve")

#exported to shiny
saveRDS(countings.WT, file = file.path("/Users/kshab/Documents/shiny/ShinyScRNASeq-v2/data","Jackson2020_WT.Rds"))



