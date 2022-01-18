library(tidyverse)
library(fastTopics) 
library(Seurat)
library(scales)
library(funfuns)


# functions defined in this script
source('./scripts/analyze_fits.R')


### Jayne changed her mind on these....
# CD4 only: topics k = 3, spatial DGE modules = 3
# gd/CD8 only: topics = 3, spatial DGE modules = 4 and modules = 6
# ILC only: topics = 3, spatial DGE modules = 3
# B only: topics = 6, spatial DGE modules = 9
# Myeloid only: topics = 5, spatial = don't know yet; skip for now
# Non-immune only: topics = 6, spatial DGE modules = 9


results <- list()

#######
# CD4 only: topics k = 3,
CD4only_SEURAT <- read_rds('./data/Ileum_CD4Tonly.rds')
CD4only_FITS <- read_rds('./topic_model_fits/Ileum_CD4Tonly_K_3_K_10_topic_model_fits.rds')
# CD4only_GROUP <- 'CD4'
CD4Tonly_results <- analyze_fit(SEURAT = CD4only_SEURAT, SINGLE_FIT = CD4only_FITS$k3, GROUP='CD4', SEED=2)

results[['CD4']] <- CD4Tonly_results
######
# gd/CD8 only: topics = 3
gdCD8only_SEURAT <- read_rds('./data/Ileum_gdCD8Tonly.rds')
gdCD8only_FITS <- read_rds('./topic_model_fits/Ileum_gdCD8Tonly_K_3_K_10_topic_model_fits.rds')
# gdCD8only_GROUP <- 'Ileum_gdCD8only'
gdCD8only_results <- analyze_fit(SEURAT = gdCD8only_SEURAT, SINGLE_FIT = gdCD8only_FITS$k3, GROUP='gdCD8', SEED=2)

results[['gdCD8']] <- gdCD8only_results

####
# ILC only: topics = 3
ILConly_SEURAT <- read_rds('./data/Ileum_ILConly.rds')
ILConly_FITS <- read_rds('./topic_model_fits/Ileum_ILConly_K_3_K_10_topic_model_fits.rds')
# ILConly_GROUP <- 'Ileum_ILConly'
ILConly_results <- analyze_fit(SEURAT = ILConly_SEURAT, SINGLE_FIT = ILConly_FITS$k3, GROUP='ILC', SEED=2)

results[['ILC']] <- ILConly_results

### Bcell
Bonly_SEURAT <- read_rds('./data/Ileum_Bonly.rds')
Bonly_FITS <- read_rds('./topic_model_fits/Ileum_Bonly_K_3_K_10_topic_model_fits.rds')
# Bonly_GROUP <- 'Ileum_Bonly'

# topic K=3
Bonly_results <- analyze_fit(SEURAT = Bonly_SEURAT, SINGLE_FIT = Bonly_FITS$k3, GROUP='Bcells', SEED=2)

results[['BCell']] <- Bonly_results
####
# Myeloid only: topics = 5

Myeloidonly_SEURAT <- read_rds('./data/Ileum_MyeloidOnly.rds')
Myeloidonly_FITS <- read_rds('./topic_model_fits/Ileum_MyeloidOnly_K_3_K_10_topic_model_fits.rds')
# Myeloidonly_GROUP <- 'Ileum_Myeloidonly'
Myeloidonly_results <- analyze_fit(SEURAT = Myeloidonly_SEURAT, SINGLE_FIT = Myeloidonly_FITS$k5, GROUP='Myeloid', SEED=2)

results[['Myeloid']] <- Myeloidonly_results

####
# Non-immune only: topics = 6

NonImmumeonly_SEURAT <- read_rds('./data/Ileum_NonImmuneOnly.rds')
NonImmumeonly_FITS <- read_rds('./topic_model_fits/Ileum_NonImmuneOnly_K_3_K_10_topic_model_fits.rds')
# NonImmumeonly_GROUP <- 'Ileum_NonImmuneonly'
NonImmumeonly_results <- analyze_fit(SEURAT = NonImmumeonly_SEURAT, SINGLE_FIT = NonImmumeonly_FITS$k6, GROUP='NonImmume', SEED=2)

results[['NonImmune']] <- NonImmumeonly_results




