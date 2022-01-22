library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(Seurat)
library(tidyverse)

set.seed(1)

# Defines a function to fit multiple topic models from a Seurat object
# Saves the result in an rds file

fit_topic_models <- 
  function(seurat_rds_path, kstart, kend){
    set.seed(1)
    dat_name <- sub('./data/(.*).rds', '\\1', seurat_rds_path)
    # browser()
    dat <- readRDS(seurat_rds_path)
    
    
    counts <- 
      GetAssayData(dat,assay = 'RNA', slot = 'counts') %>% 
      t() %>%   #transpose: now cells as rows, genes as columns, 
      floor() # converts to counts (round down) 
    
    counts <- counts[,colSums(counts) != 0]
    
    
    fits <- list()
    for (kval in kstart:kend){
      
      kname <- paste0('k', kval)
      print(paste('Fitting topic model with K=',kval))
      
      # fit <- 'TEST'
      fit <- fit_topic_model(counts, k = kval)
      fits[[kname]] <- fit

    }
    # print(paste(dat_name, 'K', kstart,'K', kend, 'topic_model_fits.rds', sep = '_'))
    write_rds(fits, file = paste(dat_name, 'K', kstart,'K', kend, 'topic_model_fits.rds', sep = '_'))
    return(fits)

  
}

# Fit topic models from K=3 to K=10 on each data subset

capture.output(fit_topic_models(seurat_rds_path = './data/GutBlood_IntegratedILCs.rds', kstart = 3, kend = 10),
               file = 'GutBlood_ILCs_K3_K25.log')


capture.output(fit_topic_models(seurat_rds_path = './data/BplasmaSeurat.rds', kstart = 3, kend = 10),
               file = 'Bplasma_K10_K25.log')


capture.output(fit_topic_models(seurat_rds_path = './data/MonoDCSeurat.rds', kstart = 3, kend = 10),
               file = 'MonoDC_K10_K25.log')


capture.output(fit_topic_models(seurat_rds_path = './data/IleumAllProcessed.rds', kstart = 3, kend = 10),
               file = 'IleumAll_K10_K25.log')

capture.output(fit_topic_models(seurat_rds_path = './data/IleumAtlasBonly.rds', kstart = 3, kend = 10),
               file = 'IleumAtlasBonly_K10_K25.log')

capture.output(fit_topic_models(seurat_rds_path = './data/IleumAtlasTILConly.rds', kstart = 3, kend = 10),
               file = 'IleumAtlasTILConly_K10_K25.log')

