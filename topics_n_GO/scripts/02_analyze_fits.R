library(tidyverse)
library(fastTopics) 
library(Seurat)
library(scales)
library(funfuns)
#

# This script defines functions to analyse the topic models
# produces figures and GO term enrichments to help with exploratory data analysis


# function that will run the analysis on a topic model 
analyze_fit <- function(SINGLE_FIT, SEURAT, SEED, GROUP, KMEANS_RANGE=NULL){
  # browser()
  set.seed(SEED)
  results <- list()
  NUM_TOPICS <- ncol(SINGLE_FIT$L)
  
  # Set colors for each topic
  if (NUM_TOPICS < 10){
    topic_colors <- RColorBrewer::brewer.pal(NUM_TOPICS, 'Set1')
    names(topic_colors) <- paste0('k', 1:NUM_TOPICS)
    topic_colors2 <- topic_colors
    names(topic_colors2) <- sub('k','',names(topic_colors2))
  } else {
    topic_colors <- colorRampPalette(colors = RColorBrewer::brewer.pal(9, 'Set1'))(NUM_TOPICS)
    names(topic_colors) <- paste0('k', 1:NUM_TOPICS)
    topic_colors2 <- topic_colors
    names(topic_colors2) <- sub('k','',names(topic_colors2))
  }
  
  # extract topic membership weight for cells and genes
  cellTopics <- as.data.frame(SINGLE_FIT$L) # find weighted topic membership for cells when using 7 topics
  geneTopics <- as.data.frame(SINGLE_FIT$F)
  
  # add topic membership to cel metadata in Seurat object
  noK <- paste0('k', 1:NUM_TOPICS)
  SEURAT <- AddMetaData(SEURAT, metadata = c(cellTopics)) 
  
  feature_plot <- 
    FeaturePlot(SEURAT, features = c(noK), 
                reduction = 'tsne',
                ncol = 3) & 
    scale_color_gradientn( colours = c('grey90', 'darkgreen'),  limits = c(0, 1), oob = squish) & 
    NoAxes() & NoLegend()
  
  # feature_plot$plot_env <- rlang::new_environment()
  
  Idents(SEURAT)<- SEURAT$phyloorder
  violin_plot <- 
    VlnPlot(SEURAT, features = c(noK), 
            pt.size = 0, 
            stack = TRUE) + 
    scale_fill_manual(values = topic_colors)
  
  # violin_plot$plot_env <- rlang::new_environment()
  # Now let's look at structure plots with different groupings:
  struct_plot_1 <- 
    structure_plot(SINGLE_FIT,
                   topics = 1:NUM_TOPICS,
                   colors = topic_colors2,
                   grouping = Idents(SEURAT), 
                   gap = 5)
  
  # struct_plot_1$plot_env <- rlang::new_environment()
  
  struct_plot_2 <- 
    structure_plot(SINGLE_FIT,
                   topics = 1:NUM_TOPICS,
                   colors = topic_colors2,
                   #grouping = Idents(il), 
                   gap = 5)
  # struct_plot_2$plot_env <- rlang::new_environment()
  
  # write these out?
  cellTopics <- cellTopics %>% rownames_to_column(var = 'CellBarcodes')
  geneTopics <- geneTopics %>% rownames_to_column(var = 'Gene')
  
  # kmeans clusters K=3 to K=10
  #TODO#  CHANGE TO MAP SAFELY
  if (!is.null(KMEANS_RANGE)){
    
    KMEANS_RES <- lapply(KMEANS_RANGE, kmeans_subroutine, SEURAT=SEURAT, SINGLE_FIT=SINGLE_FIT, topic_colors=topic_colors, SEED=SEED)
    names(KMEANS_RES) <- paste('KM', KMEANS_RANGE, sep = '')
    
    # clean up Kmeans results #
    
    # cell kmeans cluster membership at each level of K
    KMEANS_topic_clusters <- 
      map(KMEANS_RES, pluck, 3) %>% 
      purrr::reduce(left_join)
    
    # extract just the kmeans plots
    # removes the KM_clusters entry
    NEW_KMEANS_RES <- purrr::map(KMEANS_RES,function(x, y) x[names(x) != 'KM_clusters'])
    
    results[['KMEANS_RES']] <- NEW_KMEANS_RES
    results[['KMEANS_topic_clusters']] <- KMEANS_topic_clusters
    
  }
  
  
  ## DE and Volc plots
  counts <- 
    GetAssayData(SEURAT,assay = 'RNA', slot = 'counts') %>% 
    Matrix::t() %>%   #transpose: now cells as rows, genes as columns, 
    floor()
  
  
  # dim(SINGLE_FIT$L)
  # dim(SINGLE_FIT$F)
  # dim(counts)
  
  counts <- counts[rownames(counts) %in% rownames(SINGLE_FIT$L),colnames(counts) %in% rownames(SINGLE_FIT$F)]
  
  # dim(SINGLE_FIT$L)
  # dim(SINGLE_FIT$F)
  # dim(counts)
  
  ## Genes enriched in topics in general
  dfa_out <- diff_count_analysis(fit = SINGLE_FIT, X = counts)
  
  ### Genes enriched in topic based clusters 
  # DO this for each
  # dfc_clusters <- diff_count_clusters(clusters,counts)
  # entry of this
  # KMEANS_topic_clusters
  # 
  
  
  # print(dfa_out$F0)
  # print(dfa_out$F1)
  all_topics_volcano_plots <- function(dfa_out){
    res <- list()
    for (K in 1:NUM_TOPICS){
      p <- volcano_plot(diff_count_result = dfa_out, k=K)
      res[[K]] <- p
    }
    return(res)
  }
  
  VOLCs <- all_topics_volcano_plots(dfa_out)
  
  
  PVALS <- dfa_out$pval %>%
    as.data.frame() %>%
    rownames_to_column(var='Gene') %>%
    gather(-Gene, key = 'topic', value='pval') #%>% 
    # mutate(newP=10^(-pval))
  
  LFCS <- dfa_out$beta %>%
    as.data.frame() %>% 
    rownames_to_column(var='Gene') %>%
    gather(-Gene, key = 'topic', value='LFC')
  
  ZSCR <- dfa_out$Z %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    gather(-Gene, key='topic', value='Zscore')
  
  # TODO Select and write this out?
  # Only positive LFC, means genes enriched in topics relative to other topics
  # omits genes down in topics relative to other topics
  topic_DE_filename <- paste0('./outputs/',GROUP, '_K', NUM_TOPICS, '_Topic_DE_genes.tsv')
  
  topics_DE_genes <-
    LFCS %>% 
    left_join(PVALS) %>% 
    left_join(ZSCR) %>% 
    mutate(newP=10^(-pval))
  
  topics_sig_genes <- 
    topics_DE_genes %>% 
  filter(newP < 0.05) %>% 
    write_tsv(topic_DE_filename)
    
  topics_enriched_genes <- 
    topics_DE_genes %>% 
    dplyr:: select(Gene, topic, Zscore) %>% 
    group_by(topic) %>%
    nest() %>% #ungroup() %>% 
    mutate(enrich_val = map_dbl(data, ~ quantile(.x$Zscore, .975)), 
           enrich_genes= map2(.x=data, .y=enrich_val, ~ filter(.x, Zscore > .y))) %>% 
    dplyr::select(topic, enrich_genes) %>% 
    unnest(cols=enrich_genes)# %>% 
    #write_tsv(topic_DE_filename)
  
  # hist(topics_DE_genes$LFC)
  
  # Topic_enriched_genes <- 
  #   dfa_out$Z %>%
  #   as.data.frame() %>% 
  #   rownames_to_column(var = 'Gene') %>% 
  #   gather(-Gene, key='topic', value='Zscore') %>% 
  #   filter(LFC > 0) %>% 
  #   group_by(topic) %>%
  #   nest() %>% 
  #   mutate(enrich_val = map_dbl(data, ~ quantile(.x$Zscore, .95)), 
  #          enrich_genes= map2(.x=data, .y=enrich_val, ~ filter(.x, LFC > .y)))
  # 
  # 
  # Topic_enriched_genes %>% unnest(cols = )
  
  # NEED TO GENERATE APPROPRIATE MAPPING FILE
  # all genes detected in Seurat object
  All_GO <- read_tsv('./outputs/All_gene_to_GO.tsv') %>% 
    filter(Name %in% colnames(counts)) %>% 
    write_tsv('TEMP_GO_BACKGROUND.tsv')
  
  GO_filename <- paste0('./outputs/',GROUP, '_K', NUM_TOPICS, '_GO_terms.tsv')
  
  GO_filename
  
  GO_terms_4_topics <- 
    topics_enriched_genes %>% 
    group_by(topic) %>% 
    nest(enrich_genes=c(Gene,Zscore)) %>% 
    # dplyr::select(topic, enrich_genes) %>% 
    mutate(topGO_res=map(.x = enrich_genes, ~ topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = 'TEMP_GO_BACKGROUND.tsv', ont = 'BP')), 
           filtered_res=map(.x = topGO_res, ~ filter(.x, pval < 0.05))) %>% 
    dplyr::select(topic, filtered_res) %>%
    unnest(cols = filtered_res) %>% 
    ungroup()# %>% 
    #write_tsv(GO_filename)
  # now writing GO results in dedicated GO script.
  
  ### overview plot ###
  
  # both up and down?
  # overview_genes <- 
  #   LFCS %>% 
  #   left_join(PVALS) %>% 
  #   filter(pval < 0.05) %>% 
  #   group_by(topic) %>%
  #   nest() %>% 
  #   mutate(enrich_val = map_dbl(data, ~ quantile(.x$LFC, .99)), 
  #          enrich_genes= map2(.x=data, .y=enrich_val, ~ filter(.x, abs(LFC) > .y)))
  # 
  
  # rownames(SINGLE_FIT$F)
  # rownames(SINGLE_FIT$L)
  ###
  # build final results object
  #results[['topic_overview_plot']] <- topic_overview_plot
  results[['feature_plot']] <- feature_plot
  results[['violin_plot']] <- violin_plot
  results[['struct_plot_1']] <- struct_plot_1
  results[['struct_plot_2']] <- struct_plot_2
  
  results[['topics_sig_genes']] <- topics_sig_genes
  # results[['cellTopics']] <- cellTopics
  # results[['geneTopics']] <- geneTopics
  results[['volcano_plots']] <- VOLCs
  results[['GO_terms']] <- GO_terms_4_topics
  return(results)
}

# called by analyse_fit()
# runs KMeans clustering at K=3 to K=10 
# produces plots and cluster assignments for each Kmeans K
kmeans_subroutine <- 
  function(K, SINGLE_FIT, SEURAT, topic_colors, SEED){
    topic_colors2 <- topic_colors
    names(topic_colors2) <- sub('k','',names(topic_colors2))
    
    print(paste0('KMEANS ', K))
    set.seed(SEED)
    results <- list()
    NUM_TOPICS <- ncol(SINGLE_FIT$L)
    
    pca <- prcomp(SINGLE_FIT$L)$x
    KMEANS <- kmeans(pca,centers = K,iter.max = 10000)
    KM_clusters <- data.frame(CellBarcodes=names(KMEANS$cluster), 
                              topicClusters= KMEANS$cluster)
    
    
    STRUCT_PLOT <- 
      structure_plot(SINGLE_FIT,topics = 1:NUM_TOPICS, colors = topic_colors2,
                     grouping = factor(KM_clusters$topicClusters), gap = 25)
    # STRUCT_PLOT$plot_env <- rlang::new_environment()
    
    SEURAT <- AddMetaData(SEURAT, metadata = c(KM_clusters))
    DIM_PLOT <- 
      DimPlot(SEURAT, group.by = 'topicClusters', reduction = 'tsne', label = TRUE)
    # DIM_PLOT$plot_env <- rlang::new_environment()
    
    # didn't do here, but we should also export the cluster membership with between 3 and 10 clusters used.... maybe collate into single file and include cell barcode IDs??
    colnames(KM_clusters)[2] <- paste('Kmeans', K, 'topicClusters', sep = '_')
    # KM_clusters
    results[['STRUCT_PLOT']] <- STRUCT_PLOT
    results[['DIM_PLOT']] <- DIM_PLOT
    results[['KM_clusters']] <- KM_clusters
    # results[['DE_clusters']] <- DE_clusters
    
    return(results)
    
  }


analyze_fits <- function(SEURAT, FITS, SEED, GROUP){
  results_name <- paste0('./results/', GROUP, '_results.rds')
  safe_analyse_fit <- safely(.f=analyze_fit)
  res <- map(FITS, safe_analyse_fit, SEURAT=SEURAT, SEED=SEED, GROUP=GROUP)
  # print('Saving results...')
  # write_rds(res, results_name)
  # print('Done saving results!')
  return(res) # 
}



# # This is having problems, everything NULL....
# # 
# SEURAT <- read_rds('./data/GutBlood_IntegratedILCs.rds')
# FITS <- read_rds('./topic_model_fits/GutBlood_IntegratedILCs_K_3_K_25_topic_model_fits.rds')
# 
# GROUP <- 'GutBlood_ILC' # probably can find a way to cut this
# 
# results <- analyze_fits(SEURAT = SEURAT, FITS = FITS[1], GROUP = 'GutBlood_ILC', SEED=2)
