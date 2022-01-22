library(readxl)
library(tidyverse)
library(Seurat)
library(funfuns)


### TODO

# copy needed inputs to data folder
# DE from topics?

#/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE

# system('mkdir ./outputs/trad_DE_GO')


#### TRADITIONAL DEG #####

all_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/AllCells_OverallDE.xlsx') %>%
  filter(avg_logFC > 0)


Overall_results <-
  all_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/All_gene_to_GO.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Overall_results %>% write_tsv('./outputs/trad_DE_GO/all_types_GO.tsv')

all_results <- read_tsv('tmp_all_overall_results.tsv')
all_results[3354,]
##

B_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/B_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)


Bresults <-
  B_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/Bcell_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Bresults %>% write_tsv('./outputs/trad_DE_GO/B_cell_GO.tsv')

CD4T_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/CD4Tonly_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)

CD4Tresults <-
  CD4T_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/CD4_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

CD4Tresults %>% write_tsv('./outputs/trad_DE_GO/CD4_GO.tsv')


###

gdCD8_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/gdCD8TOnly_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)

gdCD8results <-
  gdCD8_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/gdCD8T_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

gdCD8results %>% write_tsv('./outputs/trad_DE_GO/gdCD8_GO.tsv')

###

ILC_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/ILConly_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)

ILCresults <-
  ILC_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/ILC_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

ILCresults %>% write_tsv('./outputs/trad_DE_GO/ILC_GO.tsv')


#


Myeloid_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/MyeloidOnly_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)

Myeloidresults <-
  Myeloid_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/Myeloid_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Myeloidresults %>% write_tsv('./outputs/trad_DE_GO/Myeloid_GO.tsv')

NonImmune_overall <- readxl::read_xlsx('/home/Julian.Trachsel/scRNAseqIleumAtlas/Ileum/DE/NonImmuneOnly_OverallDE.xlsx')%>%
  filter(avg_logFC > 0)

NonImmuneresults <-
  NonImmune_overall %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate(enriched_genes=purrr::map(.x = data,.f = unlist)) %>% 
  mutate(GO_results=map(.x = enriched_genes,
                        .f= ~topGO_wrapper(myInterestingGenes = .x,
                                           mapping_file ='./outputs/NonImmune_GO_universe.tsv', ont = 'BP'))) %>% 
  mutate(filt_results=
           purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

NonImmuneresults %>% write_tsv('./outputs/trad_DE_GO/NonImmune_GO.tsv')


###  HAYSTACK GENE MODULES

# Spatial gene module lists for B (k=4), CD4 only (k=3), gd/CD8 only (k=4), and ILC only (k=3).

system('mkdir ./outputs/haystack_GO/')

Bcell_modules <- read_xlsx('./data/B_k4_GeneModules.xlsx')


B_cell_results <- 
  Bcell_modules %>% 
  group_by(res.hc.clusters) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, unlist)) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/Bcell_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(res.hc.clusters, filt_results) %>% 
  unnest(cols=filt_results)

B_cell_results %>% write_tsv('./outputs/haystack_GO/Bcell_k4_GO.tsv')



######

CD4T_modules <- read_xlsx('./data/CD4Tonly_k3_GeneModules.xlsx')


CD4T_results <- 
  CD4T_modules %>% 
  group_by(res.hc.clusters) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, unlist)) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/CD4_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(res.hc.clusters, filt_results) %>% 
  unnest(cols=filt_results)


CD4T_results %>% write_tsv('./outputs/haystack_GO/CD4T_k3_GO.tsv')

###

gdCD8_modules <- read_xlsx('./data/GDCD8Tonly_k4_GeneModules.xlsx')


gdCD8_results <- 
  gdCD8_modules %>% 
  group_by(res.hc.clusters) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, unlist)) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/gdCD8T_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(res.hc.clusters, filt_results) %>% 
  unnest(cols=filt_results)

gdCD8_results %>% write_tsv('./outputs/haystack_GO/gdCD8T_k4_GO.tsv')

###


ILC_modules <- read_xlsx('./data/ILConly_k3_GeneModules.xlsx')


ILC_results <- 
  ILC_modules %>% 
  group_by(res.hc.clusters) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, unlist)) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/ILC_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(res.hc.clusters, filt_results) %>% 
  unnest(cols=filt_results)

ILC_results %>% write_tsv('./outputs/haystack_GO/ILC_k3_GO.tsv')


######### TOPICS ###

# Topic DE gene lists for B (k=3), CD4 only (k=3), gd/CD8 only (k=3), and ILC only (k=3)
# system('mkdir ./outputs/topics_GO')


BCell_topic_results <- 
  read_tsv('./outputs/Bcells_K3_Topic_DE_genes.tsv') %>% #pull(newP) %>% hist()
  # filter(LFC >.5) %>%
  # filter(newP < 0.05) %>% 
  group_by(topic) %>% 
  nest() %>% 
  mutate(cutoff=map_dbl(data, ~quantile(.x$Zscore, .95))) %>% 
  mutate(strict_filt=map2(data, cutoff,  ~filter(.x, Zscore > .y)), 
         lax_filt   =map(data, ~filter(.x, LFC > .5))) %>% 
  mutate(strict_GO=
           map(strict_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/Bcell_GO_universe.tsv')), 
         lax_GO=
           map(lax_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/Bcell_GO_universe.tsv'))) %>% 
  mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
         lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(ends_with('sigs'))



BCell_topic_results %>% dplyr::select(topic, strict_GO_sigs) %>% unnest(cols = strict_GO_sigs) %>% write_tsv('./outputs/topics_GO/Bcell_strict_GO.tsv')
BCell_topic_results %>% dplyr::select(topic, lax_GO_sigs) %>% unnest(cols = lax_GO_sigs) %>%  write_tsv('./outputs/topics_GO/Bcell_lax_GO.tsv')
# 
# tst %>% 
#   mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
#          lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
#   dplyr::select(ends_with('sigs'))


###### CD4

# read_tsv('./outputs/CD4_K3_Topic_DE_genes.tsv')

CD4T_topic_results <- 
  read_tsv('./outputs/CD4_K3_Topic_DE_genes.tsv') %>% #pull(newP) %>% hist()
  # filter(LFC >.5) %>%
  # filter(newP < 0.05) %>% 
  group_by(topic) %>% 
  nest() %>% 
  mutate(cutoff=map_dbl(data, ~quantile(.x$Zscore, .95))) %>% 
  mutate(strict_filt=map2(data, cutoff,  ~filter(.x, Zscore > .y)), 
         lax_filt   =map(data, ~filter(.x, LFC > .5))) %>% 
  mutate(strict_GO=
           map(strict_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/CD4_GO_universe.tsv')), 
         lax_GO=
           map(lax_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/CD4_GO_universe.tsv'))) %>% 
  mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
         lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(ends_with('sigs'))


CD4T_topic_results %>% dplyr::select(topic, strict_GO_sigs) %>% unnest(cols = strict_GO_sigs) %>% write_tsv('./outputs/topics_GO/CD4_strict_GO.tsv')
CD4T_topic_results %>% dplyr::select(topic, lax_GO_sigs) %>% unnest(cols = lax_GO_sigs) %>%  write_tsv('./outputs/topics_GO/CD4_lax_GO.tsv')


###### gdCD8

# read_tsv('./outputs/CD4_K3_Topic_DE_genes.tsv')


gdCD8_topic_results <- 
  read_tsv('./outputs/gdCD8_K3_Topic_DE_genes.tsv') %>% #pull(newP) %>% hist()
  # filter(LFC >.5) %>%
  # filter(newP < 0.05) %>% 
  group_by(topic) %>% 
  nest() %>% 
  mutate(cutoff=map_dbl(data, ~quantile(.x$Zscore, .95))) %>% 
  mutate(strict_filt=map2(data, cutoff,  ~filter(.x, Zscore > .y)), 
         lax_filt   =map(data, ~filter(.x, LFC > .5))) %>% 
  mutate(strict_GO=
           map(strict_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/gdCD8T_GO_universe.tsv')), 
         lax_GO=
           map(lax_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/gdCD8T_GO_universe.tsv'))) %>% 
  mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
         lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(ends_with('sigs'))


gdCD8_topic_results %>% dplyr::select(topic, strict_GO_sigs) %>% unnest(cols = strict_GO_sigs) %>% write_tsv('./outputs/topics_GO/gdCD8_strict_GO.tsv')
gdCD8_topic_results %>% dplyr::select(topic, lax_GO_sigs) %>% unnest(cols = lax_GO_sigs) %>%  write_tsv('./outputs/topics_GO/gdCD8_lax_GO.tsv')



##### ILC

ILC_topic_results <- 
  read_tsv('./outputs/ILC_K3_Topic_DE_genes.tsv') %>% #pull(newP) %>% hist()
  # filter(LFC >.5) %>%
  # filter(newP < 0.05) %>% 
  group_by(topic) %>% 
  nest() %>% 
  mutate(cutoff=map_dbl(data, ~quantile(.x$Zscore, .95))) %>% 
  mutate(strict_filt=map2(data, cutoff,  ~filter(.x, Zscore > .y)), 
         lax_filt   =map(data, ~filter(.x, LFC > .5))) %>% 
  mutate(strict_GO=
           map(strict_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/ILC_GO_universe.tsv')), 
         lax_GO=
           map(lax_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/ILC_GO_universe.tsv'))) %>% 
  mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
         lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(ends_with('sigs'))


ILC_topic_results %>% dplyr::select(topic, strict_GO_sigs) %>% unnest(cols = strict_GO_sigs) %>% write_tsv('./outputs/topics_GO/ILC_strict_GO.tsv')
ILC_topic_results %>% dplyr::select(topic, lax_GO_sigs) %>% unnest(cols = lax_GO_sigs) %>%  write_tsv('./outputs/topics_GO/ILC_lax_GO.tsv')

#### ILC Gut Blood Haystack


ILC_topic_results <- 
  read_xlsx('./data/GutBloodILCs_k9_GeneModules.xlsx') %>% #pull(newP) %>% hist()
  # filter(LFC >.5) %>%
  # filter(newP < 0.05) %>% 
  group_by(topic) %>% 
  nest() %>% 
  mutate(cutoff=map_dbl(data, ~quantile(.x$Zscore, .95))) %>% 
  mutate(strict_filt=map2(data, cutoff,  ~filter(.x, Zscore > .y)), 
         lax_filt   =map(data, ~filter(.x, LFC > .5))) %>% 
  mutate(strict_GO=
           map(strict_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/ILC_GO_universe.tsv')), 
         lax_GO=
           map(lax_filt,
               ~topGO_wrapper(myInterestingGenes = .x$Gene, mapping_file = './outputs/ILC_GO_universe.tsv'))) %>% 
  mutate(strict_GO_sigs = map(strict_GO, ~filter(.x, pval < 0.05)), 
         lax_GO_sigs = map(lax_GO, ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(ends_with('sigs'))

