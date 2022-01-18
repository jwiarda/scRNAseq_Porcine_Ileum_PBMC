library(readxl)
library(tidyverse)
library(Seurat)
library(funfuns)

library(fastTopics)

system('mkdir ./outputs/official_annotations_GO/')

Bcells <- read_xlsx('./data/final_groups/B_FINALannot_OverallDE.xlsx')


B_cell_results <- 
  Bcells %>% 
  filter(avg_logFC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/Bcell_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

B_cell_results %>% write_tsv('./outputs/official_annotations_GO/B_cells.tsv')



######
Myeloid <- read_xlsx('./data/final_groups/Myeloid_FINALannot_OverallDE.xlsx')


Myeloid_results <- 
  Myeloid %>% 
  filter(avg_logFC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/Myeloid_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Myeloid_results %>% write_tsv('./outputs/official_annotations_GO/Myeloid.tsv')


###

NonImmune <- read_xlsx('./data/final_groups/NonImmune_FINALannot_OverallDE.xlsx')


NonImmune_results <- 
  NonImmune %>% 
  filter(avg_logFC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/NonImmune_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

NonImmune_results %>% write_tsv('./outputs/official_annotations_GO/NonImmune.tsv')

###

TILC <- read_xlsx('./data/final_groups/TILC_FINALannot_OverallDE (1).xlsx')

TILC_results <- 
  TILC %>% 
  filter(avg_logFC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/TILC_GO_universe.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

TILC_results %>% write_tsv('./outputs/official_annotations_GO/TILC.tsv')


###

All <- read_xlsx('./data/final_groups/AllCells_FINALannot_OverallDE.xlsx')

All_results <- 
  All %>% 
  filter(avg_logFC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = './outputs/All_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

All_results %>% write_tsv('./outputs/official_annotations_GO/All.tsv')

# gut blood TILCs

#DELETE
# # this is a little different because jayne wanted to know which genes
# # were contributing to the sig go enrichments in each cluster
# 
# 
# # probably unneccessary helper function
# subset_GO_gene_list <- 
#   function(GO_gene_list, enriched_genes_vector){
#     # subset each entry in the list to only contain the genes in the 
#     # enriched gene vector
#     map(GO_gene_list, ~ .x[.x %in% enriched_genes_vector])
#   }

gut_blood_TILC <- read_xlsx('./data/GutBloodILCs_k9_GeneModules.xlsx')

gut_blood_TILC_results <- 
  gut_blood_TILC %>% 
  mutate(cluster = 
           case_when(res.hc.clusters %in% c(3, 8)  ~ '3_8', 
                     TRUE                    ~ as.character(res.hc.clusters))) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes,
               ~topGO_wrapper(myInterestingGenes = .x,
                              mapping_file = './outputs/gut_blood_TILC_GO_universe.tsv')), 
         filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols = filt_results) 


gut_blood_TILC_results %>% 
  write_tsv('./outputs/official_annotations_GO/gut_blood_TILC_k9_GeneModules_GO.tsv')

  

# 
# contributing_genes <- 
#   gut_blood_TILC_results %>% 
#   mutate(sig_GOs=map(filt_results, ~pull(.x, GO.ID))) %>% 
#   dplyr::select(cluster, allGOs, sig_GOs, enriched_genes) %>% 
#   mutate(filt_allGOs=map2(.x = allGOs, .y=sig_GOs, .f = ~.x[names(.x) %in% .y])) %>% 
#   dplyr::select(cluster,filt_allGOs, enriched_genes) %>% 
#   mutate(SUBSET_ALLGOS=map2(.x=filt_allGOs, .y=enriched_genes, .f=subset_GO_gene_list)) %>% 
#   dplyr::select(cluster, SUBSET_ALLGOS) %>%
#   unnest_longer(col = SUBSET_ALLGOS,
#                 indices_to = 'GO.ID',
#                 values_to = 'involved_genes',
#                 simplify = T) %>% 
#   ungroup() %>%
#   mutate(involved_genes=map_chr(involved_genes, paste0, collapse = "_"))
# 
# final_results <- gut_blood_TILC_results %>% 
#   dplyr::select(cluster, filt_results) %>% 
#   unnest(cols=filt_results) %>%
#   left_join(contributing_genes)
# 
# 
# final_results %>%
#   write_tsv('./outputs/official_annotations_GO/gut_blood_TILC_k9_GeneModules_GO.tsv')
# 
# 
# # step0: pull sig GO terms from res_df
# # step1: reduce allGOs to only significant GO terms (sig_allGOs)
# # step2: 
# 
# 
# GOdat_TEST1 <- gut_blood_TILC_results$GO_obj[[1]]
# 
# GOdat_TEST1@allGenes
# 
# GOdat_TEST1@graph
# 
# 
# ## THESE GO TERMS
# filtered_res <- gut_blood_TILC_results$GO_res[[1]]
# 
# filtered_res$GO.ID
# 
# ###
# allGO = genesInTerm(GOdat_TEST1)
# 
# allGO[["GO:0006457"]][allGO[["GO:0006457"]] %in% gut_blood_TILC_results$enriched_genes[[1]]]
# 
# 
# TTTEEESSSTTT <- 
#   TEST 
# TTTEEESSSTTT %>% 
# paste0(TTTEEESSSTTT$involved_genes[[1]], collapse = "_")
# 
# 
# 
# TTTEEESSSTTT %>% 
#   dplyr::select(cluster, SUBSET_ALLGOS) %>% 
#   unnest(cols = SUBSET_ALLGOS) %>% 
# 
# enframe(TTTEEESSSTTT$SUBSET_ALLGOS)
# 
# 
# ## all true
# names(allGO) %in% filtered_res$GO.ID
# 
# all
# 
# ###
# # TODO
# # To implement this I need to have topGO_wrapper return the topGOdata object
# 
# 
# 
# 
# #
# # 
# # Hi Jean-Pierre,
# #   I assume you were using TopGO from the message header,
# # so you can do the following:
# #     # sample are your candidate genes and allgenes all other genes
# #   # you included in your analysis
# universe = factor(as.integer(allgenes %in% sample))
# names(universe) = allgenes
# #   # you created an GOdata object for TopGO
# #   # "GO" would be your annotation Gene2GOid list
# GOdata = new('topGOdata',ontology="BP",allGenes=universe,
#                  annot=annFUN.gene2GO,gene2GO=GO)
# #   # retrieve genes2GO list from the "expanded" annotation in GOdata
# allGO = genesInTerm(GOdata)
# allGO["GO:0000109"]
# #   #$`GO:0000109`
# #   #[1] "ENSG00000012061" "ENSG00000104472" "ENSG00000175595"
# SAM_ANOTATION = lapply(allGO,function(x) x[x %in% sample] )
#   # Your significant genes for GO:0051427
# SAM_ANOTATION[["GO:0051427"]]
# #   Hope it was helpful.
# # Best wishes,
# # Ricardo
# #     Jean-Pierre Desvignes wrote:
# #   Hello,
# #   I would known if it was possible to retrieve the genes related to a
# # GO term ?
# #     If yes, how to do that ?
# #     GO.ID  Term  Level  Annotated  "Significant"  Expected
# # classicFisher
# #   GO:0051427  hormone receptor binding  5  16  17  "2"  0.12  0.00583
# #   For example, what are the 2 Significant genes here ?
# #         Thanks a lot.
