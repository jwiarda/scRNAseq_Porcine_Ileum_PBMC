library(tidyverse)
library(biomaRt)  # used to map GOterms to Ensembl IDs

# Order #
# this should be an early script, probably after data download.



GO_background_generator <- 
  function(GO_term_universe, seurat_obj, output_path){
    # gene list is enriched genes, 
    # GO_term_universe is already prepared, needs to be filtered
    # seurat obj will be used to filter GO_term_universe
    # browser()
    # need to make sure all genes in seurat object are detected
    counts <- seurat_obj@assays$RNA@counts
    if (!all(rowSums(counts) > 0)){
      print('some genes not detected, filtering genes with 0 counts')
      counts <- counts[rowSums(counts) > 0,]
    }
    
    all_GO <- read_tsv('./outputs/All_gene_to_GO.tsv')
    new_universe <-  all_GO %>% 
      filter(Name %in% rownames(counts)) %>% 
      write_tsv(output_path)
    print(paste('number of genes in all universe',nrow(all_GO)))
    print(paste('number of genes in new universe',nrow(new_universe)))
    
  }




# wrapper function that makes using topGO easier (I think)

# you can also access this function by loading my R package 'funfuns'
# that way you don't need to read in this function every time
# remotes::install_github('jtrachsel/funfuns')
# library(funfuns)


topGO_wrapper <- function(myInterestingGenes, #vector
                          mapping_file,       # two column file
                          ont='BP',
                          algor = 'elim',
                          statistic='Fisher',
                          nodeSize=10, 
                          return_GOdata=F){

  require(topGO)
  # browser()
  geneID2GO <- readMappings(mapping_file)
  geneNames <- names(geneID2GO)

  # Get the list of genes of interest
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames

  #initialize topGOdata object
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO,
                nodeSize=nodeSize)
  
  # contstruct a tibble that maps interesting genes to GO terms
  # this can add a decent amount of time to the function call...
  interesting_genes_in_GOs <- 
    genesInTerm(GOdata) %>% 
    enframe(name='GO.ID', 
            value='genes_in_GO') %>% 
    mutate(involved_genes=map(.x=genes_in_GO, ~ .x[.x %in% myInterestingGenes]), 
           involved_genes=map_chr(.x=involved_genes, ~paste(.x, collapse = '_'))) %>% 
    dplyr::select(GO.ID, involved_genes) 
  
  # Run topGO test
  resultTopGO.elim <- runTest(GOdata, algorithm = algor, statistic = statistic )
  allRes <- GenTable(GOdata, pval = resultTopGO.elim,
                     orderBy = "pval",
                     topNodes = length(GOdata@graph@nodes), #include all nodes
                     numChar=1000)
  
  # clean up results and add in extra info
  allRes <- allRes %>%
    mutate(ont=ifelse(ont=='BP', 'Biological Process',
                      ifelse(ont=='MF', 'Molecular Function', "Cellular Component"))) %>%
    mutate(GO_aspect = ont,
           algorithm = algor,
           statistic = statistic) %>%
    dplyr::select(-ont) %>% 
    left_join(interesting_genes_in_GOs)
  
  if (return_GOdata == TRUE){
    return(list(allRes, GOdata))
  } else {
    return(allRes)
  }
  

}








#### Get ensemble gene_IDs for pig genome ###
# This section uses the internet to map these IDs, so sometimes it is very slow
# in fact it is not working right now.  I think the USDA's internet tubes get all 
# blocked up during business hours, what with all the video calls

# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("sscrofa_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id'))

# examine result
head(EG2GO,15)

# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

### format mapping file for use with topGO wrapper function
EG2GO <- 
  EG2GO %>%
  group_by(ensembl_gene_id) %>% 
  summarise(GO=paste(go_id, sep = ' ', collapse = ',')) %>% 
  transmute(EnsemblID=ensembl_gene_id, 
            GO=GO) 


#### Read in IleumAtlas data

# to map gene 'Name' to Ensembl gene_IDs
gene_IDs <- readxl::read_xlsx('UnfilteredGeneInfo.xlsx')

library(Seurat)
Ileum_all <- readRDS(file = 'data/IleumAtlasAll.rds')

detected_genes <- rownames(Ileum_all@assays$RNA@counts)

#

#



# all GO terms detected in this Ileum Atlas data

GO_gene_universe <- 
  gene_IDs %>% 
  filter(Name %in% detected_genes) %>% 
  left_join(EG2GO) %>% 
  filter(!is.na(GO))


GO_gene_universe %>%
  dplyr::select(Name, GO) %>% 
  write_tsv('./outputs/All_gene_to_GO.tsv')



# Ileum_all <- read_rds('data/IleumAtlasAll.rds')
# GO universe is allready set up for this.  
#

## setup background GO terms for each cellular subset.

Ileum_Bcells <- read_rds('data/Ileum_Bonly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_Bcells, 
                        output_path = './outputs/Bcell_GO_universe.tsv')


Ileum_CD4T <- read_rds('data/Ileum_CD4Tonly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_CD4T, 
                        output_path = './outputs/CD4_GO_universe.tsv')


Ileum_gdCD8T <- read_rds('data/Ileum_gdCD8Tonly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_gdCD8T, 
                        output_path = './outputs/gdCD8T_GO_universe.tsv')


Ileum_ILC <- read_rds('data/Ileum_ILConly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_ILC, 
                        output_path = './outputs/ILC_GO_universe.tsv')


Ileum_Myeloid <- read_rds('data/Ileum_MyeloidOnly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_Myeloid, 
                        output_path = './outputs/Myeloid_GO_universe.tsv')


Ileum_NonImmune <- read_rds('data/Ileum_NonImmuneOnly.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = Ileum_NonImmune, 
                        output_path = './outputs/NonImmune_GO_universe.tsv')


#

TILC <- read_rds('data/Ileum_TILConly.rds')

GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = TILC, 
                        output_path = './outputs/TILC_GO_universe.tsv')

#

gut_blood_TILC <- read_rds('data/GutBlood_IntegratedILCs.rds')


GO_background_generator(GO_term_universe = './outputs/All_gene_to_GO.tsv',
                        seurat_obj = gut_blood_TILC, 
                        output_path = './outputs/gut_blood_TILC_GO_universe.tsv')

