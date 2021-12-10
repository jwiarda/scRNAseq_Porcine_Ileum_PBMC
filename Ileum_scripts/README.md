Only .Rmd files for scripts are available in this GitHub repository. Knit markdown files for scripts listed, including output images and step-by-step process descriptions, are available at _______ . Additional input and output data from analyses can also be found at ______ as follows:

*
*
*

# File Key:

| File name | File description |
|-----------|:----------------:|
| 00_Ileum_B_Annotations.Rmd | Hybrid method to create final annotations for B lineage lymphocytes after considering cell clustering, topic modeling, and spatial DGE results |
| 00_Ileum_CD4T_Annotations.Rmd | Hybrid method to create final annotations for CD4 alpha beta T cells after considering cell clustering, topic modeling, and spatial DGE results |
| 00_Ileum_GDCD8_Annotation.Rmd | Hybrid method to create final annotations for gamma delta & CD8 alpha beta T cells after considering cell clustering, topic modeling, and spatial DGE results |
| 00_Ileum_ILC_Annotations.Rmd | Hybrid method to create final annotations for innate lymphoid cells after considering cell clustering, topic modeling, and spatial DGE results |
| 00_Ileum_Myeloid_Annotations.Rmd | Method to create final annotations for myeloid lineage leukocytes based on cell clustering results |
| 00_Ileum_NonImmune_Annotations.Rmd | Method to create final annotations for non-leukocytes based on cell clustering results |
| 02_Ileum_Ambient_RNA_Removal.Rmd | Calculation and removal of contaminating ambient RNA from cells |
| 03_Ileum_Gene_Cell_Filtering.Rmd | Filtering of non-expressed genes and poor quality cells | 
| 04b_Ileum_Doublet_Removal.Rmd | Removal of cell doublets |
| 05_Ileum_Normalize_Integrate_DimReduc_Cluster.Rmd | SCT and log normalization, integration, calculation of PCA, UMAP, and t-SNE dimensionality reductions, and cell clustering |
| 06a_Ileum_AllCells_GenCharacterization.Rmd | General characterization of the data, where we visualize various metrics for our cells and/or clusters and query canonical genes to assign cell lineage identities to cel clusters |
| 06b_Ileum_AllCells_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst cell clusters |
| 06c_Ileum_AllCells_OverallDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a cluster relative to the rest of the dataset |
| 06d_Ileum_AllCells_PairwiseDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a cluster relative to another cluster |
| 06e_Ileum_AllCells_Pseudobulk_Correlations.Rmd | Cluster-based pseudobulk correlation as another method to assess similarity/dissimilarity between pairwise cluster combinations |
| 06f_Ileum_AllCells_Lineage_Tissue_Comps.Rmd | Tissue composition comparisons, to determine differences in the compositions of different cells and/or transcripts across different sample types |
| 06g_Ileum_AllCells_Sample_Comps.Rmd | Sample-specific quality comparisons, where we assess the whole transcriptomic profiles from each sample |
| 06h_Ileum_AllCells_GeneSetEnrichment.Rmd | Cluster-based gene set enrichment analysis, where we calculate enrichment of signatures for sorted immune cell populations from bulk RNA-seq datasets within our single-cell dataset |
| 07a | __ |
| 07b | __ |
| 07c | __ |
| 07d | __ |
| 07e | __ |
| 07f | __ |
| 07g | __ |
| 07h | __ |
| 07i | __ |
| 07j_Ileum_AllCells_PredictionMapping_Visualize.Rmd | Visualize the results of all cell mapping and label prediction |
| 08a_Ileum_TILC_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as T/ILC lineage lymphocytes |
| 08b_Ileum_TILC_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only T/ILC cell clusters |
| 08c_Ileum_TILC_OverallDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single T/ILC cell cluster relative to the rest of the T/ILC lineage lymphocytes in the dataset |
| 08d_Ileum_TILC_PredictionMapping_Visualize.Rmd | Visualize the results of cell mapping and label prediction for only T/ILC lineage lymphocytes |
| 08e_Ileum_TILC_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only T/ILC lineage lymphocytes in multidimensional data space |
| 08f | __ |
| 08g | __ |
| 08h_Ileum_TILC_Cell_Calling.Rmd | Classification of cell clusters as CD4 alpha beta T cells, CD8 alpha beta T cells, gamma delta T cells, or ILCs based on expression of marker genes |
| 09a_Ileum_Tnaive_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as naive CD4 or CD8 alpha beta T cells |
| 09b_Ileum_Tnaive_PredictionMapping_Visualize.Rmd | Visualize the results of cell mapping and label prediction for only naive CD4 or CD8 alpha beta T cells |
| 09c_Ileum_NaiveT_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only naive CD4 or CD8 alpha beta T cells in multidimensional data space |
| 10a_Ileum_CD4T_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as non-naive CD4 alpha beta T cells |
| 10b_Ileum_CD4T_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only non-naive CD4 alpha beta T cell clusters |
| 10c_Ileum_CD4T_ClusteredDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single non-naive CD4 alpha beta T cell cluster relative to the rest of the non-naive CD4 alpha beta T cells |
| 10d_Ileum_CD4T_PredictionMapping_Visualize.Rmd | Visualize the results of cell mapping and label prediction for only non-naive CD4 alpha beta T cells |
| 10e_Ileum_CD4T_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only non-naive CD4 alpha beta T cells in multidimensional data space |
| 11a_Ileum_GDCD8T_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as gamma delta or non-naive CD8 alpha beta T cells |
| 11b_Ileum_GDCD8T_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only gamma delta/non-naive CD8 alpha beta T cell clusters |
| 11c_Ileum_GDCD8T_ClusteredDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single gamma delta/non-naive CD8 alpha beta T cell cluster relative to the rest of the gamma delta/non-naive CD8 alpha beta T cells |
| 11d_Ileum_GDCD8T_PredictionMapping_Visualization.Rmd | Visualize the results of cell mapping and label prediction for only gamma delta/non-naive CD8 alpha beta T cells |
| 11e_Ileum_GDCD8T_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only gamma detla/non-naive CD8 alpha beta T cells in multidimensional data space |
| 12a_Ileum_ILC_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as ILCs |
| 12b_Ileum_ILC_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only ILC clusters |
| 12c_Ileum_ILC_ClusteredDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single ILC cluster relative to the rest of the ILCs |
| 12d_Ileum_ILC_PredictionMapping_Visualize.Rmd |  Visualize the results of cell mapping and label prediction for only ILCs |
| 12e_Ileum_ILC_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only ILCs in multidimensional data space |
| 13a_Ileum_B_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as B lineage lymphocytes |
13b_Ileum_B_HierarchicalClustering.Rmd
13c_Ileum_B_OverallDGE.Rmd
13d_Ileum_B_PredictionMapping_Visualization.Rmd
13e_Ileum_B_SpatialDGE.Rmd
14a_Ileum_Myeloid_Subset_DimReduc.Rmd
14b_Ileum_Myeloid_HierarchicalClustering.Rmd
14c_Ileum_Myeloid_OverallDGE.Rmd
14d_Ileum_Myeloid_PredictionMapping_Visualize.Rmd
14e_Ileum_Myeloid_SpatialDGE.Rmd
15a_Ileum_NonImmune_Subset_DimReduc.Rmd
15b_Ileum_NonImmune_HierarchicalClustering.Rmd
15c_Ileum_NonImmune_OverallDGE.Rmd
15d_Ileum_NonImmune_PredictionMapping_Visualization.Rmd
15e_Ileum_NonImmune_SpatialDGE.Rmd
16a_Ileum_AllCellsAnnot_GenCharacterization.Rmd
16b_Ileum_AllCellsAnnot_HierarchichalClustering.Rmd
16c_Ileum_AllCellsAnnot_OverallDGE.Rmd
16d_Ileum_AllCellsAnnot_PairwiseDGE.Rmd
16e_Ileum_AllCellsAnnot_Pseudobulk_Correlations.Rmd
16f_Ileum_AllCellsAnnot_VisualizeMapPredict.Rmd
17a_Ileum_TILCAnnot_GenCharacterization.Rmd
17b_Ileum_TILCAnnot_HierarchicalClustering.Rmd
17c_Ileum_TILCAnnot_OverallDGE.Rmd
18a_Ileum_BAnnot_GenCharacterization.Rmd
18b_Ileum_BAnnot_HierarchicalClustering.Rmd
18c_Ileum_BAnnot_OverallDGE.Rmd
19a_Ileum_MyeloidAnnot_GenCharacterization.Rmd
19b_Ileum_MyeloidAnnot_HierarchicalClustering.Rmd
19c_Ileum_MyeloidAnnot_OverallDGE.Rmd
20a_Ileum_NonImmuneAnnot_GenCharacterization.Rmd
20b_Ileum_NonImmuneAnnot_OverallDGE.Rmd
21a_Ileum_GDT_Subset_DimReduc.Rmd
21a_Ileum_PPcomp_Subset_DimReduc.Rmd
21b_Ileum_PPcomp_DA.Rmd