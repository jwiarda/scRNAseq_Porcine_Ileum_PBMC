Only .Rmd files for scripts are available in this GitHub repository. Knit markdown files for scripts listed, including output images and step-by-step process descriptions, are available at _______ . Additional input and output data from analyses can also be found at ______. Cell annotations created in 16a-f and used in subseqent scripts were made by also incorporating results from topic modeling and GO analyses included in a separate subdirectory.

# File Key:

| File name | File description |
|-----------|:----------------:|
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
| 07a_Ileum_AllCells_CreateReference_SsPBMC.Rmd | Create reference dataset from previously-published and annotated porcine PBMC scRNA-seq data |
| 07b_Ileum_AllCells_CreateReference_HsIleum.Rmd | Create reference dataset from previously-published and annotated human ileum scRNA-seq data |
| 07c_Ileum_AllCells_CreateReference_MmIleum.Rmd | Create reference dataset from previously-published and annotated murine ileum scRNA-seq data |
| 07d_Ileum_AllCells_CreateQuery_SsIleum.Rmd | Convert porcine ileum scRNA-seq data into query datasets for comparison to previously-published porcine PBMCs |
| 07e_Ileum_AllCells_CreateQuery_SsIleum_Humanized.Rmd | Convert porcine ileum scRNA-seq data into query datasets for comparison to previously-published human ileum |
| 07f_Ileum_AllCells_CreateQuery_SsIleum_Murinized.Rmd | Convert porcine ileum scRNA-seq data into query datasets for comparison to previously-published murine ileum |
| 07g_Ileum_AllCells_PredictionMapping_RefSsPBMC.Rmd | Perform reference-based mapping and label transfer of porcine ileum query data onto porcine PBMC reference data |
| 07h_Ileum_AllCells_PredictionMapping_RefHsIleum.Rmd | Perform reference-based mapping and label transfer of porcine ileum query data onto human ileum reference data |
| 07i_Ileum_AllCells_PredictionMapping_RefMmIleum.Rmd | Perform reference-based mapping and label transfer of porcine ileum query data onto murine ileum reference data |
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
| 13b_Ileum_B_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only B lineage lymphocyte clusters |
| 13c_Ileum_B_OverallDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single B lineage lymphocyte cluster relative to the rest of the B lineage lymphocytes |
| 13d_Ileum_B_PredictionMapping_Visualization.Rmd | Visualize the results of cell mapping and label prediction for only B lineage lymphocytes |
| 13e_Ileum_B_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only B lineage lymphocytes in multidimensional data space |
| 14a_Ileum_Myeloid_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as myeloid lineage leukocytes |
| 14b_Ileum_Myeloid_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only myeloid lineage leukocyte clusters |
| 14c_Ileum_Myeloid_OverallDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single myeloid lineage leukocyte cluster relative to the rest of the myeloid lineage leukocytes |
| 14d_Ileum_Myeloid_PredictionMapping_Visualize.Rmd | Visualize the results of cell mapping and label prediction for only myeloid lineage leukocytes |
| 14e_Ileum_Myeloid_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only myeloid lineage leukocytes in multidimensional data space |
| 15a_Ileum_NonImmune_Subset_DimReduc.Rmd | Subset and re-process cells that were identified as non-leukocytes |
| 15b_Ileum_NonImmune_HierarchicalClustering.Rmd | Hierarchical clustering of the data to determine phylogenetic relationship amongst only non-leukocyte clusters |
| 15c_Ileum_NonImmune_OverallDGE.Rmd | Cluster-based differential gene expression analysis to determine genes differentially expressed in a single non-leukocyte cluster relative to the rest of the non-leukocytes |
| 15d_Ileum_NonImmune_PredictionMapping_Visualization.Rmd | Visualize the results of cell mapping and label prediction for only non-leukocytes |
| 15e_Ileum_NonImmune_SpatialDGE.Rmd | Cluster-independent calculation of differentially expressed genes and gene modules for only non-leukocytes in multidimensional data space |
| 16a_Ileum_B_Annotations.Rmd | Hybrid method to create final annotations for B lineage lymphocytes after considering cell clustering, topic modeling, and spatial DGE results |
| 16b_Ileum_CD4T_Annotations.Rmd | Hybrid method to create final annotations for CD4 alpha beta T cells after considering cell clustering, topic modeling, and spatial DGE results |
| 16c_Ileum_GDCD8_Annotation.Rmd | Hybrid method to create final annotations for gamma delta & CD8 alpha beta T cells after considering cell clustering, topic modeling, and spatial DGE results |
| 16d_Ileum_ILC_Annotations.Rmd | Hybrid method to create final annotations for innate lymphoid cells after considering cell clustering, topic modeling, and spatial DGE results |
| 16e_Ileum_Myeloid_Annotations.Rmd | Method to create final annotations for myeloid lineage leukocytes based on cell clustering results |
| 16f_Ileum_NonImmune_Annotations.Rmd | Method to create final annotations for non-leukocytes based on cell clustering results |
| 17a_Ileum_AllCellsAnnot_GenCharacterization.Rmd | General characterization of the data similar to 06a but with official cell type/lineage annotations instead of cell clusters |
| 17b_Ileum_AllCellsAnnot_HierarchichalClustering.Rmd | Hierarchical clustering of the data similar to 06b but with official cell type/lineage annotations instead of cell clusters |
| 17c_Ileum_AllCellsAnnot_OverallDGE.Rmd | Differential gene expression analysis similar to 06c but with official cell type/lineage annotations instead of cell clusters |
| 17d_Ileum_AllCellsAnnot_PairwiseDGE.Rmd | Differential gene expression analysis similar to 06d but with official cell type/lineage annotations instead of cell clusters |
| 17e_Ileum_AllCellsAnnot_Pseudobulk_Correlations.Rmd | Pseudobulk correlation analysis similar to 06e but with official cell type/lineage annotations instead of cell clusters |
| 17f_Ileum_AllCellsAnnot_VisualizeMapPredict.Rmd | Visualizing mapping and prediction scores, now also considering official cell type/lineage annotations |
| 18a_Ileum_TILCAnnot_GenCharacterization.Rmd | Add final cell type annotations to T/ILC lineage lymphocytes and visualize results |
| 18b_Ileum_TILCAnnot_HierarchicalClustering.Rmd |  Hierarchical clustering similar to 16b but for only T/ILC lineage lymphocytes with official cell type annotations |
| 18c_Ileum_TILCAnnot_OverallDGE.Rmd | Differential gene expression analysis similar to 16c but for only T/ILC lineage lymphocytes with official cell type annotations |
| 19a_Ileum_BAnnot_GenCharacterization.Rmd | Add final cell type annotations to B lineage lymphocytes and visualize results |
| 19b_Ileum_BAnnot_HierarchicalClustering.Rmd | Hierarchical clustering similar to 16b but for only B lineage lymphocytes with official cell type annotations |
| 19c_Ileum_BAnnot_OverallDGE.Rmd | Differential gene expression analysis similar to 16c but for only B lineage lymphocytes with official cell type annotations |
| 20a_Ileum_MyeloidAnnot_GenCharacterization.Rmd | Add final cell type annotations to myeloid lineage lymphocytes and visualize results |
| 20b_Ileum_MyeloidAnnot_HierarchicalClustering.Rmd | Hierarchical clustering similar to 16b but for only myeloid lineage lymphocytes with official cell type annotations |
| 20c_Ileum_MyeloidAnnot_OverallDGE.Rmd | Differential gene expression analysis similar to 16c but for only myeloid lineage lymphocytes with official cell type annotations |
| 21a_Ileum_NonImmuneAnnot_GenCharacterization.Rmd | Add final cell type annotations to non-leukocytes and visualize results |
| 21b_Ileum_NonImmuneAnnot_OverallDGE.Rmd | Differential gene expression analysis similar to 16c but for only non-leukocytes with official cell type annotations |
| 22a_Ileum_PPcomp_Subset_DimReduc.Rmd | Create a data subset consisting of only cells from PP and non-PP samples (excluding whole ileum samples) |
| 22b_Ileum_PPcomp_DA.Rmd | Perform differential abundance analysis between cells derived from PP and non-PP samples (from data subset generated in 21a) |
