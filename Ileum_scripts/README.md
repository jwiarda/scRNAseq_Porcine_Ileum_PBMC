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
| 06a_Ileum_AllCells_GenCharacterization.Rmd | Assessment of cell metrics, canonical gene query, cell lineage assignment of cell clusters, assessing cell/meta data distributions |
| 06b_Ileum_AllCells_HierarchicalClustering.Rmd | Determine relatedness of cell clusters through hierarchical clustering |
| 06c_Ileum_AllCells_OverallDGE.Rmd | Determine genes differentially expressed between cell clusters of entire dataset |
06d_Ileum_AllCells_PairwiseDGE.Rmd
06e_Ileum_AllCells_Pseudobulk_Correlations.Rmd
06f_Ileum_AllCells_Lineage_Tissue_Comps.Rmd
06g_Ileum_AllCells_Sample_Comps.Rmd
06h_Ileum_AllCells_GeneSetEnrichment.Rmd
07j_Ileum_AllCells_PredictionMapping_Visualize.Rmd
08a_Ileum_TILC_Subset_DimReduc.Rmd
08b_Ileum_TILC_HierarchicalClustering.Rmd
08c_Ileum_TILC_OverallDGE.Rmd
08d_Ileum_TILC_PredictionMapping_Visualize.Rmd
08e_Ileum_TILC_Cell_Calling_NOTDONE.Rmd
08e_Ileum_TILC_Cell_Calling.Rmd
08e_Ileum_TILC_SpatialDGE.Rmd
08h_Ileum_TILC_Cell_Calling.Rmd
09a_Ileum_Tnaive_Subset_DimReduc.Rmd
09b_Ileum_Tnaive_PredictionMapping_Visualize.Rmd
09c_Ileum_NaiveT_SpatialDGE.Rmd
10a_Ileum_CD4T_Subset_DimReduc.Rmd
10b_Ileum_CD4T_HierarchicalClustering.Rmd
10c_Ileum_CD4T_ClusteredDGE.Rmd
10d_Ileum_CD4T_PredictionMapping_Visualize.Rmd
10e_Ileum_CD4T_SpatialDGE.Rmd
11a_Ileum_GDCD8T_Subset_DimReduc.Rmd
11b_Ileum_GDCD8T_HierarchicalClustering.Rmd
11c_Ileum_GDCD8T_ClusteredDGE.Rmd
11d_Ileum_GDCD8T_PredictionMapping_Visualization.Rmd
11e_Ileum_GDCD8T_SpatialDGE.Rmd
12a_Ileum_ILC_Subset_DimReduc.Rmd
12b_Ileum_ILC_HierarchicalClustering.Rmd
12c_Ileum_ILC_ClusteredDGE.Rmd
12d_Ileum_ILC_PredictionMapping_Visualize.Rmd
12e_Ileum_ILC_SpatialDGE.Rmd
13a_Ileum_B_Subset_DimReduc.Rmd
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
