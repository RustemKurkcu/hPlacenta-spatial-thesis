hPlacenta-architecture
├── config
│   ├── celltype_refinement_map.csv
│   ├── config.R
│   ├── config_01_gene_signatures.R
│   ├── gene_sets.R
│   ├── ligand_receptor_pairs.csv
│   └── vulnerability_config.R
├── data
│   ├── processed
│   │   ├── multiome_rna_seurat.rds
│   │   ├── seurat_spatial_prepared.qs
│   │   ├── slidetags_mapped_to_multiome.rds
│   │   └── starmap_spatial_raw_plus_imputed_seurat.rds
│   └── raw
│       ├── Broad_SCP2601human-placenta-architecture
│       │   ├── barcodes-slidetag-atac-js40-filt.tsv.gz
│       │   ├── barcodes_atac.tsv
│       │   ├── barcodes_rna.tsv
│       │   ├── barcodes_rna2.tsv
│       │   ├── genes_rna.tsv
│       │   ├── genes_rna2.tsv
│       │   ├── humanplacenta_cluster.csv
│       │   ├── humanplacenta_expression.csv.gz
│       │   ├── humanplacenta_expression_raw.csv.gz
│       │   ├── humanplacenta_genescore.csv.gz
│       │   ├── humanplacenta_motifzscore.csv.gz
│       │   ├── humanplacenta_spatial.csv
│       │   ├── manifest.txt
│       │   ├── matrix.mtx
│       │   ├── metadata.csv
│       │   ├── OriginalDataFromSingleCellAndSpatialPaper.zip
│       │   └── peakAnnotations_atac.tsv
│       ├── zenodo_spatial
│       │   ├── STARmap-ISH_sample_W11_G1.png
│       │   ├── STARmap-ISH_sample_W11_G2.png
│       │   ├── STARmap-ISH_sample_W11_G3.png
│       │   ├── STARmap-ISH_sample_W11_G4.png
│       │   ├── STARmap-ISH_sample_W11_G5.png
│       │   ├── STARmap-ISH_sample_W11_G6.png
│       │   ├── STARmap-ISH_sample_W7-2_G1.png
│       │   ├── STARmap-ISH_sample_W7-2_G3.png
│       │   ├── STARmap-ISH_sample_W7-2_G4.png
│       │   ├── STARmap-ISH_sample_W7-2_G5.png
│       │   ├── STARmap-ISH_sample_W7-2_G6.png
│       │   ├── STARmap-ISH_sample_W8-2_G1.png
│       │   ├── STARmap-ISH_sample_W8-2_G2.png
│       │   ├── STARmap-ISH_sample_W8-2_G3.png
│       │   ├── STARmap-ISH_sample_W8-2_G4.png
│       │   ├── STARmap-ISH_sample_W8-2_G5.png
│       │   ├── STARmap-ISS_sample_W11_imputed_expression.csv
│       │   ├── STARmap-ISS_sample_W11_raw_expression.csv
│       │   ├── STARmap-ISS_sample_W11_spots_metadata.csv
│       │   ├── STARmap-ISS_sample_W7_cell_metadata.csv
│       │   ├── STARmap-ISS_sample_W7_imputed_DORC.csv
│       │   ├── STARmap-ISS_sample_W7_imputed_expression.csv
│       │   ├── STARmap-ISS_sample_W7_imputed_gene_activity.csv
│       │   ├── STARmap-ISS_sample_W7_raw_expression.csv
│       │   ├── STARmap-ISS_sample_W7_spots_metadata.csv
│       │   ├── STARmap-ISS_sample_W8-2_cell_metadata.csv
│       │   ├── STARmap-ISS_sample_W8-2_imputed_DORC.csv
│       │   ├── STARmap-ISS_sample_W8-2_imputed_expression.csv
│       │   ├── STARmap-ISS_sample_W8-2_imputed_gene_activity.csv
│       │   ├── STARmap-ISS_sample_W8-2_raw_expression.csv
│       │   ├── STARmap-ISS_sample_W8-2_spots_metadata.csv
│       │   ├── STARmap-ISS_sample_W9_cell_metadata.csv
│       │   ├── STARmap-ISS_sample_W9_imputed_DORC.csv
│       │   ├── STARmap-ISS_sample_W9_imputed_expression.csv
│       │   ├── STARmap-ISS_sample_W9_imputed_gene_activity.csv
│       │   ├── STARmap-ISS_sample_W9_raw_expression.csv
│       │   └── STARmap-ISS_sample_W9_spots_metadata.csv
│       ├── Dataset_Summary.txt
│       └── STARmap-ISS_sample_W7_raw_expression.csv
├── docs
│   ├── Additional
│   │   ├── ANALYSIS_REPORT.md
│   │   ├── Re-OrganizeScripts
│   │   ├── README.md
│   │   └── SOURCE_MANIFEST.md
│   ├── 00_COUNTS_AND_PLATFORMS.md
│   ├── 01_ANALYSIS_GUARDRAILS.md
│   ├── 02_FIGURE_PLAN_AND_METHODS_SNIPPETS.md
│   ├── 03_RUN_ORDER.md
│   ├── 04_GREENBAUM_BRIDGE.md
│   ├── 05_METHODS_WRITEUP_TEMPLATE.md
│   ├── 06C_cellchat_optional.R
│   ├── ANALYSIS_REPORT.md
│   ├── directory_tree.txt
│   ├── Force_Push.wsl
│   ├── SOURCE_MANIFEST.bib
│   └── SOURCE_MANIFEST.md
├── jian-shu-lab
│   ├── chromatin-potential
│   │   ├── 1. calculate_dorc.ipynb
│   │   └── 2. dorc_analysis.ipynb
│   ├── GWAS_analysis
│   │   ├── ALL_bedgraph_calc_tissue2.R
│   │   ├── build_module_annotations.R
│   │   ├── build_module_annotations1.sh
│   │   ├── build_module_annotations2.sh
│   │   ├── clean_bedgraphs.sh
│   │   ├── clean_bedgraphs_placentas2g.sh
│   │   ├── create_annot_from_bedgraph.sh
│   │   ├── create_annot_from_bedgraph_placentas2g.sh
│   │   ├── get_sd_annot.R
│   │   ├── get_sd_annot1.sh
│   │   ├── get_sd_annot1_placentas2g.sh
│   │   ├── gwas_hits_celltypes_gene.R
│   │   ├── ldsc_mega.sh
│   │   ├── ldsc_mega_placentas2g.sh
│   │   ├── ldsc_postprocess.R
│   │   ├── ldsc_postprocess1.sh
│   │   ├── ldsc_postprocess1_placentas2g.sh
│   │   ├── ldsc_reg.sh
│   │   ├── ldsc_reg_placentas2g.sh
│   │   ├── make_annot_combine_from_bedgraph.py
│   │   ├── placenta_get_bedgraphs.R
│   │   ├── placenta_magma_sets.R
│   │   ├── placenta_magma_sets.R~
│   │   ├── placenta_s2g_generate.R
│   │   ├── placenta_traits_enrichment_eachtrait.R
│   │   ├── placenta_traits_enrichments.R
│   │   └── process_modules_workflow.R
│   ├── multiome-analysis
│   │   ├── cluster_labels_res03.tsv
│   │   ├── cluster_labels_res03_EVT.tsv
│   │   ├── cluster_labels_res03_mac.tsv
│   │   └── preprocessing-clustering-peakanalysis.R
│   ├── STARmap_ISH_data_processing
│   │   ├── README.md
│   │   ├── round_alignment.py
│   │   ├── stitch_gene_channels.py
│   │   ├── stitch_nuclei_images.py
│   │   └── visualize_multiplexed_ISH.py
│   ├── STARmap_ISS_data_processing
│   │   ├── integration_label_transfer
│   │   │   ├── build_reference_data.R
│   │   │   └── seurat_reference_mapping.R
│   │   ├── visualization
│   │   │   ├── visualize_cell_types_boundaries.py
│   │   │   ├── visualize_markers_whole_slice.py
│   │   │   └── visualize_markers_zoom-in.py
│   │   ├── cell_segmentation.py
│   │   ├── create_starfish_format_data.py
│   │   ├── image_registration.py
│   │   ├── infer_cell_interaction.R
│   │   ├── README.md
│   │   ├── run_starfish_pipeline.py
│   │   ├── stitch_decoded_spots.py
│   │   ├── stitch_nuclei_image.py
│   │   └── stitch_spots_after_segmentation.py
│   └── STARmap_ISS_Imputation
│       ├── ATAC_profiles
│       │   ├── final_mapping_dorc.py
│       │   └── final_mapping_geneactivity.py
│       ├── final_mapping.py
│       ├── intermediate_mapping.py
│       ├── README.md
│       ├── visualize_benchmark.py
│       └── visualize_imputation_results.py
├── output
│   ├── figures
│   │   ├── 05_spatial
│   │   │   ├── slidetags
│   │   │   │   ├── W11
│   │   │   │   │   ├── slidetags_W11_3panel_global.png
│   │   │   │   │   ├── slidetags_W11_3panel_global_caption.txt
│   │   │   │   │   ├── slidetags_W11_3panel_global_white_red.png
│   │   │   │   │   ├── slidetags_W11_3panel_global_white_red_caption.txt
│   │   │   │   │   ├── slidetags_W11_permissiveness_global_week_11.png
│   │   │   │   │   ├── slidetags_W11_permissiveness_global_week_11_caption.txt
│   │   │   │   │   ├── slidetags_W11_protected_and_diff.png
│   │   │   │   │   └── slidetags_W11_protected_and_diff_caption.txt
│   │   │   │   ├── W8-2
│   │   │   │   │   ├── slidetags_W8-2_3panel_global.png
│   │   │   │   │   ├── slidetags_W8-2_3panel_global_caption.txt
│   │   │   │   │   ├── slidetags_W8-2_3panel_global_white_red.png
│   │   │   │   │   ├── slidetags_W8-2_3panel_global_white_red_caption.txt
│   │   │   │   │   ├── slidetags_W8-2_permissiveness_global_week_8.png
│   │   │   │   │   ├── slidetags_W8-2_permissiveness_global_week_8_caption.txt
│   │   │   │   │   ├── slidetags_W8-2_protected_and_diff.png
│   │   │   │   │   └── slidetags_W8-2_protected_and_diff_caption.txt
│   │   │   │   └── W9
│   │   │   │       ├── slidetags_W9_3panel_global.png
│   │   │   │       ├── slidetags_W9_3panel_global_caption.txt
│   │   │   │       ├── slidetags_W9_3panel_global_white_red.png
│   │   │   │       ├── slidetags_W9_3panel_global_white_red_caption.txt
│   │   │   │       ├── slidetags_W9_permissiveness_global_week_9.png
│   │   │   │       ├── slidetags_W9_permissiveness_global_week_9_caption.txt
│   │   │   │       ├── slidetags_W9_protected_and_diff.png
│   │   │   │       └── slidetags_W9_protected_and_diff_caption.txt
│   │   │   ├── slidetags_harmonized
│   │   │   │   ├── W11
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel.pdf
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel.png
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel_white_red.pdf
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel_white_red.png
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel_with_celltypes.pdf
│   │   │   │   │   ├── slidetags_harmonized_W11_3panel_with_celltypes.png
│   │   │   │   │   ├── slidetags_harmonized_W11_highlighted_celltypes.png
│   │   │   │   │   └── slidetags_harmonized_W11_permissiveness_week_11.png
│   │   │   │   └── W8-2
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel.pdf
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel.png
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel_white_red.pdf
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel_white_red.png
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel_with_celltypes.pdf
│   │   │   │       ├── slidetags_harmonized_W8-2_3panel_with_celltypes.png
│   │   │   │       ├── slidetags_harmonized_W8-2_highlighted_celltypes.png
│   │   │   │       └── slidetags_harmonized_W8-2_permissiveness_week_8.png
│   │   │   ├── slidetags_with_permissiveness
│   │   │   │   ├── W11
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel.pdf
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel.png
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_white_red.pdf
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_white_red.png
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_with_celltypes.pdf
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_with_celltypes.png
│   │   │   │   │   ├── slidetags_with_permissiveness_W11_highlighted_celltypes.png
│   │   │   │   │   └── slidetags_with_permissiveness_W11_permissiveness_week_11.png
│   │   │   │   └── W8-2
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel.pdf
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel.png
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_white_red.pdf
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_white_red.png
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_with_celltypes.pdf
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_with_celltypes.png
│   │   │   │       ├── slidetags_with_permissiveness_W8-2_highlighted_celltypes.png
│   │   │   │       └── slidetags_with_permissiveness_W8-2_permissiveness_week_8.png
│   │   │   ├── starmap
│   │   │   │   ├── W8-2
│   │   │   │   │   ├── starmap_W8-2_3panel_global.png
│   │   │   │   │   ├── starmap_W8-2_3panel_global_caption.txt
│   │   │   │   │   ├── starmap_W8-2_3panel_global_white_red.png
│   │   │   │   │   ├── starmap_W8-2_3panel_global_white_red_caption.txt
│   │   │   │   │   ├── starmap_W8-2_permissiveness_global_week_8.png
│   │   │   │   │   ├── starmap_W8-2_permissiveness_global_week_8_caption.txt
│   │   │   │   │   ├── starmap_W8-2_protected_and_diff.png
│   │   │   │   │   └── starmap_W8-2_protected_and_diff_caption.txt
│   │   │   │   └── W9
│   │   │   │       ├── starmap_W9_3panel_global.png
│   │   │   │       ├── starmap_W9_3panel_global_caption.txt
│   │   │   │       ├── starmap_W9_3panel_global_white_red.png
│   │   │   │       ├── starmap_W9_3panel_global_white_red_caption.txt
│   │   │   │       ├── starmap_W9_permissiveness_global_week_9.png
│   │   │   │       ├── starmap_W9_permissiveness_global_week_9_caption.txt
│   │   │   │       ├── starmap_W9_protected_and_diff.png
│   │   │   │       └── starmap_W9_protected_and_diff_caption.txt
│   │   │   ├── starmap_harmonized
│   │   │   │   └── W8-2
│   │   │   │       ├── starmap_harmonized_W8-2_3panel.pdf
│   │   │   │       ├── starmap_harmonized_W8-2_3panel.png
│   │   │   │       ├── starmap_harmonized_W8-2_3panel_white_red.pdf
│   │   │   │       ├── starmap_harmonized_W8-2_3panel_white_red.png
│   │   │   │       ├── starmap_harmonized_W8-2_3panel_with_celltypes.pdf
│   │   │   │       ├── starmap_harmonized_W8-2_3panel_with_celltypes.png
│   │   │   │       ├── starmap_harmonized_W8-2_highlighted_celltypes.png
│   │   │   │       └── starmap_harmonized_W8-2_permissiveness_week_8.png
│   │   │   ├── starmap_with_permissiveness
│   │   │   │   └── W8-2
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel.pdf
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel.png
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel_white_red.pdf
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel_white_red.png
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel_with_celltypes.pdf
│   │   │   │       ├── starmap_with_permissiveness_W8-2_3panel_with_celltypes.png
│   │   │   │       ├── starmap_with_permissiveness_W8-2_highlighted_celltypes.png
│   │   │   │       └── starmap_with_permissiveness_W8-2_permissiveness_week_8.png
│   │   │   ├── README_plot_interpretation
│   │   │   ├── README_plot_interpretation.txt
│   │   │   ├── SlideTags_hotspot_vs_coldspot_volcano.png
│   │   │   └── STARmap_hotspot_vs_coldspot_volcano.png
│   │   ├── New folder
│   │   │   ├── 05_spatial
│   │   │   │   ├── slidetags_harmonized
│   │   │   │   │   ├── W11
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel.pdf
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel.png
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel_white_red.pdf
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel_white_red.png
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel_with_celltypes.pdf
│   │   │   │   │   │   ├── slidetags_harmonized_W11_3panel_with_celltypes.png
│   │   │   │   │   │   ├── slidetags_harmonized_W11_highlighted_celltypes.png
│   │   │   │   │   │   └── slidetags_harmonized_W11_permissiveness_week_11.png
│   │   │   │   │   └── W8-2
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel.pdf
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel.png
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel_white_red.pdf
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel_white_red.png
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel_with_celltypes.pdf
│   │   │   │   │       ├── slidetags_harmonized_W8-2_3panel_with_celltypes.png
│   │   │   │   │       ├── slidetags_harmonized_W8-2_highlighted_celltypes.png
│   │   │   │   │       └── slidetags_harmonized_W8-2_permissiveness_week_8.png
│   │   │   │   ├── slidetags_with_permissiveness
│   │   │   │   │   ├── W11
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel.pdf
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel.png
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_white_red.pdf
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_white_red.png
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_with_celltypes.pdf
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_3panel_with_celltypes.png
│   │   │   │   │   │   ├── slidetags_with_permissiveness_W11_highlighted_celltypes.png
│   │   │   │   │   │   └── slidetags_with_permissiveness_W11_permissiveness_week_11.png
│   │   │   │   │   └── W8-2
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel.pdf
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel.png
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_white_red.pdf
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_white_red.png
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_with_celltypes.pdf
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_3panel_with_celltypes.png
│   │   │   │   │       ├── slidetags_with_permissiveness_W8-2_highlighted_celltypes.png
│   │   │   │   │       └── slidetags_with_permissiveness_W8-2_permissiveness_week_8.png
│   │   │   │   ├── starmap_harmonized
│   │   │   │   │   └── W8-2
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel.pdf
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel.png
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel_white_red.pdf
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel_white_red.png
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel_with_celltypes.pdf
│   │   │   │   │       ├── starmap_harmonized_W8-2_3panel_with_celltypes.png
│   │   │   │   │       ├── starmap_harmonized_W8-2_highlighted_celltypes.png
│   │   │   │   │       └── starmap_harmonized_W8-2_permissiveness_week_8.png
│   │   │   │   └── starmap_with_permissiveness
│   │   │   │       └── W8-2
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel.pdf
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel.png
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel_white_red.pdf
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel_white_red.png
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel_with_celltypes.pdf
│   │   │   │           ├── starmap_with_permissiveness_W8-2_3panel_with_celltypes.png
│   │   │   │           ├── starmap_with_permissiveness_W8-2_highlighted_celltypes.png
│   │   │   │           └── starmap_with_permissiveness_W8-2_permissiveness_week_8.png
│   │   │   ├── celltype_proportions_by_week_and_version.png
│   │   │   ├── gene_coordination_heatmap.png
│   │   │   ├── Multiome_immune_subtype_umap.png
│   │   │   ├── multiome_QC_violin.png
│   │   │   ├── multiome_reference_tsne_celltype.png
│   │   │   ├── multiome_reference_umap_celltype.png
│   │   │   ├── Slide-tags_niche_center_celltype_composition.png
│   │   │   ├── Slide-tags_niche_geneset_heatmap.png
│   │   │   ├── Slide-tags_niche_week_composition.png
│   │   │   ├── SlideTags_immune_subtype_umap.png
│   │   │   ├── SlideTags_permissiveness_top10pct_composition.png
│   │   │   ├── SlideTags_permissiveness_week_11.png
│   │   │   ├── SlideTags_permissiveness_week_8.png
│   │   │   ├── SlideTags_permissiveness_week_9.png
│   │   │   ├── slidetags_spatial_celltype.png
│   │   │   ├── slidetags_spatial_celltype_facet_by_sample.png
│   │   │   ├── slidetags_spatial_celltype_W11.png
│   │   │   ├── slidetags_spatial_celltype_W8-2.png
│   │   │   ├── slidetags_spatial_celltype_W9.png
│   │   │   ├── SlideTags_spatial_lr_top20.png
│   │   │   ├── SlideTags_spatial_week_11.png
│   │   │   ├── SlideTags_spatial_week_8.png
│   │   │   ├── SlideTags_spatial_week_9.png
│   │   │   ├── slidetags_umap_author.png
│   │   │   ├── slidetags_umap_confidence.png
│   │   │   ├── slidetags_umap_predicted.png
│   │   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_11_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_8_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap.png
│   │   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_clustered.png
│   │   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_fixed.png
│   │   │   ├── SlideTags_week_9_neighbor_z_heatmap.png
│   │   │   ├── STARmap_immune_subtype_umap.png
│   │   │   ├── STARmap_niche_center_celltype_composition.png
│   │   │   ├── STARmap_niche_clusters_spatial.png
│   │   │   ├── STARmap_niche_geneset_heatmap.png
│   │   │   ├── STARmap_niche_week_composition.png
│   │   │   ├── STARmap_permissiveness_top10pct_composition.png
│   │   │   ├── STARmap_permissiveness_week_8.png
│   │   │   ├── STARmap_permissiveness_week_9.png
│   │   │   ├── starmap_spatial_celltype.png
│   │   │   ├── starmap_spatial_celltype_8.png
│   │   │   ├── starmap_spatial_celltype_9.png
│   │   │   ├── starmap_spatial_celltype_facet_by_timepoint.png
│   │   │   ├── STARmap_spatial_lr_top20.png
│   │   │   ├── STARmap_spatial_week_8.png
│   │   │   ├── STARmap_spatial_week_9.png
│   │   │   ├── starmap_umap_celltype.png
│   │   │   ├── starmap_umap_week.png
│   │   │   ├── STARmap_week_8_k15_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_8_k25_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_8_k40_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_8_k8_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_8_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_9_k15_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_9_k25_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_9_k40_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_9_k8_neighbor_z_heatmap.png
│   │   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_clustered.png
│   │   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_fixed.png
│   │   │   ├── STARmap_week_9_neighbor_z_heatmap.png
│   │   │   ├── timecourse_Cytotoxic_NK_enhanced.png
│   │   │   ├── timecourse_Ethanolamine_Metabolism_enhanced.png
│   │   │   ├── timecourse_FLT1_mean.png
│   │   │   ├── timecourse_GALNT1_enhanced.png
│   │   │   ├── timecourse_HLA-G_enhanced.png
│   │   │   ├── timecourse_HLA-G_mean.png
│   │   │   ├── timecourse_Immune_Tolerance_enhanced.png
│   │   │   ├── timecourse_MMP_ECM_Remodeling_enhanced.png
│   │   │   ├── timecourse_MMP2_enhanced.png
│   │   │   ├── timecourse_MMP2_mean.png
│   │   │   ├── timecourse_MMP9_enhanced.png
│   │   │   ├── timecourse_MMP9_mean.png
│   │   │   ├── timecourse_NKG7_enhanced.png
│   │   │   ├── timecourse_PLD1_enhanced.png
│   │   │   ├── timecourse_PLD1_mean.png
│   │   │   ├── timecourse_score_Cytotoxic_NK_mean.png
│   │   │   ├── timecourse_score_Ethanolamine_Metabolism_mean.png
│   │   │   ├── timecourse_score_Immune_Tolerance_mean.png
│   │   │   └── timecourse_score_MMP_ECM_Remodeling_mean.png
│   │   ├── celltype_proportions_by_week_and_version.png
│   │   ├── gene_coordination_dotplot_all.png
│   │   ├── gene_coordination_heatmap.png
│   │   ├── gene_coordination_top_hits_barplot.png
│   │   ├── How to interpret your plot correctly04A.txt
│   │   ├── posthoc_timecourse_HLA.G_all.png
│   │   ├── posthoc_timecourse_HLA.G_significant_only.png
│   │   ├── posthoc_timecourse_MMP2_all.png
│   │   ├── posthoc_timecourse_MMP2_significant_only.png
│   │   ├── posthoc_timecourse_MMP9_all.png
│   │   ├── posthoc_timecourse_NKG7_all.png
│   │   ├── posthoc_timecourse_PLD1_all.png
│   │   ├── posthoc_timecourse_PLD1_significant_only.png
│   │   ├── posthoc_timecourse_score_Cytotoxic_NK_all.png
│   │   ├── posthoc_timecourse_score_Ethanolamine_Metabolism_all.png
│   │   ├── posthoc_timecourse_score_Ethanolamine_Metabolism_significant_only.png
│   │   ├── posthoc_timecourse_score_Immune_Tolerance_all.png
│   │   ├── posthoc_timecourse_score_Immune_Tolerance_significant_only.png
│   │   ├── posthoc_timecourse_score_MMP_ECM_Remodeling_all.png
│   │   ├── posthoc_timecourse_score_MMP_ECM_Remodeling_significant_only.png
│   │   ├── Slide-tags_niche_center_celltype_composition.png
│   │   ├── Slide-tags_niche_clusters_spatial.png
│   │   ├── Slide-tags_niche_geneset_heatmap.png
│   │   ├── Slide-tags_niche_week_composition.png
│   │   ├── SlideTags_adjacency_followup_Erythroblasts_to_FIB1_week_11.png
│   │   ├── SlideTags_adjacency_followup_Erythroblasts_to_FIB1_week_8.png
│   │   ├── SlideTags_adjacency_followup_Erythroblasts_to_FIB1_week_9.png
│   │   ├── SlideTags_adjacency_followup_EVT_progenitor_to_Hofbauer_cells_week_11.png
│   │   ├── SlideTags_adjacency_followup_EVT_progenitor_to_Hofbauer_cells_week_8.png
│   │   ├── SlideTags_adjacency_followup_EVT_progenitor_to_Hofbauer_cells_week_9.png
│   │   ├── SlideTags_adjacency_followup_Hofbauer_cells_to_Endothelial_week_11.png
│   │   ├── SlideTags_adjacency_followup_Hofbauer_cells_to_Endothelial_week_8.png
│   │   ├── SlideTags_adjacency_followup_Hofbauer_cells_to_Endothelial_week_9.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_progenitor_week_11.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_progenitor_week_8.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_progenitor_week_9.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_week_11.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_week_8.png
│   │   ├── SlideTags_adjacency_followup_STB_progenitor_to_STB_week_9.png
│   │   ├── SlideTags_adjacency_followup_STB_to_Endothelial_week_11.png
│   │   ├── SlideTags_adjacency_followup_STB_to_Endothelial_week_8.png
│   │   ├── SlideTags_adjacency_followup_STB_to_Endothelial_week_9.png
│   │   ├── SlideTags_adjacency_followup_STB_to_STB_progenitor_week_11.png
│   │   ├── SlideTags_adjacency_followup_STB_to_STB_progenitor_week_8.png
│   │   ├── SlideTags_adjacency_followup_STB_to_STB_progenitor_week_9.png
│   │   ├── SlideTags_adjacency_followup_vCTB_to_STB_progenitor_week_11.png
│   │   ├── SlideTags_adjacency_followup_vCTB_to_STB_progenitor_week_8.png
│   │   ├── SlideTags_adjacency_followup_vCTB_to_STB_progenitor_week_9.png
│   │   ├── SlideTags_neighbor_enrichment_trend_top.png
│   │   ├── SlideTags_permissiveness_top10pct_composition.png
│   │   ├── SlideTags_permissiveness_week_11.png
│   │   ├── SlideTags_permissiveness_week_8.png
│   │   ├── SlideTags_permissiveness_week_9.png
│   │   ├── SlideTags_spatial_lr_top20.png
│   │   ├── SlideTags_spatial_week_11.png
│   │   ├── SlideTags_spatial_week_8.png
│   │   ├── SlideTags_spatial_week_9.png
│   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_11_k15_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_11_k25_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_11_k40_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_11_k8_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_11_neighbor_z_heatmap.png
│   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_8_k15_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_8_k25_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_8_k40_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_8_k8_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_8_neighbor_z_heatmap.png
│   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_9_k15_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_9_k25_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_9_k40_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_clustered.png
│   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_fixed.png
│   │   ├── SlideTags_week_9_k8_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── SlideTags_week_9_neighbor_z_heatmap.png
│   │   ├── STARmap_niche_center_celltype_composition.png
│   │   ├── STARmap_niche_clusters_spatial.png
│   │   ├── STARmap_niche_geneset_heatmap.png
│   │   ├── STARmap_niche_week_composition.png
│   │   ├── STARmap_permissiveness_top10pct_composition.png
│   │   ├── STARmap_permissiveness_week_8.png
│   │   ├── STARmap_permissiveness_week_9.png
│   │   ├── STARmap_spatial_lr_top20.png
│   │   ├── STARmap_spatial_week_8.png
│   │   ├── STARmap_spatial_week_9.png
│   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_8_k15_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_8_k25_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_8_k40_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_8_k8_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_8_neighbor_z_heatmap.png
│   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_9_k15_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_9_k25_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_9_k40_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_clustered.png
│   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_clustered_bwr.png
│   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_fixed.png
│   │   ├── STARmap_week_9_k8_neighbor_z_heatmap_fixed_bwr.png
│   │   ├── STARmap_week_9_neighbor_z_heatmap.png
│   │   ├── timecourse_Cytotoxic_NK_enhanced.png
│   │   ├── timecourse_Ethanolamine_Metabolism_enhanced.png
│   │   ├── timecourse_GALNT1_enhanced.png
│   │   ├── timecourse_HLA-G_enhanced.png
│   │   ├── timecourse_Immune_Tolerance_enhanced.png
│   │   ├── timecourse_MMP_ECM_Remodeling_enhanced.png
│   │   ├── timecourse_MMP2_enhanced.png
│   │   ├── timecourse_MMP9_enhanced.png
│   │   ├── timecourse_NKG7_enhanced.png
│   │   └── timecourse_PLD1_enhanced.png
│   ├── logs
│   │   ├── 01_preprocess_harmony_embeddings.log
│   │   ├── 02A_preprocess_multiome_reference.log
│   │   ├── 03A_map_slidetags_to_multiome.log
│   │   ├── 03B_map_starmap_to_multiome.log
│   │   ├── 03C_harmonize_celltype_labels.log
│   │   ├── 03D_compare_celltype_proportions_across_versions.log
│   │   ├── 04A_gene_of_interest_timecourse.log
│   │   ├── 04A_gene_of_interest_timecourse_enhanced.log
│   │   ├── 04A_timecourse_posthoc_significance_plots.log
│   │   ├── 04B_immune_subsets_refinement.log
│   │   ├── 04C_gene_coordination_posthoc_plots.log
│   │   ├── 04C_gene_coordination_score.log
│   │   ├── 05A_spatial_overview_plots.log
│   │   ├── 05B_load_and_format_SCP2601.log
│   │   ├── 05B_neighborhood_enrichment.log
│   │   ├── 05C_permissiveness_score_maps.log
│   │   ├── 05D_neighborhood_DE.log
│   │   ├── 05E_spatial_lr_proximity.log
│   │   ├── 05F_interaction_adjacency_followup.log
│   │   ├── 05G_niche_susceptibility_scoring.log
│   │   ├── 05H_spatial_permissiveness_panels_global.log
│   │   ├── 05I_hotspot_vs_coldspot_DE.log
│   │   ├── 06A_cellchat_spatial_constrained.log
│   │   ├── 06C_cellchat_optional.log
│   │   └── 21_rpac_spatial_routes.log
│   ├── objects
│   │   ├── cellchat_slidetags.rds
│   │   ├── multiome_SCP2601_formatted.rds
│   │   ├── SCP2601_combined_spatial.rds
│   │   ├── slidetags_harmonized.rds
│   │   ├── slidetags_mapped.rds
│   │   ├── slidetags_mapped_to_multiome.rds
│   │   ├── slidetags_with_immune_subtypes.rds
│   │   ├── slidetags_with_permissiveness.rds
│   │   ├── starmap_W11_SCP2601.rds
│   │   ├── starmap_W7_SCP2601.rds
│   │   ├── starmap_W8-2_SCP2601.rds
│   │   └── starmap_W9_SCP2601.rds
│   ├── output
│   │   └── tables
│   │       ├── Slide-tags_niche_pseudobulk_de.csv
│   │       ├── Slide-tags_niche_week_center_celltype_composition.csv
│   │       ├── Slide-tags_niche_week_chisq_stdres.csv
│   │       ├── Slide-tags_niche_week_composition.csv
│   │       ├── STARmap_week_8_k15_neighbor_expected.csv
│   │       ├── STARmap_week_8_k15_neighbor_observed.csv
│   │       ├── STARmap_week_8_k15_neighbor_z.csv
│   │       ├── STARmap_week_8_k25_neighbor_expected.csv
│   │       ├── STARmap_week_8_k25_neighbor_observed.csv
│   │       ├── STARmap_week_8_k25_neighbor_z.csv
│   │       ├── STARmap_week_8_k40_neighbor_expected.csv
│   │       ├── STARmap_week_8_k40_neighbor_observed.csv
│   │       ├── STARmap_week_8_k40_neighbor_z.csv
│   │       ├── STARmap_week_8_k8_neighbor_z.csv
│   │       ├── STARmap_week_8_neighbor_expected.csv
│   │       ├── STARmap_week_8_neighbor_observed.csv
│   │       ├── STARmap_week_8_neighbor_z.csv
│   │       ├── STARmap_week_9_k15_neighbor_expected.csv
│   │       ├── STARmap_week_9_k15_neighbor_observed.csv
│   │       ├── STARmap_week_9_k15_neighbor_z.csv
│   │       ├── STARmap_week_9_k25_neighbor_expected.csv
│   │       ├── STARmap_week_9_k25_neighbor_observed.csv
│   │       ├── STARmap_week_9_k25_neighbor_z.csv
│   │       ├── STARmap_week_9_k40_neighbor_expected.csv
│   │       ├── STARmap_week_9_k40_neighbor_observed.csv
│   │       ├── STARmap_week_9_k40_neighbor_z.csv
│   │       ├── STARmap_week_9_k8_neighbor_expected.csv
│   │       ├── STARmap_week_9_k8_neighbor_observed.csv
│   │       ├── STARmap_week_9_k8_neighbor_z.csv
│   │       ├── STARmap_week_9_neighbor_expected.csv
│   │       ├── STARmap_week_9_neighbor_observed.csv
│   │       ├── STARmap_week_9_neighbor_z.csv
│   │       ├── timecourse_gene_module_summaries.csv
│   │       ├── timecourse_gene_module_summaries_enhanced.csv
│   │       ├── timecourse_trend_significance_posthoc.csv
│   │       └── timecourse_trend_significance_posthoc_significant_only.csv
│   ├── reports
│   │   └── 01_methods_and_provenance.md
│   ├── tables
│   │   ├── 05_spatial
│   │   │   ├── DE_results
│   │   │   │   ├── slidetags
│   │   │   │   │   ├── W11_DE_hotspot_vs_background.csv
│   │   │   │   │   ├── W11_volcano.png
│   │   │   │   │   ├── W8-2_DE_hotspot_vs_background.csv
│   │   │   │   │   ├── W8-2_volcano.png
│   │   │   │   │   ├── W9_DE_hotspot_vs_background.csv
│   │   │   │   │   └── W9_volcano.png
│   │   │   │   └── starmap
│   │   │   │       ├── W8-2_DE_hotspot_vs_background.csv
│   │   │   │       ├── W8-2_volcano.png
│   │   │   │       ├── W9_DE_hotspot_vs_background.csv
│   │   │   │       └── W9_volcano.png
│   │   │   ├── slidetags
│   │   │   │   ├── W11
│   │   │   │   │   ├── meta.json
│   │   │   │   │   ├── README_sample.txt
│   │   │   │   │   ├── slidetags_W11_hotspot_cells_top10pct.csv
│   │   │   │   │   ├── slidetags_W11_mapped_coords_global.csv
│   │   │   │   │   └── slidetags_W11_protected_cells_bottom10pct.csv
│   │   │   │   ├── W8-2
│   │   │   │   │   ├── meta.json
│   │   │   │   │   ├── README_sample.txt
│   │   │   │   │   ├── slidetags_W8-2_hotspot_cells_top10pct.csv
│   │   │   │   │   ├── slidetags_W8-2_mapped_coords_global.csv
│   │   │   │   │   └── slidetags_W8-2_protected_cells_bottom10pct.csv
│   │   │   │   └── W9
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       ├── slidetags_W9_hotspot_cells_top10pct.csv
│   │   │   │       ├── slidetags_W9_mapped_coords_global.csv
│   │   │   │       └── slidetags_W9_protected_cells_bottom10pct.csv
│   │   │   ├── slidetags_harmonized
│   │   │   │   ├── W11
│   │   │   │   │   ├── meta.json
│   │   │   │   │   ├── README_sample.txt
│   │   │   │   │   └── slidetags_harmonized_W11_mapped_coords.csv
│   │   │   │   └── W8-2
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       └── slidetags_harmonized_W8-2_mapped_coords.csv
│   │   │   ├── slidetags_with_permissiveness
│   │   │   │   ├── W11
│   │   │   │   │   ├── meta.json
│   │   │   │   │   ├── README_sample.txt
│   │   │   │   │   └── slidetags_with_permissiveness_W11_mapped_coords.csv
│   │   │   │   └── W8-2
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       └── slidetags_with_permissiveness_W8-2_mapped_coords.csv
│   │   │   ├── starmap
│   │   │   │   ├── W8-2
│   │   │   │   │   ├── meta.json
│   │   │   │   │   ├── README_sample.txt
│   │   │   │   │   ├── starmap_W8-2_hotspot_cells_top10pct.csv
│   │   │   │   │   ├── starmap_W8-2_mapped_coords_global.csv
│   │   │   │   │   └── starmap_W8-2_protected_cells_bottom10pct.csv
│   │   │   │   └── W9
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       ├── starmap_W9_hotspot_cells_top10pct.csv
│   │   │   │       ├── starmap_W9_mapped_coords_global.csv
│   │   │   │       └── starmap_W9_protected_cells_bottom10pct.csv
│   │   │   ├── starmap_harmonized
│   │   │   │   └── W8-2
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       └── starmap_harmonized_W8-2_mapped_coords.csv
│   │   │   ├── starmap_with_permissiveness
│   │   │   │   └── W8-2
│   │   │   │       ├── meta.json
│   │   │   │       ├── README_sample.txt
│   │   │   │       └── starmap_with_permissiveness_W8-2_mapped_coords.csv
│   │   │   ├── permissiveness_global_allcells.csv
│   │   │   ├── permissiveness_global_stats.json
│   │   │   ├── PROCESS_LOG.txt
│   │   │   ├── QC_hotspot_protected_counts_per_sample.csv
│   │   │   ├── QC_permissiveness_global_finiteness.csv
│   │   │   ├── SlideTags_hotspot_vs_coldspot_DE.csv
│   │   │   ├── STARmap_hotspot_vs_coldspot_DE.csv
│   │   │   └── summary_table_with_highlights.csv
│   │   ├── 01_qc_summary_by_week.csv
│   │   ├── 06A_cellchat_run_summary.csv
│   │   ├── cellchat_slidetags_interactions.csv
│   │   ├── celltype_proportions_by_week_and_version.csv
│   │   ├── celltype_ratio_to_fibro_by_week_and_version.csv
│   │   ├── celltype_week_trends_by_version.csv
│   │   ├── gene_coordination_scores.csv
│   │   ├── gene_coordination_scores_with_stats.csv
│   │   ├── gene_coordination_significant_only.csv
│   │   ├── Multiome_immune_subtype_counts.csv
│   │   ├── SCP2601_gene_inventory.csv
│   │   ├── SCP2601_spatial_qc_summary.csv
│   │   ├── Slide-tags_niche_center_celltype_composition.csv
│   │   ├── Slide-tags_niche_center_celltype_top5.csv
│   │   ├── Slide-tags_niche_composition.csv
│   │   ├── Slide-tags_niche_fraction_trends.csv
│   │   ├── Slide-tags_niche_geneset_scores.csv
│   │   ├── Slide-tags_niche_markers.csv
│   │   ├── Slide-tags_niche_pseudobulk_de.csv
│   │   ├── Slide-tags_niche_week_center_celltype_composition.csv
│   │   ├── Slide-tags_niche_week_chisq_stdres.csv
│   │   ├── Slide-tags_niche_week_composition.csv
│   │   ├── Slide-tags_niche_week_geneset_scores.csv
│   │   ├── Slide-tags_niche_week_rewiring_jsd.csv
│   │   ├── SlideTags_adjacency_followup_summary.csv
│   │   ├── SlideTags_celltype_harmonization_summary.csv
│   │   ├── SlideTags_immune_subtype_counts.csv
│   │   ├── SlideTags_neighbor_enrichment_effect_summary.csv
│   │   ├── SlideTags_neighbor_enrichment_k_robustness.csv
│   │   ├── SlideTags_neighbor_enrichment_long.csv
│   │   ├── SlideTags_neighbor_enrichment_trend_shortlist.csv
│   │   ├── SlideTags_neighbor_enrichment_trends.csv
│   │   ├── SlideTags_NK_module_coverage_qc.csv
│   │   ├── SlideTags_permissiveness_cell_level.csv
│   │   ├── SlideTags_permissiveness_top10pct_composition.csv
│   │   ├── SlideTags_spatial_lr_edges.csv
│   │   ├── SlideTags_spatial_lr_k_robustness.csv
│   │   ├── SlideTags_spatial_lr_summary.csv
│   │   ├── SlideTags_week_11_k15_neighbor_expected.csv
│   │   ├── SlideTags_week_11_k15_neighbor_observed.csv
│   │   ├── SlideTags_week_11_k15_neighbor_z.csv
│   │   ├── SlideTags_week_11_k25_neighbor_expected.csv
│   │   ├── SlideTags_week_11_k25_neighbor_observed.csv
│   │   ├── SlideTags_week_11_k25_neighbor_z.csv
│   │   ├── SlideTags_week_11_k40_neighbor_expected.csv
│   │   ├── SlideTags_week_11_k40_neighbor_observed.csv
│   │   ├── SlideTags_week_11_k40_neighbor_z.csv
│   │   ├── SlideTags_week_11_k8_neighbor_expected.csv
│   │   ├── SlideTags_week_11_k8_neighbor_observed.csv
│   │   ├── SlideTags_week_11_k8_neighbor_z.csv
│   │   ├── SlideTags_week_11_neighbor_expected.csv
│   │   ├── SlideTags_week_11_neighbor_observed.csv
│   │   ├── SlideTags_week_11_neighbor_z.csv
│   │   ├── SlideTags_week_8_k15_neighbor_expected.csv
│   │   ├── SlideTags_week_8_k15_neighbor_observed.csv
│   │   ├── SlideTags_week_8_k15_neighbor_z.csv
│   │   ├── SlideTags_week_8_k25_neighbor_expected.csv
│   │   ├── SlideTags_week_8_k25_neighbor_observed.csv
│   │   ├── SlideTags_week_8_k25_neighbor_z.csv
│   │   ├── SlideTags_week_8_k40_neighbor_expected.csv
│   │   ├── SlideTags_week_8_k40_neighbor_observed.csv
│   │   ├── SlideTags_week_8_k40_neighbor_z.csv
│   │   ├── SlideTags_week_8_k8_neighbor_expected.csv
│   │   ├── SlideTags_week_8_k8_neighbor_observed.csv
│   │   ├── SlideTags_week_8_k8_neighbor_z.csv
│   │   ├── SlideTags_week_8_neighbor_expected.csv
│   │   ├── SlideTags_week_8_neighbor_observed.csv
│   │   ├── SlideTags_week_8_neighbor_z.csv
│   │   ├── SlideTags_week_9_k15_neighbor_expected.csv
│   │   ├── SlideTags_week_9_k15_neighbor_observed.csv
│   │   ├── SlideTags_week_9_k15_neighbor_z.csv
│   │   ├── SlideTags_week_9_k25_neighbor_expected.csv
│   │   ├── SlideTags_week_9_k25_neighbor_observed.csv
│   │   ├── SlideTags_week_9_k25_neighbor_z.csv
│   │   ├── SlideTags_week_9_k40_neighbor_expected.csv
│   │   ├── SlideTags_week_9_k40_neighbor_observed.csv
│   │   ├── SlideTags_week_9_k40_neighbor_z.csv
│   │   ├── SlideTags_week_9_k8_neighbor_expected.csv
│   │   ├── SlideTags_week_9_k8_neighbor_observed.csv
│   │   ├── SlideTags_week_9_k8_neighbor_z.csv
│   │   ├── SlideTags_week_9_neighbor_expected.csv
│   │   ├── SlideTags_week_9_neighbor_observed.csv
│   │   ├── SlideTags_week_9_neighbor_z.csv
│   │   ├── STARmap_celltype_harmonization_summary.csv
│   │   ├── STARmap_immune_subtype_counts.csv
│   │   ├── STARmap_neighbor_enrichment_effect_summary.csv
│   │   ├── STARmap_neighbor_enrichment_k_robustness.csv
│   │   ├── STARmap_neighbor_enrichment_long.csv
│   │   ├── STARmap_neighbor_enrichment_trend_shortlist.csv
│   │   ├── STARmap_neighbor_enrichment_trends.csv
│   │   ├── STARmap_niche_center_celltype_composition.csv
│   │   ├── STARmap_niche_center_celltype_top5.csv
│   │   ├── STARmap_niche_composition.csv
│   │   ├── STARmap_niche_fraction_trends.csv
│   │   ├── STARmap_niche_geneset_scores.csv
│   │   ├── STARmap_niche_markers.csv
│   │   ├── STARmap_niche_pseudobulk_de.csv
│   │   ├── STARmap_niche_susceptibility_classification.csv
│   │   ├── STARmap_niche_susceptibility_scores.csv
│   │   ├── STARmap_niche_susceptibility_week_composition.csv
│   │   ├── STARmap_niche_susceptibility_week_trends.csv
│   │   ├── STARmap_niche_week_center_celltype_composition.csv
│   │   ├── STARmap_niche_week_chisq_stdres.csv
│   │   ├── STARmap_niche_week_composition.csv
│   │   ├── STARmap_niche_week_geneset_scores.csv
│   │   ├── STARmap_niche_week_rewiring_jsd.csv
│   │   ├── STARmap_NK_module_coverage_qc.csv
│   │   ├── STARmap_permissiveness_cell_level.csv
│   │   ├── STARmap_permissiveness_top10pct_composition.csv
│   │   ├── STARmap_spatial_lr_edges.csv
│   │   ├── STARmap_spatial_lr_k_robustness.csv
│   │   ├── STARmap_spatial_lr_summary.csv
│   │   ├── STARmap_week_8_k15_neighbor_expected.csv
│   │   ├── STARmap_week_8_k15_neighbor_observed.csv
│   │   ├── STARmap_week_8_k15_neighbor_z.csv
│   │   ├── STARmap_week_8_k25_neighbor_expected.csv
│   │   ├── STARmap_week_8_k25_neighbor_observed.csv
│   │   ├── STARmap_week_8_k25_neighbor_z.csv
│   │   ├── STARmap_week_8_k40_neighbor_expected.csv
│   │   ├── STARmap_week_8_k40_neighbor_observed.csv
│   │   ├── STARmap_week_8_k40_neighbor_z.csv
│   │   ├── STARmap_week_8_k8_neighbor_expected.csv
│   │   ├── STARmap_week_8_k8_neighbor_observed.csv
│   │   ├── STARmap_week_8_k8_neighbor_z.csv
│   │   ├── STARmap_week_8_neighbor_expected.csv
│   │   ├── STARmap_week_8_neighbor_observed.csv
│   │   ├── STARmap_week_8_neighbor_z.csv
│   │   ├── STARmap_week_9_k15_neighbor_expected.csv
│   │   ├── STARmap_week_9_k15_neighbor_observed.csv
│   │   ├── STARmap_week_9_k15_neighbor_z.csv
│   │   ├── STARmap_week_9_k25_neighbor_expected.csv
│   │   ├── STARmap_week_9_k25_neighbor_observed.csv
│   │   ├── STARmap_week_9_k25_neighbor_z.csv
│   │   ├── STARmap_week_9_k40_neighbor_expected.csv
│   │   ├── STARmap_week_9_k40_neighbor_observed.csv
│   │   ├── STARmap_week_9_k40_neighbor_z.csv
│   │   ├── STARmap_week_9_k8_neighbor_expected.csv
│   │   ├── STARmap_week_9_k8_neighbor_observed.csv
│   │   ├── STARmap_week_9_k8_neighbor_z.csv
│   │   ├── STARmap_week_9_neighbor_expected.csv
│   │   ├── STARmap_week_9_neighbor_observed.csv
│   │   ├── STARmap_week_9_neighbor_z.csv
│   │   ├── timecourse_gene_module_summaries.csv
│   │   ├── timecourse_gene_module_summaries_enhanced.csv
│   │   ├── timecourse_trend_significance_posthoc.csv
│   │   ├── timecourse_trend_significance_posthoc_significant_only.csv
│   │   └── vulnerability_output_manifest.csv
│   └── troubleshooting
│       └── 01_preprocess_harmony_embeddings
│           ├── week_W7
│           │   ├── diagnostics.txt
│           │   ├── expression_preview_head200x50.csv
│           │   ├── expression_raw_first20lines.txt
│           │   ├── Made-Manual-expression_raw_first20lines.txt
│           │   └── spots_preview_head500.csv
│           └── week_W8-2
│               ├── diagnostics.txt
│               └── spots_preview_head500.csv
├── Papers
│   ├── Acute response to pathogens in the early human development main.pdf
│   ├── Bioengineering   Transla Med - 2022 - Cui - Engineering placenta‐like organoids containing endogenous vascular cells from.pdf
│   ├── RPAC-Deep_Pathway_Analysis_V2.0_A_Pathway_Analysis_Framework_Incorporating_Multi-Dimensional_Omics_Data.pdf
│   ├── RPAC-pj_phd_defense_slides_20220421_1000.pdf
│   └── Single-cell reconstruction of the early maternal–fetal interface in.pdf
├── Projects-to-consider-for-alt-aproach
│   ├── BO'sWork
│   │   ├── Bo Thesis final.pdf
│   │   └── Metagene paper text.pdf
│   └── RPAC
│       ├── pj_phd_defense_slides_20220421_1000.pdf
│       └── Rpac.pdf
├── scripts
│   ├── 00_archive
│   │   ├── 02_preprocess
│   │   │   └── 02A_preprocess_multiome_reference.R
│   │   ├── 03_mapping
│   │   │   ├── 03A_map_slidetags_to_multiome.R
│   │   │   ├── 03B_map_starmap_to_multiome.R
│   │   │   ├── 03B_map_starmap_to_multiome_original.R
│   │   │   ├── 03C_harmonize_celltype_labels.R
│   │   │   └── 03D_compare_celltype_proportions_across_versions.R
│   │   ├── 04_timecourse
│   │   │   ├── 04A_gene_of_interest_timecourse-2.R
│   │   │   ├── 04A_gene_of_interest_timecourse.R
│   │   │   ├── 04A_gene_of_interest_timecourse_ENHANCED.R
│   │   │   ├── 04A_timecourse_posthoc_significance_plots.R.R
│   │   │   ├── 04B_immune_subsets_refinement.R
│   │   │   ├── 04C_gene_coordination_score.R
│   │   │   └── 4C_gene_coordination_posthoc_plots.R.R
│   │   ├── 05_spatial
│   │   │   ├── 05A_spatial_overview_plots.R
│   │   │   ├── 05B_neighborhood_enrichment.R
│   │   │   ├── 05C_permissiveness_score_maps.R
│   │   │   ├── 05D_neighborhood_DE.R
│   │   │   ├── 05E_spatial_lr_proximity.R
│   │   │   ├── 05F_and_05H_description.md
│   │   │   ├── 05F_description
│   │   │   ├── 05F_description.txt
│   │   │   ├── 05F_interaction_adjacency_followup-2.R
│   │   │   ├── 05F_interaction_adjacency_followup.R
│   │   │   ├── 05H_description
│   │   │   ├── 05H_description.txt
│   │   │   ├── 05H_qc_summary.R
│   │   │   ├── 05H_spatial_permissiveness_panels_global-original.R
│   │   │   ├── 05H_spatial_permissiveness_panels_global.R
│   │   │   ├── 05H_spatial_permissiveness_panels_global_W_FigurLegends.R
│   │   │   ├── 05I_DE_hotspot_vs_background.R
│   │   │   ├── 5g.R
│   │   │   └── UntitledR.R
│   │   ├── 06_cell_communication
│   │   │   ├── 06A_cellchat_spatial_constrained.R
│   │   │   ├── 06B_simple_LR_scoring.R
│   │   │   └── 06C_cellchat_optional.R
│   │   ├── 07_export
│   │   │   └── 07A_export_shareable_outputs.R
│   │   ├── 08_metagene
│   │   │   ├── 08A_housekeeping_qc.R
│   │   │   └── 08B_metagene_nmf.R
│   │   ├── 08_metagenes
│   │   │   ├── 08A_housekeeping_diagnostics.R
│   │   │   ├── 08B_metagene_module_discovery.R
│   │   │   └── 08C_metagene_spatiotemporal_maps.R
│   │   ├── AltAnalysisScripts
│   │   │   ├── 05B_load_and_format_SCP2601.R
│   │   │   ├── 16_compute_misi_SCP2601_spatial.R
│   │   │   └── 21_rpac_spatial_routes.R
│   │   ├── 06C_cellchat_optional.R
│   │   ├── 16_compute_misi.R
│   │   ├── 16b_misi_v2_subindices.R
│   │   ├── 16c_plot_misi_embeddings_V2.R
│   │   ├── 16d_plot_misi_violins.R
│   │   ├── 17_cellchat_by_condition.R
│   │   ├── 18_cellphonedb_export.R
│   │   ├── 21_rpac_v2_corrected_routes.R
│   │   └── utils.r
│   ├── 01_active_pipeline
│   │   └── 01_preprocess_harmony_embeddings.R
│   └── R
│       ├── utils-old.R
│       ├── utils-old2.R
│       ├── utils-old3.R
│       └── utils.R
├── Archiver_Script.R
├── Create-DorectoryTree.r
├── CRITICAL_FIXES.md
├── data_inventory_report.txt
├── directory_tree.md
├── directory_tree.txt
├── hPlacenta-architecture.Rproj
├── MASTER_PROJECT_MANIFEST.txt
├── Protect-GitHub.ps1
├── README-Project.md
├── README-updated-pipline.md
├── README.md
├── README_aRCHIVE.md
├── RUN_PIPELINE_ADVANCED.R
└── RUN_PIPELINE_MINIMAL.R
