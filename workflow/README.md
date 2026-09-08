# Canonical public workflow

This is the shortest retained public analysis lineage. It is not a development history.

`00_revision_plot_style.R` is the shared plotting helper. The remaining scripts are grouped below in execution/conceptual order.

## Taxonomic profiles

- `analysis/07F_PHASE4_merge_metaphlan_profiles.py`

## Public microbiome context

- `analysis/08B_PHASE5_cMD_metadata_and_raw_accession_inventory.R`
- `analysis/08F_PHASE5_ExperimentHub_cMD_resource_inventory.R`
- `analysis/08G_PHASE5_extract_EH_cMD_SE_tables.R`
- `analysis/08H_PHASE5_extract_matrix_only_public_tables.R`
- `analysis/08I_PHASE5_public_preharmonization_and_waiting_ronneby.R`

## Functional analysis

- `analysis/09A_PHASE6_build_humann_manifest.py`
- `analysis/10_PHASE6_merge_humann_outputs.sh`
- `analysis/11_PHASE5_harmonize_public_HUMAnN_pathways.py`
- `analysis/13_PHASE6_prepare_HUMAnN_analysis_matrices.py`
- `analysis/14_PHASE6_HUMAnN_pathway_ordination.R`
- `analysis/16_PHASE6_HUMAnN_association_tests.R`

## Genome-resolved analysis

- `analysis/19_PHASE7A_MAG_quality_summary.sh`
- `analysis/24_PHASE7A_MAG_CoverM_abundance_ordination.sh`
- `analysis/96A_build_public_MAG_catalogue.py`
- `analysis/98A_build_MAG_sample_read_provenance.py`

## PFAS context

- `analysis/85A_rebuild_item3_toxicokinetic_summary.R`

## Publication figures

- `final_figures/build_figure_1_excellence.R` through `build_figure_6_excellence.R`: final main figures.
- `final_figures/supplementary/build_figure_S01.R` through `build_figure_S06.R`: final supplementary figures.

Raw/private preprocessing and internal HPC orchestration are excluded because their underlying restricted inputs are not distributed here.
