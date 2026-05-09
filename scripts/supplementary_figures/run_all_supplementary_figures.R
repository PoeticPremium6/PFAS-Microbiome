# Run all supplementary figure scripts from repository root
scripts <- c(
  "scripts/supplementary_figures/00_supplementary_pcoa_elimination_rate.R",
  "scripts/supplementary_figures/01_suppfig_s1_species_deseq2_log2fc.R",
  "scripts/supplementary_figures/02_suppfig_s2_family_and_elimination_boxplots.R",
  "scripts/supplementary_figures/03_suppfig_s3_compound_specific_taxa_correlations.R",
  "scripts/supplementary_figures/04_suppfig_s4_species_relative_abundance_dotplot.R",
  "scripts/supplementary_figures/05_suppfig_s5_s6_rf_qc_stability.R",
  "scripts/supplementary_figures/06_suppfig_s7_subsystem_cluster_diagnostics.R"
)

for (script in scripts) {
  message("\nRunning: ", script)
  source(script)
}
