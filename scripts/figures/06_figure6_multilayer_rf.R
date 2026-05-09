# FIGURE 6 MASTER SCRIPT — COMPLETE UPDATED VERSION
#
# Includes:
#   - Taxa features from vst_physeq
#   - Reaction features from ReactionAbundance.csv
#   - Subsystem features from SubsystemAbundance.csv
#   - Metabolite features from Metabolite_Clean_Summarized_Test.csv
#
# Main analyses:
#   - Spearman correlations with average PFAS elimination rate
#   - Feature network
#   - Repeated random forest prediction and stability analysis
#   - RF stable-feature correlation heatmap
#   - Supplementary RF QC and stability plots
#
# Important fixes:
#   - Handles mismatched sample columns across feature matrices.
#   - Robustly builds heatmap table without relying on Var1/Var2/Freq.
#   - Adds metabolite-label map for VMH-like metabolite IDs.
#   - Exports untranslated metabolite labels for manual curation.
#
# Requirements:
#   - vst_physeq already loaded in the R environment
############################################################
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(scales)
  library(patchwork)
  library(readr)
  library(stringr)
  library(forcats)
  library(caret)
  library(ranger)
})

# ============================================================
# 1. PATHS
# ============================================================

reaction_path <- file.path("data", "model_outputs", "ReactionAbundance.csv")
subsystem_path <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
metabolic_data_path <- file.path("data", "model_outputs", "Metabolite_Clean_Summarized_Test.csv")

out_dir <- file.path("results", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_fig6_png <- file.path(out_dir, "Figure6_Final_A_to_D_WithMetabolites.png")
out_fig6_pdf <- file.path(out_dir, "Figure6_Final_A_to_D_WithMetabolites.pdf")

out_csv_all <- file.path(out_dir, "Figure6_PFAS_Elimination_Correlations_AllFeatures_WithMetabolites.csv")
out_csv_sig <- file.path(out_dir, "Figure6_PFAS_Elimination_Correlations_FilteredForPlot_WithMetabolites.csv")
out_csv_nodes <- file.path(out_dir, "Figure6_NetworkNodes_TopNPerType_WithMetabolites.csv")
out_csv_edges <- file.path(out_dir, "Figure6_NetworkEdges_TopNPerType_WithMetabolites.csv")

csv_perf <- file.path(out_dir, "Supp_RF_Performance_ByRun_WithMetabolites.csv")
csv_resid_stats <- file.path(out_dir, "Supp_RF_ResidualStats_ByRun_WithMetabolites.csv")
csv_importance_all <- file.path(out_dir, "Supp_RF_Importance_AllRuns_WithMetabolites.csv")
csv_selected_topk <- file.path(out_dir, "Supp_RF_Selected_Top25_ByRun_WithMetabolites.csv")
csv_stability <- file.path(out_dir, "Supp_RF_Stability_Summary_Top25Freq_WithMetabolites.csv")
csv_topstable_type <- file.path(out_dir, "Supp_RF_TopStableFeatures_PerType_Top25_WithMetabolites.csv")
csv_topdisplay <- file.path(out_dir, "Supp_RF_TopFeatures_DisplayedInStabilityFig_Top25_WithMetabolites.csv")
csv_stablecand <- file.path(out_dir, "Supp_RF_StableCandidates_Top25Freq_ge0.25_WithMetabolites.csv")
csv_stab_overall <- file.path(out_dir, "Supp_RF_Stability_OverallSummary_WithMetabolites.csv")
csv_stab_bytype <- file.path(out_dir, "Supp_RF_Stability_ByType_WithMetabolites.csv")
csv_rf_filtering <- file.path(out_dir, "Supp_RF_FeatureFiltering_Diagnostics_WithMetabolites.csv")
csv_rf_feature_heatmap <- file.path(out_dir, "Supp_RF_TopStableFeature_CorrelationHeatmapData_WithMetabolites.csv")

csv_metabolite_untranslated <- file.path(out_dir, "Figure6_MetaboliteLabels_StillNeedTranslation.csv")
txt_metabolite_template <- file.path(out_dir, "Figure6_MetaboliteLabels_ManualMap_Template.txt")

fig_qc_4panel <- file.path(out_dir, "SuppFig_RF_QC_100Runs_4Panel_WithMetabolites.png")
fig_qc_4panel_pdf <- file.path(out_dir, "SuppFig_RF_QC_100Runs_4Panel_WithMetabolites.pdf")

fig_stability_4panel <- file.path(out_dir, "SuppFig_RF_Stability_100Runs_4Panel_WithMetabolites.png")
fig_stability_4panel_pdf <- file.path(out_dir, "SuppFig_RF_Stability_100Runs_4Panel_WithMetabolites.pdf")

fig_rf_extra <- file.path(out_dir, "SuppFig_RF_TopStable_Features_Extra_WithMetabolites.png")
fig_rf_extra_pdf <- file.path(out_dir, "SuppFig_RF_TopStable_Features_Extra_WithMetabolites.pdf")

# ============================================================
# 2. PARAMETERS
# ============================================================

rho_threshold <- 0.12
pval_threshold <- 0.20
max_nodes_per_type <- 5
top_summary_per_type <- 10

min_shared_samples <- 4
min_pairs_cor <- min_shared_samples

seed_base <- 42
train_fraction <- 0.8
n_runs <- 100
top_k <- 25

num_trees <- 1000
min_node_size <- 5
cv_folds <- 5
cv_repeats <- 1

top_rf_features_per_type_main <- 5
top_rf_features_heatmap <- 18

debug_fast_rf <- FALSE
if (debug_fast_rf) {
  n_runs <- 5
  num_trees <- 250
}

show_plots_in_rstudio <- FALSE

feature_type_levels <- c("Taxa", "Reaction", "Subsystem", "Metabolite")
feature_type_levels_with_unknown <- c(feature_type_levels, "Unknown")

# ============================================================
# 3. MANUAL LABEL MAPS
# ============================================================

reaction_label_map <- c(
  "3OAACPR2" = "3-oxoacyl-ACP reductase",
  "3HACPR1" = "3-hydroxyacyl-ACP reductase",
  "4HTHRS" = "4-hydroxy-L-threonine synthase",
  "AAGAAH" = "N-acetylglucosamine aminohydrolase",
  "AMY2e" = "Extracellular amylase",
  "BACCL" = "Biotin carboxyl carrier protein ligase",
  "CEPA" = "Cephalosporin biosynthesis reaction",
  "EX core3(e)" = "Exchange of O-glycan core 3",
  "EX core4(e)" = "Exchange of O-glycan core 4",
  "EX core6(e)" = "Exchange of O-glycan core 6",
  "EX core7(e)" = "Exchange of O-glycan core 7",
  "EX glycogen2(e)" = "Exchange of glycogen structure 2",
  "GACPCD" = "Glutaryl-ACP decarboxylase",
  "GALASE OGLYCAN2epp" = "Periplasmic O-glycan galactosidase 2",
  "GALASE OGLYCAN3epp" = "Periplasmic O-glycan galactosidase 3",
  "GAMYe" = "Extracellular glucoamylase",
  "GLCNACASE OGLYCAN1epp" = "Periplasmic O-glycan N-acetylglucosaminidase 1",
  "GLYt4r" = "Reversible glycine transport",
  "Glycine transport" = "Glycine transport",
  "L-aspartate:nad[c]P+ oxidoreductase (deaminating)" = "L-aspartate NAD(P)+ oxidoreductase",
  "MALCOAMT" = "Malonyl-CoA methyltransferase",
  "RBK Dr" = "Ribokinase",
  "r1116" = "Uncharacterised reaction r1116",
  "SALCS2" = "Salicylate CoA synthase",
  "SIAASE OGLYCAN2epp" = "Periplasmic O-glycan sialidase 2",
  "sink s" = "Sink reaction",
  "SQLS" = "Squalene synthase",
  "S-Propane-1,2-diol facilitated transport" = "S-propane-1,2-diol facilitated transport",
  "T ANTIGENtex" = "T-antigen transport",
  "DM 4d 5" = "Demand reaction DM 4d 5",
  "DM 2HYM2EPH" = "Demand reaction DM 2HYM2EPH",
  "1,4-dihydroxy-2-naphthoate octaprenyltransferase" =
    "1,4-dihydroxy-2-naphthoate octaprenyltransferase",
  "1,4-dihydroxy-2-naphthoate Periplasm Transport" =
    "1,4-dihydroxy-2-naphthoate periplasmic transport",
  "Exchange of 1,4-Dihydroxy-2-naphthoate" =
    "Exchange of 1,4-dihydroxy-2-naphthoate",
  "10DMMCNFDOR" = "Dimethylmenaquinone:ferricytochrome oxidoreductase",
  "12DGR180t" = "1,2-diacylglycerol transport",
  "12PPD Stex" = "S-propane-1,2-diol extracellular transport",
  "12PPDtpp" = "S-propane-1,2-diol periplasmic transport",
  "13PPDH" = "Propanediol dehydrogenase",
  "13PPDtex" = "R-propane-1,3-diol extracellular transport",
  "13PPDtpp" = "R-propane-1,3-diol periplasmic transport",
  "15DAPtex" = "Diaminopentane extracellular transport",
  "1H2NPTH" = "1-hydroxy-2-naphthoate hydrolase",
  "1P4H2CBXLAH" = "1-pyrroline-4-hydroxy-2-carboxylate dehydrogenase",
  "22IDPOR" = "2,2-iminodipropanoate oxidoreductase",
  "23DAPAL" = "2,3-diaminopropionate ammonia-lyase",
  "24DCOAR" = "2,4-dienoyl-CoA reductase",
  "26DAPt2r" = "2,6-diaminopimelate reversible transport",
  "2DDPENTHL" = "2-dehydro-3-deoxy-D-pentonate hydrolase",
  "2DHPL" = "2-dehydropantoate aldolase",
  "2FBZOR" = "2-fluorobenzoate oxidoreductase",
  "2FBZTOR" = "2-fluorobenzoyl-CoA reductase",
  "2IMZS" = "2-isopropylmalate synthase",
  "2INSD" = "Inosine deaminase",
  "2IPDPIPT" = "Dimethylallyltranstransferase",
  "2OBUTFDXOR" = "2-oxobutyrate:ferredoxin oxidoreductase",
  "2OXOADOX" = "2-oxoadipate dehydrogenase",
  "33HPHACt2r" = "3-(3-hydroxyphenyl)propionate reversible transport",
  "34DHOXPEGOX" = "3,4-dihydroxyphenylglycol oxidase",
  "34DHPHAt2r" = "3,4-dihydroxyphenylacetate reversible transport",
  "35CGMPt2" = "3,5-cyclic GMP transport",
  "36DAHXI" = "3,6-dideoxy-D-arabino-hexose isomerase",
  "3DGUR" = "3-dehydroglucuronate reductase",
  "3FBZOR" = "3-fluorobenzoate oxidoreductase",
  "3FBZTOR" = "3-fluorobenzoyl-CoA reductase",
  "3FCHLOR" = "3-fluorocatechol chlorohydrolase",
  "3HACPR2" = "3-hydroxyacyl-ACP reductase 2",
  "3HAD100" = "3-hydroxyacyl-CoA dehydrogenase C10:0",
  "3HAD120" = "3-hydroxyacyl-CoA dehydrogenase C12:0",
  "3HAD140" = "3-hydroxyacyl-CoA dehydrogenase C14:0",
  "3HAD160" = "3-hydroxyacyl-CoA dehydrogenase C16:0",
  "3HAD161" = "3-hydroxyacyl-CoA dehydrogenase C16:1",
  "3HAD180" = "3-hydroxyacyl-CoA dehydrogenase C18:0",
  "3HAD40" = "3-hydroxyacyl-CoA dehydrogenase C4:0",
  "3HAD60" = "3-hydroxyacyl-CoA dehydrogenase C6:0",
  "3HAD80" = "3-hydroxyacyl-CoA dehydrogenase C8:0",
  "3MEACMPte" = "3-methyladenine cAMP transport",
  "3MOBDC2" = "3-methyl-2-oxobutanoate decarboxylase",
  "3MOPDC2" = "3-methyl-2-oxopentanoate decarboxylase",
  "3MOPLPAMO" = "3-methyl-2-oxopentanoate lipoamide oxidoreductase",
  "3OAACPR1" = "3-oxoacyl-ACP reductase 1",
  "3OADPCOAT" = "3-oxoadipyl-CoA thiolase",
  "3OPCPOOR" = "3-oxopropionyl-CoA oxidoreductase",
  "4ABUTtex" = "4-aminobutyrate extracellular transport",
  "4AHMMPtex" = "4-amino-5-hydroxymethyl-2-methylpyrimidine extracellular transport",
  "4FBZOR" = "4-fluorobenzoate oxidoreductase",
  "4FCHLOR" = "4-fluorocatechol chlorohydrolase",
  "4GBTAH" = "4-guanidinobutanoate amidinohydrolase",
  "4HBHYOXy" = "4-hydroxybenzoate hydroxylase",
  "4HBZCL" = "4-hydroxybenzoate CoA ligase"
)

# Fill or correct these manually as needed after checking:
# Figure6_MetaboliteLabels_StillNeedTranslation.csv
metabolite_label_map <- c(
  "galt[fe]" = "Galactitol exchange",
  "so4[fe]" = "Sulfate exchange",
  "actn R[u]" = "(R)-Acetoin",
  "cmp[fe]" = "CMP exchange",
  "glcn[fe]" = "Gluconate exchange",
  "h[fe]" = "Proton exchange",
  "co2528[ve]" = "Unknown Metabolite (co2528)",
  "glygly2[fe]" = "Glycylglycine exchange",
  "xylottr[u]" = "Xylotriose transport",
  "zn2[fe]" = "Zinc ion exchange",
  "dcholf[fe]" = "Deoxycholate exchange",
  "gbbtn[fe]" = "Gamma-butyrobetaine exchange",
  "mqn8[fe]" = "Menaquinone-8 exchange",
  "ca2[d]" = "Calcium ion exchange",
  "rblfrd[fe]" = "Riboflavin-derived metabolite exchange",
  "Tn antigen[fe]" = "Tn antigen exchange"
)

# ============================================================
# 4. HELPERS
# ============================================================

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "nd", "ND", "na", "n/a", "N/A", "<LOD", "<LOQ", "< LOQ")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

clean_label <- function(x, max_chars = Inf) {
  x <- as.character(x)
  x <- str_replace_all(x, "^Taxa\\.Taxa\\:\\:", "")
  x <- str_replace_all(x, "^Reaction\\.Reaction\\:\\:", "")
  x <- str_replace_all(x, "^Subsystem\\.Subsystem\\:\\:", "")
  x <- str_replace_all(x, "^Metabolite\\.Metabolite\\:\\:", "")
  x <- str_replace_all(x, "^Taxa\\:\\:", "")
  x <- str_replace_all(x, "^Reaction\\:\\:", "")
  x <- str_replace_all(x, "^Subsystem\\:\\:", "")
  x <- str_replace_all(x, "^Metabolite\\:\\:", "")
  x <- str_replace_all(x, "^(Taxa|OTU|Reaction|Subsystem|Metabolite)[\\.\\:\\_\\-\\s]+", "")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  
  if (is.finite(max_chars)) {
    x <- str_trunc(x, width = max_chars)
  }
  
  x
}

standardise_reaction_key <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  str_trim(x)
}

humanise_feature_label <- function(label, type = NULL, max_chars = Inf) {
  label_clean <- clean_label(label, max_chars = Inf)
  
  if (!is.null(type)) {
    type_chr <- as.character(type)
  } else {
    type_chr <- rep(NA_character_, length(label_clean))
  }
  
  out <- label_clean
  
  reaction_key <- standardise_reaction_key(label_clean)
  reaction_mapped <- reaction_label_map[reaction_key]
  
  reaction_idx <- is.na(type_chr) | type_chr == "Reaction" | str_detect(label_clean, "^Reaction")
  out[reaction_idx & !is.na(reaction_mapped)] <- unname(
    reaction_mapped[reaction_idx & !is.na(reaction_mapped)]
  )
  
  metabolite_idx <- !is.na(type_chr) & type_chr == "Metabolite"
  metabolite_key <- label_clean
  metabolite_mapped <- metabolite_label_map[metabolite_key]
  out[metabolite_idx & !is.na(metabolite_mapped)] <- unname(
    metabolite_mapped[metabolite_idx & !is.na(metabolite_mapped)]
  )
  
  out <- str_replace_all(out, "^EX ([A-Za-z0-9]+)\\(e\\)$", "Exchange of \\1")
  out <- str_replace_all(out, "^DM ([A-Za-z0-9 ]+)$", "Demand reaction \\1")
  out <- str_replace_all(out, "^sink ([A-Za-z0-9 ]+)$", "Sink reaction \\1")
  
  clean_label(out, max_chars = max_chars)
}

make_unique_safe <- function(x) {
  make.unique(as.character(x), sep = "_dup")
}

theme_clean <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      panel.grid.minor = element_blank()
    )
}

safe_print <- function(p) {
  if (isTRUE(show_plots_in_rstudio)) {
    tryCatch(
      print(p),
      error = function(e) {
        message("Plot was saved, but RStudio plot pane could not display it: ", e$message)
      }
    )
  }
}

get_taxa_labels <- function(physeq_obj, otu_rownames) {
  if (!is.null(tax_table(physeq_obj, errorIfNULL = FALSE))) {
    tax_df <- as.data.frame(tax_table(physeq_obj))
    tax_df$TaxonID <- rownames(tax_df)
    
    if ("Species" %in% colnames(tax_df)) {
      labels <- as.character(tax_df$Species)
      bad <- is.na(labels) | labels == "" | labels == "NA"
      labels[bad] <- tax_df$TaxonID[bad]
    } else {
      labels <- tax_df$TaxonID
    }
    
    names(labels) <- tax_df$TaxonID
    labels <- labels[otu_rownames]
    labels[is.na(labels)] <- otu_rownames[is.na(labels)]
    return(labels)
  }
  
  otu_rownames
}

is_vmh_like_label <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  
  short_code <- nchar(x) <= 16 &
    str_detect(x, "^[A-Za-z0-9_\\-\\(\\)\\[\\]\\+]+$") &
    !str_detect(x, "\\s")
  
  vmh_exchange <- str_detect(x, "^EX[_\\s]")
  vmh_demand <- str_detect(x, "^DM[_\\s]")
  vmh_sink <- str_detect(x, "^sink[_\\s]")
  r_code <- str_detect(x, "^r[0-9]+$")
  epp_or_tex <- str_detect(x, "epp$|tex$")
  
  short_code | vmh_exchange | vmh_demand | vmh_sink | r_code | epp_or_tex
}

is_vmh_metabolite_like_label <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  
  str_detect(x, "^[A-Za-z0-9_]+\\[[a-zA-Z]+\\]$") |
    str_detect(x, "^[A-Za-z0-9_]+\\([a-zA-Z]+\\)$") |
    (
      nchar(x) <= 14 &
        str_detect(x, "^[A-Za-z0-9_]+$") &
        !str_detect(x, "\\s")
    )
}

get_plotted_reaction_labels_needing_translation <- function(plotted_df, out_dir = NULL) {
  if (is.null(plotted_df) || nrow(plotted_df) == 0) {
    return(tibble())
  }
  
  out <- plotted_df %>%
    mutate(
      type = as.character(type),
      original_label = as.character(feature_label),
      human_label = humanise_feature_label(original_label, type = type, max_chars = Inf),
      still_cryptic = type == "Reaction" & is_vmh_like_label(human_label)
    ) %>%
    filter(type == "Reaction", still_cryptic) %>%
    distinct(source, type, model_feature, feature_key, original_label, human_label) %>%
    arrange(source, original_label)
  
  if (!is.null(out_dir)) {
    write_csv(
      out,
      file.path(out_dir, "Figure6_PLOTTED_ReactionLabels_StillNeedTranslation_WithMetabolites.csv")
    )
    
    if (nrow(out) > 0) {
      writeLines(
        paste0('"', out$original_label, '"', collapse = ",\n"),
        file.path(out_dir, "Figure6_PLOTTED_ReactionLabels_StillNeedTranslation_CopyPaste_WithMetabolites.txt")
      )
    }
  }
  
  out
}

compute_feature_cor <- function(feature_matrix, rates_vec, min_pairs = 4) {
  features <- rownames(feature_matrix)
  
  res_list <- lapply(features, function(f) {
    x <- to_numeric_safe(feature_matrix[f, ])
    y <- rates_vec[colnames(feature_matrix)]
    
    good <- !is.na(x) & !is.na(y)
    n_pairs <- sum(good)
    
    if (n_pairs < min_pairs) return(NULL)
    if (var(x[good], na.rm = TRUE) == 0) return(NULL)
    
    ct <- suppressWarnings(
      tryCatch(
        cor.test(x[good], y[good], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
    )
    
    if (is.null(ct)) return(NULL)
    
    data.frame(
      feature_key = f,
      rho = as.numeric(ct$estimate),
      p = ct$p.value,
      n_pairs = n_pairs,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(res_list)
}

guess_id_col <- function(dt, preferred_names = NULL) {
  if (!is.null(preferred_names)) {
    hit <- preferred_names[preferred_names %in% colnames(dt)]
    if (length(hit) > 0) return(hit[1])
  }
  colnames(dt)[1]
}

prepare_feature_df <- function(dt, id_col = NULL, keep_samples, type_label, min_shared_samples = 4) {
  if (is.null(id_col) || !(id_col %in% colnames(dt))) {
    id_col <- colnames(dt)[1]
  }
  
  samples <- intersect(colnames(dt), keep_samples)
  
  if (length(samples) < min_shared_samples) {
    warning(
      "Skipping ", type_label,
      ": fewer than ", min_shared_samples,
      " shared samples found."
    )
    return(NULL)
  }
  
  df <- dt[, c(id_col, samples), drop = FALSE]
  
  feature_label_raw <- as.character(df[[id_col]])
  feature_label <- clean_label(feature_label_raw, max_chars = Inf)
  feature_label <- humanise_feature_label(feature_label, type = type_label, max_chars = Inf)
  
  feature_key <- paste0(type_label, "::", make_unique_safe(feature_label))
  
  rownames(df) <- feature_key
  df[[id_col]] <- NULL
  
  df <- as.data.frame(df, check.names = FALSE)
  df[] <- lapply(df, to_numeric_safe)
  
  lookup <- tibble(
    feature_key = feature_key,
    feature_label = feature_label,
    type = type_label
  )
  
  list(matrix = df[, samples, drop = FALSE], lookup = lookup)
}

median_impute_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (all(is.na(v))) {
      v[] <- 0
    } else {
      med <- median(v, na.rm = TRUE)
      v[is.na(v)] <- med
    }
    x[, j] <- v
  }
  
  x
}

build_rf_matrix <- function(feature_mats, feature_lookup, Y, min_non_missing = 2) {
  
  feature_mats <- feature_mats[!vapply(feature_mats, is.null, logical(1))]
  
  feature_mats <- lapply(feature_mats, function(m) {
    m <- as.data.frame(m, check.names = FALSE)
    m[] <- lapply(m, to_numeric_safe)
    m
  })
  
  shared_samples <- Reduce(intersect, lapply(feature_mats, colnames))
  shared_samples <- intersect(shared_samples, names(Y))
  
  if (length(shared_samples) < 4) {
    stop(
      "Fewer than 4 shared samples across feature matrices and Y. ",
      "Check sample names in taxa/reaction/subsystem/metabolite matrices."
    )
  }
  
  message("\nRF build: shared samples across all feature classes and Y: ", length(shared_samples))
  
  feature_mats <- lapply(feature_mats, function(m) {
    m[, shared_samples, drop = FALSE]
  })
  
  features_df <- do.call(rbind, unname(feature_mats))
  features_df <- as.data.frame(features_df, check.names = FALSE)
  features_df[] <- lapply(features_df, to_numeric_safe)
  
  rf_lookup <- tibble(feature_key = rownames(features_df)) %>%
    left_join(feature_lookup, by = "feature_key") %>%
    mutate(
      feature_label = coalesce(feature_label, feature_key),
      type = coalesce(type, case_when(
        str_detect(feature_key, "^Taxa::") ~ "Taxa",
        str_detect(feature_key, "^Reaction::") ~ "Reaction",
        str_detect(feature_key, "^Subsystem::") ~ "Subsystem",
        str_detect(feature_key, "^Metabolite::") ~ "Metabolite",
        TRUE ~ "Unknown"
      )),
      feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
      model_feature = paste0("F", sprintf("%05d", row_number()))
    )
  
  rownames(features_df) <- rf_lookup$model_feature
  
  X_raw <- t(as.matrix(features_df))
  storage.mode(X_raw) <- "numeric"
  
  Y2 <- Y[rownames(X_raw)]
  
  valid_y <- !is.na(Y2) & is.finite(Y2)
  X_raw <- X_raw[valid_y, , drop = FALSE]
  Y2 <- Y2[valid_y]
  
  feature_non_missing <- colSums(!is.na(X_raw))
  feature_sd_pre <- apply(X_raw, 2, sd, na.rm = TRUE)
  
  keep_pre <- feature_non_missing >= min_non_missing
  X_pre <- X_raw[, keep_pre, drop = FALSE]
  
  if (ncol(X_pre) < 2) {
    warning("RF pre-filter retained fewer than 2 features. Relaxing to any non-missing feature.")
    keep_pre <- feature_non_missing >= 1
    X_pre <- X_raw[, keep_pre, drop = FALSE]
  }
  
  X_imp <- median_impute_matrix(X_pre)
  
  feature_sd_post <- apply(X_imp, 2, sd, na.rm = TRUE)
  keep_post <- is.finite(feature_sd_post) & feature_sd_post > 0
  
  X_imp <- X_imp[, keep_post, drop = FALSE]
  
  rf_lookup <- rf_lookup %>%
    filter(model_feature %in% colnames(X_imp)) %>%
    mutate(
      feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
      type = factor(type, levels = feature_type_levels_with_unknown)
    )
  
  filter_diag <- tibble(
    model_feature = colnames(X_raw),
    non_missing_n = feature_non_missing,
    sd_pre_imputation = feature_sd_pre,
    kept_pre_imputation = colnames(X_raw) %in% colnames(X_pre),
    kept_post_imputation = colnames(X_raw) %in% colnames(X_imp)
  ) %>%
    left_join(
      rf_lookup %>% select(model_feature, feature_key, feature_label, type),
      by = "model_feature"
    )
  
  list(
    X = as.data.frame(X_imp, check.names = FALSE),
    Y = Y2,
    rf_lookup = rf_lookup,
    filter_diag = filter_diag,
    n_features_raw = ncol(X_raw),
    n_features_after = ncol(X_imp),
    n_samples = nrow(X_imp),
    shared_samples = rownames(X_imp)
  )
}

# ============================================================
# 5. CHECK INPUT OBJECT
# ============================================================

if (!exists("vst_physeq")) {
  stop("vst_physeq not found. Load or compute it first.")
}

vst_physeq_subset <- subset_samples(
  vst_physeq,
  Condition %in% c("PFAS_Exposed", "Reference")
)

# ============================================================
# 6. COMPUTE AVERAGE PFAS ELIMINATION RATE
# ============================================================

pfas_cols <- c(
  "k_PFOA", "k_PFPeS", "k_PFHxS", "k_PFHpS",
  "k_LPFOS", "k_PFOS_MP1", "k_PFOS_MP345", "k_PFOS_MP26"
)

sdat <- as.data.frame(sample_data(vst_physeq_subset), check.names = FALSE)

present_pfas_cols <- intersect(pfas_cols, colnames(sdat))
if (length(present_pfas_cols) == 0) {
  stop("No PFAS elimination-rate columns found in sample_data().")
}

for (cc in present_pfas_cols) {
  sdat[[cc]] <- to_numeric_safe(sdat[[cc]])
}

sdat$average_k <- rowMeans(
  sdat[, present_pfas_cols, drop = FALSE],
  na.rm = TRUE
)

sdat$average_k[
  rowSums(!is.na(sdat[, present_pfas_cols, drop = FALSE])) == 0
] <- NA_real_

sample_data(vst_physeq_subset)$average_k <- sdat$average_k

pfas_rates <- as.numeric(sample_data(vst_physeq_subset)$average_k)
names(pfas_rates) <- sample_names(vst_physeq_subset)
pfas_rates <- pfas_rates[!is.na(pfas_rates)]

if (length(pfas_rates) < min_shared_samples) {
  stop("Too few samples with non-missing average_k.")
}

keep_samples <- names(pfas_rates)

message("\nSamples with non-missing average_k: ", length(keep_samples))

# ============================================================
# 7. FEATURE MATRICES
# ============================================================

otu_tab <- as.data.frame(otu_table(vst_physeq_subset), check.names = FALSE)

if (!taxa_are_rows(vst_physeq_subset)) {
  otu_tab <- t(otu_tab)
  otu_tab <- as.data.frame(otu_tab, check.names = FALSE)
}

otu_tab <- otu_tab[, keep_samples, drop = FALSE]
otu_tab[] <- lapply(otu_tab, to_numeric_safe)

taxa_labels <- get_taxa_labels(vst_physeq_subset, rownames(otu_tab))
taxa_labels <- humanise_feature_label(taxa_labels, type = "Taxa", max_chars = Inf)

taxa_feature_key <- paste0("Taxa::", make_unique_safe(taxa_labels))
rownames(otu_tab) <- taxa_feature_key

taxa_lookup <- tibble(
  feature_key = taxa_feature_key,
  feature_label = taxa_labels,
  type = "Taxa"
)

reaction_dt <- fread(reaction_path, data.table = FALSE, check.names = FALSE)
subsystem_dt <- fread(subsystem_path, data.table = FALSE, check.names = FALSE)
metabolite_dt <- fread(metabolic_data_path, data.table = FALSE, check.names = FALSE)

reaction_id_col <- guess_id_col(
  reaction_dt,
  preferred_names = c("Reaction", "Reactions", "reaction", "reaction_id")
)

subsystem_id_col <- guess_id_col(
  subsystem_dt,
  preferred_names = c("Subsystems", "Subsystem", "subsystem")
)

metabolite_id_col <- guess_id_col(
  metabolite_dt,
  preferred_names = c("Metabolite", "Metabolites", "metabolite", "metabolite_id", "Name", "ID")
)

reaction_prepared <- prepare_feature_df(
  reaction_dt,
  id_col = reaction_id_col,
  keep_samples = keep_samples,
  type_label = "Reaction",
  min_shared_samples = min_shared_samples
)

subsystem_prepared <- prepare_feature_df(
  subsystem_dt,
  id_col = subsystem_id_col,
  keep_samples = keep_samples,
  type_label = "Subsystem",
  min_shared_samples = min_shared_samples
)

metabolite_prepared <- prepare_feature_df(
  metabolite_dt,
  id_col = metabolite_id_col,
  keep_samples = keep_samples,
  type_label = "Metabolite",
  min_shared_samples = min_shared_samples
)

feature_mats <- list(Taxa = otu_tab)
feature_lookup <- taxa_lookup

if (!is.null(reaction_prepared)) {
  feature_mats$Reaction <- reaction_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, reaction_prepared$lookup)
}

if (!is.null(subsystem_prepared)) {
  feature_mats$Subsystem <- subsystem_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, subsystem_prepared$lookup)
}

if (!is.null(metabolite_prepared)) {
  feature_mats$Metabolite <- metabolite_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, metabolite_prepared$lookup)
}

message("\nFeature matrix sample-column diagnostics before RF build:")
for (nm in names(feature_mats)) {
  message("  ", nm, ": ", ncol(feature_mats[[nm]]), " sample columns")
}
message("  Shared across all feature matrices: ", length(Reduce(intersect, lapply(feature_mats, colnames))))

feature_lookup <- feature_lookup %>%
  distinct(feature_key, .keep_all = TRUE) %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = factor(type, levels = feature_type_levels)
  )

message("\nFeatures loaded:")
message("  Taxa: ", nrow(feature_mats$Taxa))
if (!is.null(feature_mats$Reaction)) message("  Reactions: ", nrow(feature_mats$Reaction))
if (!is.null(feature_mats$Subsystem)) message("  Subsystems: ", nrow(feature_mats$Subsystem))
if (!is.null(feature_mats$Metabolite)) message("  Metabolites: ", nrow(feature_mats$Metabolite))

# ============================================================
# 8. SPEARMAN CORRELATION ANALYSIS
# ============================================================

cor_list <- lapply(names(feature_mats), function(tp) {
  mat_i <- feature_mats[[tp]]
  
  samples_i <- intersect(colnames(mat_i), names(pfas_rates))
  mat_i <- mat_i[, samples_i, drop = FALSE]
  
  compute_feature_cor(mat_i, pfas_rates, min_pairs = min_pairs_cor) %>%
    mutate(type = tp)
})

all_cor <- bind_rows(cor_list) %>%
  left_join(feature_lookup, by = c("feature_key", "type")) %>%
  mutate(
    feature_label = coalesce(feature_label, feature_key),
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    abs_rho = abs(rho),
    p_adj_BH = p.adjust(p, method = "BH"),
    direction = ifelse(rho > 0, "Positive", "Negative"),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  arrange(type, desc(abs_rho))

write_csv(all_cor, out_csv_all)

sig_all <- all_cor %>%
  filter(abs_rho >= rho_threshold, p <= pval_threshold)

write_csv(sig_all, out_csv_sig)

if (nrow(sig_all) == 0) {
  stop("No features passed the rho/p-value thresholds for Figure 6.")
}

message("\nCorrelation features passing plot filter: ", nrow(sig_all))

# ============================================================
# 9. NETWORK NODES AND EDGES
# ============================================================

net_nodes <- sig_all %>%
  group_by(type) %>%
  slice_max(abs_rho, n = max_nodes_per_type, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(node_name = make_unique_safe(feature_label))

nodes <- distinct(bind_rows(
  tibble(
    name = net_nodes$node_name,
    label = net_nodes$feature_label,
    type = as.character(net_nodes$type)
  ),
  tibble(
    name = "PFAS_Elimination_Rate",
    label = "PFAS elimination rate",
    type = "PFAS"
  )
))

edges <- tibble(
  from = "PFAS_Elimination_Rate",
  to = net_nodes$node_name,
  weight = net_nodes$abs_rho,
  rho = net_nodes$rho,
  p = net_nodes$p,
  p_adj_BH = net_nodes$p_adj_BH,
  n_pairs = net_nodes$n_pairs,
  sign = ifelse(net_nodes$rho > 0, "pos", "neg")
)

write_csv(nodes, out_csv_nodes)
write_csv(edges, out_csv_edges)

node_stats <- edges %>%
  group_by(to) %>%
  summarise(
    mean_abs_rho = mean(abs(rho), na.rm = TRUE),
    mean_p = mean(p, na.rm = TRUE),
    node_size = -log10(mean_p + 1e-10),
    .groups = "drop"
  )

net <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

set.seed(123)
layout_network <- create_layout(net, layout = "fr") %>%
  left_join(node_stats, by = c("name" = "to")) %>%
  mutate(node_size = ifelse(is.na(node_size), 2.5, node_size))

# ============================================================
# 10. COLORS
# ============================================================

col_type <- c(
  "Taxa" = "purple3",
  "Reaction" = "deepskyblue3",
  "Subsystem" = "slateblue3",
  "Metabolite" = "darkorange2",
  "PFAS" = "red3",
  "Unknown" = "grey60"
)

col_edge_sign <- c(
  "pos" = "#B22222",
  "neg" = "#2E8B57"
)

# ============================================================
# 11. MAIN FIGURE 6A/B
# ============================================================

p_summary <- sig_all %>%
  group_by(type) %>%
  slice_max(abs_rho, n = top_summary_per_type, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    type = factor(type, levels = feature_type_levels),
    feature_label_plot = humanise_feature_label(feature_label, type = type, max_chars = 65),
    feature_label_plot = fct_reorder(feature_label_plot, abs_rho)
  ) %>%
  ggplot(aes(x = feature_label_plot, y = abs_rho, fill = type)) +
  geom_col(alpha = 0.95, colour = "black", width = 0.54) +
  coord_flip() +
  facet_wrap(~ type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = "Feature",
    y = expression("|Spearman " * rho * "|")
  ) +
  theme_clean(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 11.2, face = "bold", margin = margin(r = 6)),
    axis.text.x = element_text(size = 12.5, face = "bold"),
    axis.title.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.6, "lines"),
    plot.margin = margin(12, 10, 12, 10)
  )

p_network <- ggraph(layout_network) +
  geom_edge_link(
    aes(width = weight, colour = sign),
    alpha = 0.8,
    lineend = "round"
  ) +
  scale_edge_color_manual(
    values = col_edge_sign,
    name = "Correlation direction",
    labels = c(
      "neg" = "Negative (ρ < 0)",
      "pos" = "Positive (ρ > 0)"
    )
  ) +
  geom_node_point(
    aes(color = type, size = node_size),
    alpha = 0.95
  ) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 3.7,
    color = "black",
    fontface = "bold",
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = col_type,
    name = "Feature class",
    labels = c(
      "PFAS" = "PFAS elimination rate",
      "Reaction" = "Metabolic reaction",
      "Subsystem" = "Functional subsystem",
      "Taxa" = "Microbial taxon",
      "Metabolite" = "Model-predicted metabolite"
    )
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(p))),
    range = c(3.2, 9.0)
  ) +
  scale_edge_width(
    range = c(0.60, 3.4),
    name = "|ρ|"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 5)),
    size = guide_legend(order = 2),
    edge_colour = guide_legend(order = 3),
    edge_width = guide_legend(order = 4)
  ) +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11.5),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 12. BUILD RF MATRIX
# ============================================================

rf_build <- build_rf_matrix(
  feature_mats = feature_mats,
  feature_lookup = feature_lookup,
  Y = pfas_rates,
  min_non_missing = 2
)

X <- rf_build$X
Y <- rf_build$Y
rf_lookup <- rf_build$rf_lookup

write_csv(rf_build$filter_diag, csv_rf_filtering)

message("\nRF matrix diagnostic:")
message("  Samples: ", rf_build$n_samples)
message("  Raw features before RF filtering: ", rf_build$n_features_raw)
message("  Features after imputation and variance filtering: ", rf_build$n_features_after)
message("  RF refits planned: ", n_runs)

if (ncol(X) < 2) {
  stop(
    paste0(
      "RF still has fewer than 2 non-zero variance features. ",
      "Check diagnostic CSV: ", csv_rf_filtering
    )
  )
}

# ============================================================
# 13. RF RUN FUNCTION
# ============================================================

run_rf_once <- function(run_id, X, Y, train_fraction, seed,
                        cv_folds, cv_repeats,
                        num_trees, min_node_size) {
  
  set.seed(seed)
  
  train_idx <- createDataPartition(Y, p = train_fraction, list = FALSE)
  
  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  
  X_test <- X[-train_idx, , drop = FALSE]
  Y_test <- Y[-train_idx]
  
  train_ctrl <- trainControl(
    method = if (cv_repeats > 1) "repeatedcv" else "cv",
    number = cv_folds,
    repeats = cv_repeats
  )
  
  tgrid <- expand.grid(
    mtry = max(1, floor(sqrt(ncol(X_train)))),
    splitrule = "variance",
    min.node.size = min_node_size
  )
  
  rf_model <- train(
    x = X_train,
    y = Y_train,
    method = "ranger",
    trControl = train_ctrl,
    tuneGrid = tgrid,
    importance = "permutation",
    num.trees = num_trees,
    metric = "RMSE"
  )
  
  pred_test <- predict(rf_model, X_test)
  resid_test <- as.numeric(Y_test - pred_test)
  
  rmse_test <- sqrt(mean(resid_test^2, na.rm = TRUE))
  r2_test <- suppressWarnings(cor(Y_test, pred_test, use = "complete.obs")^2)
  
  sd_resid <- sd(resid_test, na.rm = TRUE)
  std_resid <- if (!is.finite(sd_resid) || sd_resid == 0) {
    rep(NA_real_, length(resid_test))
  } else {
    resid_test / sd_resid
  }
  
  shapiro_p <- if (length(resid_test) >= 3 && length(resid_test) <= 5000) {
    suppressWarnings(shapiro.test(resid_test)$p.value)
  } else {
    NA_real_
  }
  
  perf <- tibble(
    run = run_id,
    seed = seed,
    n_train = nrow(X_train),
    n_test = nrow(X_test),
    rmse = rmse_test,
    r2 = as.numeric(r2_test)
  )
  
  resid_stats <- tibble(
    run = run_id,
    resid_mean = mean(resid_test, na.rm = TRUE),
    resid_sd = sd(resid_test, na.rm = TRUE),
    shapiro_p = shapiro_p
  )
  
  var_imp <- varImp(rf_model, scale = TRUE)$importance
  
  imp_tbl <- tibble(
    run = run_id,
    model_feature = rownames(var_imp),
    importance = as.numeric(var_imp[[1]])
  ) %>%
    arrange(desc(importance))
  
  test_df <- tibble(
    run = run_id,
    sample = rownames(X_test),
    obs = as.numeric(Y_test),
    pred = as.numeric(pred_test),
    resid = resid_test,
    std_resid = std_resid
  )
  
  list(
    perf = perf,
    resid_stats = resid_stats,
    importance = imp_tbl,
    test_df = test_df
  )
}

# ============================================================
# 14. RUN RF REPEATED REFITS
# ============================================================

message("\n============================================================")
message("Starting RF repeated-refit analysis")
message("Total RF refits: ", n_runs)
message("============================================================")
flush.console()

run_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  message(sprintf("RF refit %03d / %03d", i, n_runs))
  flush.console()
  
  run_results[[i]] <- run_rf_once(
    run_id = i,
    X = X,
    Y = Y,
    train_fraction = train_fraction,
    seed = seed_base + i,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    num_trees = num_trees,
    min_node_size = min_node_size
  )
}

message("Finished RF repeated-refit analysis.")
flush.console()

perf_df <- bind_rows(lapply(run_results, `[[`, "perf"))
resid_stats_df <- bind_rows(lapply(run_results, `[[`, "resid_stats"))
test_all_df <- bind_rows(lapply(run_results, `[[`, "test_df"))

imp_all_df <- bind_rows(lapply(run_results, `[[`, "importance")) %>%
  left_join(
    rf_lookup %>%
      select(model_feature, feature_key, feature_label, type),
    by = "model_feature"
  ) %>%
  mutate(
    feature_label = coalesce(feature_label, model_feature),
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = coalesce(as.character(type), "Unknown"),
    type = factor(type, levels = feature_type_levels_with_unknown)
  )

# ============================================================
# 15. RF STABILITY SUMMARY
# ============================================================

imp_topk_df <- imp_all_df %>%
  group_by(run) %>%
  arrange(desc(importance), .by_group = TRUE) %>%
  slice_head(n = top_k) %>%
  ungroup() %>%
  mutate(selected_topk = 1L)

stability_df <- imp_all_df %>%
  group_by(model_feature, feature_key, feature_label, type) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    sd_importance = sd(importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    imp_topk_df %>%
      group_by(model_feature) %>%
      summarise(topk_count = sum(selected_topk), .groups = "drop"),
    by = "model_feature"
  ) %>%
  mutate(
    topk_count = replace_na(topk_count, 0L),
    topk_freq = topk_count / n_runs,
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = factor(type, levels = feature_type_levels_with_unknown)
  ) %>%
  arrange(desc(topk_freq), desc(mean_importance))

top_stable_per_type <- stability_df %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = 10, with_ties = FALSE) %>%
  ungroup()

top_display <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = top_rf_features_per_type_main, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(type, desc(topk_freq), desc(mean_importance))

stable_candidates <- stability_df %>%
  filter(topk_freq >= 0.25) %>%
  arrange(desc(topk_freq), desc(mean_importance))

# ============================================================
# 15B. METABOLITE TRANSLATION DIAGNOSTIC
# ============================================================

metabolite_labels_needing_translation <- stability_df %>%
  filter(type == "Metabolite") %>%
  mutate(
    current_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    needs_translation = is_vmh_metabolite_like_label(current_label)
  ) %>%
  filter(needs_translation) %>%
  distinct(
    model_feature,
    feature_key,
    current_label,
    topk_freq,
    mean_importance
  ) %>%
  arrange(desc(topk_freq), desc(mean_importance), current_label)

write_csv(
  metabolite_labels_needing_translation,
  csv_metabolite_untranslated
)

if (nrow(metabolite_labels_needing_translation) > 0) {
  writeLines(
    paste0('"', metabolite_labels_needing_translation$current_label, '" = "",'),
    txt_metabolite_template
  )
} else {
  writeLines(
    "No VMH-like metabolite labels detected after current mapping.",
    txt_metabolite_template
  )
}

cat("\n============================================================\n")
cat("METABOLITE LABELS STILL NEEDING TRANSLATION\n")
cat("============================================================\n")
print(metabolite_labels_needing_translation, n = Inf)

# ============================================================
# 16. WRITE RF CSV OUTPUTS
# ============================================================

write_csv(perf_df, csv_perf)
write_csv(resid_stats_df, csv_resid_stats)
write_csv(imp_all_df, csv_importance_all)
write_csv(imp_topk_df, csv_selected_topk)
write_csv(stability_df, csv_stability)
write_csv(top_stable_per_type, csv_topstable_type)
write_csv(top_display, csv_topdisplay)
write_csv(stable_candidates, csv_stablecand)

supp_stability_summary <- tibble(
  n_runs = n_runs,
  train_fraction = train_fraction,
  top_k = top_k,
  rmse_median = median(perf_df$rmse, na.rm = TRUE),
  rmse_iqr25 = quantile(perf_df$rmse, 0.25, na.rm = TRUE),
  rmse_iqr75 = quantile(perf_df$rmse, 0.75, na.rm = TRUE),
  r2_median = median(perf_df$r2, na.rm = TRUE),
  r2_iqr25 = quantile(perf_df$r2, 0.25, na.rm = TRUE),
  r2_iqr75 = quantile(perf_df$r2, 0.75, na.rm = TRUE),
  n_features_total = nrow(stability_df),
  n_features_selected_at_least_once = sum(stability_df$topk_count > 0, na.rm = TRUE),
  n_features_selected_in_ge10pct_runs = sum(stability_df$topk_freq >= 0.10, na.rm = TRUE),
  n_features_selected_in_ge25pct_runs = sum(stability_df$topk_freq >= 0.25, na.rm = TRUE),
  n_features_selected_in_ge50pct_runs = sum(stability_df$topk_freq >= 0.50, na.rm = TRUE)
)

write_csv(supp_stability_summary, csv_stab_overall)

supp_stability_by_type <- stability_df %>%
  group_by(type) %>%
  summarise(
    n_features = n(),
    n_selected_ge10pct = sum(topk_freq >= 0.10, na.rm = TRUE),
    n_selected_ge25pct = sum(topk_freq >= 0.25, na.rm = TRUE),
    n_selected_ge50pct = sum(topk_freq >= 0.50, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(supp_stability_by_type, csv_stab_bytype)

# ============================================================
# 17. MAIN FIGURE 6C — TOP STABLE RF PREDICTORS PER TYPE
# ============================================================

main_rf_plot_df <- top_display %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    feature_wrapped = str_wrap(feature_label, width = 34),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  group_by(type) %>%
  arrange(topk_freq, .by_group = TRUE) %>%
  mutate(feature_wrapped = factor(feature_wrapped, levels = unique(feature_wrapped))) %>%
  ungroup()

p_rf_stable_main <- ggplot(
  main_rf_plot_df,
  aes(x = topk_freq, y = feature_wrapped, fill = type)
) +
  geom_col(width = 0.70, color = "black", linewidth = 0.25) +
  facet_wrap(~ type, scales = "free_y", ncol = 1, drop = TRUE) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    x = sprintf("Selection frequency\n(top-%d across %d RF refits)", top_k, n_runs),
    y = "Predictive feature"
  ) +
  theme_clean(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13.5),
    axis.text.y = element_text(size = 10.2, face = "bold", lineheight = 0.90),
    axis.text.x = element_text(size = 12.0, face = "bold"),
    axis.title.x = element_text(size = 13.5),
    axis.title.y = element_text(size = 13.5),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.0, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 18. MAIN FIGURE 6D — TOP STABLE FEATURE HEATMAP
# ============================================================

readable_stability <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    readable_label = !(type == "Reaction" & is_vmh_like_label(feature_label)),
    priority = case_when(
      type != "Reaction" ~ 2,
      readable_label ~ 1,
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(priority), desc(topk_freq), desc(mean_importance))

top_heatmap_features <- readable_stability %>%
  slice_head(n = top_rf_features_heatmap) %>%
  pull(model_feature)

top_heatmap_features <- intersect(top_heatmap_features, colnames(X))

message("\nTop heatmap features retained: ", length(top_heatmap_features))
print(top_heatmap_features)

if (length(top_heatmap_features) >= 2) {
  
  X_heat <- X[, top_heatmap_features, drop = FALSE]
  
  cor_mat <- suppressWarnings(
    cor(
      X_heat,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  rownames(cor_mat) <- top_heatmap_features
  colnames(cor_mat) <- top_heatmap_features
  
  heatmap_labels <- rf_lookup %>%
    filter(model_feature %in% top_heatmap_features) %>%
    mutate(
      feature_label_full = humanise_feature_label(
        feature_label,
        type = type,
        max_chars = Inf
      ),
      feature_label_wrapped_raw = str_wrap(feature_label_full, width = 22),
      type = factor(type, levels = feature_type_levels),
      model_feature = factor(model_feature, levels = top_heatmap_features)
    ) %>%
    arrange(model_feature) %>%
    mutate(
      model_feature = as.character(model_feature),
      feature_label_wrapped = make.unique(
        as.character(feature_label_wrapped_raw),
        sep = " "
      )
    )
  
  if (nrow(heatmap_labels) != length(top_heatmap_features)) {
    warning(
      "Some top_heatmap_features were not found in rf_lookup. ",
      "Proceeding with matched features only."
    )
    
    top_heatmap_features <- intersect(
      top_heatmap_features,
      heatmap_labels$model_feature
    )
    
    cor_mat <- cor_mat[
      top_heatmap_features,
      top_heatmap_features,
      drop = FALSE
    ]
  }
  
  heatmap_order <- heatmap_labels$feature_label_wrapped[
    match(top_heatmap_features, heatmap_labels$model_feature)
  ]
  
  heatmap_df_raw <- as.data.frame(
    as.table(cor_mat),
    stringsAsFactors = FALSE
  )
  
  names(heatmap_df_raw)[seq_len(3)] <- c(
    "Feature1",
    "Feature2",
    "Spearman_rho"
  )
  
  heatmap_df <- heatmap_df_raw %>%
    mutate(
      Feature1 = as.character(Feature1),
      Feature2 = as.character(Feature2),
      Spearman_rho = as.numeric(Spearman_rho)
    ) %>%
    left_join(
      heatmap_labels %>%
        select(model_feature, Label1 = feature_label_wrapped),
      by = c("Feature1" = "model_feature")
    ) %>%
    left_join(
      heatmap_labels %>%
        select(model_feature, Label2 = feature_label_wrapped),
      by = c("Feature2" = "model_feature")
    ) %>%
    filter(!is.na(Label1), !is.na(Label2)) %>%
    mutate(
      Label1 = factor(Label1, levels = heatmap_order),
      Label2 = factor(Label2, levels = rev(heatmap_order)),
      x_num = as.numeric(Label1),
      y_num = as.numeric(Label2)
    ) %>%
    filter(!is.na(x_num), !is.na(y_num))
  
  label_position_df <- heatmap_labels %>%
    mutate(
      Label = factor(feature_label_wrapped, levels = heatmap_order),
      x_num = as.numeric(Label),
      y_num = length(heatmap_order) - x_num + 1
    ) %>%
    filter(!is.na(x_num), !is.na(y_num))
  
  row_type_boxes <- label_position_df %>%
    transmute(
      x = 0.25,
      y = y_num,
      type = type
    )
  
  col_type_boxes <- label_position_df %>%
    transmute(
      x = x_num,
      y = length(heatmap_order) + 0.75,
      type = type
    )
  
  write_csv(heatmap_df, csv_rf_feature_heatmap)
  
  p_rf_heatmap_main <- ggplot() +
    geom_tile(
      data = heatmap_df,
      aes(
        x = x_num,
        y = y_num,
        fill = Spearman_rho
      ),
      color = "white",
      linewidth = 0.25
    ) +
    geom_tile(
      data = row_type_boxes,
      aes(
        x = x,
        y = y,
        color = type
      ),
      fill = "white",
      width = 0.34,
      height = 0.78,
      linewidth = 1.35,
      show.legend = FALSE
    ) +
    geom_tile(
      data = col_type_boxes,
      aes(
        x = x,
        y = y,
        color = type
      ),
      fill = "white",
      width = 0.78,
      height = 0.34,
      linewidth = 1.35,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman ρ"
    ) +
    scale_color_manual(
      values = col_type,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq_along(heatmap_order),
      labels = heatmap_order,
      expand = expansion(mult = c(0.08, 0.03))
    ) +
    scale_y_continuous(
      breaks = seq_along(rev(heatmap_order)),
      labels = rev(heatmap_order),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Predictive feature",
      y = "Predictive feature"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 48,
        hjust = 1,
        vjust = 1,
        face = "bold",
        size = 9.2,
        lineheight = 0.86
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 9.2,
        lineheight = 0.86
      ),
      axis.title = element_text(
        face = "bold",
        size = 14
      ),
      legend.title = element_text(
        face = "bold",
        size = 12
      ),
      legend.text = element_text(
        size = 11
      ),
      panel.grid = element_blank(),
      plot.margin = margin(12, 12, 12, 16)
    )
  
} else {
  
  p_rf_heatmap_main <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Too few stable features for heatmap",
      fontface = "bold"
    ) +
    theme_void()
}

# ============================================================
# 19. CHECK ONLY PLOTTED FIGURE 6C/D LABELS FOR TRANSLATION
# ============================================================

plotted_rf_labels_for_translation <- bind_rows(
  main_rf_plot_df %>%
    mutate(source = "Figure6C") %>%
    select(source, type, model_feature, feature_key, feature_label),
  readable_stability %>%
    filter(model_feature %in% top_heatmap_features) %>%
    mutate(source = "Figure6D") %>%
    select(source, type, model_feature, feature_key, feature_label)
) %>%
  distinct(source, type, model_feature, feature_key, feature_label)

plotted_untranslated <- get_plotted_reaction_labels_needing_translation(
  plotted_rf_labels_for_translation,
  out_dir = out_dir
)

cat("\n============================================================\n")
cat("PLOTTED reaction labels still needing translation\n")
cat("============================================================\n")
print(plotted_untranslated, n = Inf)

# ============================================================
# 20. SAVE MAIN FIGURE 6 A-D
# ============================================================

figure6_final <- (p_summary | p_network) / (p_rf_stable_main | p_rf_heatmap_main) +
  plot_layout(
    widths = c(1.18, 1.22),
    heights = c(1.08, 1.18),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 30)
    )
  ) &
  theme(legend.position = "right")

ggsave(
  out_fig6_png,
  plot = figure6_final,
  width = 22,
  height = 19.2,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  out_fig6_pdf,
  plot = figure6_final,
  width = 22,
  height = 19.2,
  bg = "white",
  limitsize = FALSE
)

safe_print(figure6_final)

message("\nSaved main Figure 6:")
message("  ", out_fig6_png)
message("  ", out_fig6_pdf)

# ============================================================
# 21. SUPPLEMENTARY FIGURE 1 — RF QC
# ============================================================

rmse_med <- median(perf_df$rmse, na.rm = TRUE)
r2_med <- median(perf_df$r2, na.rm = TRUE)

resid_mean_med <- median(resid_stats_df$resid_mean, na.rm = TRUE)
resid_sd_med <- median(resid_stats_df$resid_sd, na.rm = TRUE)
shap_med <- median(resid_stats_df$shapiro_p, na.rm = TRUE)

std_resid_cor <- suppressWarnings(
  cor.test(
    test_all_df$obs,
    test_all_df$std_resid,
    method = "spearman",
    exact = FALSE
  )
)

std_resid_label <- paste0(
  "Spearman ρ = ", round(as.numeric(std_resid_cor$estimate), 2),
  "\np = ", signif(std_resid_cor$p.value, 2)
)

p_qc1 <- ggplot(test_all_df, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.12, size = 1, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE, colour = "red") +
  theme_bw(base_size = 12) +
  labs(
    title = "Observed vs predicted",
    x = "Observed PFAS elimination",
    y = "Predicted PFAS elimination"
  ) +
  annotate(
    "text",
    x = min(test_all_df$obs, na.rm = TRUE),
    y = max(test_all_df$pred, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0(
      "Median across runs:\n",
      "RMSE = ", round(rmse_med, 3), "\n",
      "R² = ", round(r2_med, 3)
    )
  )

p_qc2 <- ggplot(test_all_df, aes(x = pred, y = resid)) +
  geom_point(alpha = 0.12, size = 1, colour = "darkorange") +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 12) +
  labs(
    title = "Residuals vs predicted",
    x = "Predicted PFAS elimination",
    y = "Residual (observed − predicted)"
  ) +
  annotate(
    "text",
    x = min(test_all_df$pred, na.rm = TRUE),
    y = max(test_all_df$resid, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0(
      "Median across runs:\n",
      "Mean(resid) = ", round(resid_mean_med, 3), "\n",
      "SD(resid) = ", round(resid_sd_med, 3)
    )
  )

p_qc3 <- ggplot(test_all_df, aes(x = obs, y = std_resid)) +
  geom_point(alpha = 0.12, size = 1, colour = "purple3") +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate(
    "text",
    x = min(test_all_df$obs, na.rm = TRUE),
    y = max(test_all_df$std_resid, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = std_resid_label
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Standardized residuals vs observed",
    x = "Observed PFAS elimination",
    y = "Standardized residual"
  )

p_qc4 <- ggplot(test_all_df, aes(x = resid)) +
  geom_density(fill = "grey70", alpha = 0.7, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  theme_bw(base_size = 12) +
  labs(
    title = "Residual density",
    x = "Residual",
    y = "Density"
  ) +
  annotate(
    "text",
    x = min(test_all_df$resid, na.rm = TRUE),
    y = Inf,
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0("Median Shapiro p = ", signif(shap_med, 3))
  )

qc_4panel <- (p_qc1 | p_qc2) / (p_qc3 | p_qc4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_qc_4panel,
  plot = qc_4panel,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_qc_4panel_pdf,
  plot = qc_4panel,
  width = 12,
  height = 10,
  bg = "white",
  limitsize = FALSE
)

safe_print(qc_4panel)

# ============================================================
# 22. SUPPLEMENTARY FIGURE 2 — RF STABILITY DIAGNOSTICS
# ============================================================

rmse_iqr <- quantile(perf_df$rmse, c(0.25, 0.75), na.rm = TRUE)
r2_iqr <- quantile(perf_df$r2, c(0.25, 0.75), na.rm = TRUE)

n_selected_once <- sum(stability_df$topk_count > 0, na.rm = TRUE)
n_selected_10 <- sum(stability_df$topk_freq >= 0.10, na.rm = TRUE)
n_selected_25 <- sum(stability_df$topk_freq >= 0.25, na.rm = TRUE)

p_stab1 <- ggplot(perf_df, aes(x = "", y = rmse)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey85", colour = "black") +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.2) +
  annotate(
    "text",
    x = 1.28,
    y = max(perf_df$rmse, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Median = ", round(rmse_med, 3),
      "\nIQR = ", round(rmse_iqr[1], 3), "–", round(rmse_iqr[2], 3),
      "\nRuns = ", n_runs
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "RMSE across refits",
    x = NULL,
    y = "RMSE"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_stab2 <- ggplot(perf_df, aes(x = "", y = r2)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey85", colour = "black") +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.2) +
  annotate(
    "text",
    x = 1.28,
    y = max(perf_df$r2, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Median = ", round(r2_med, 3),
      "\nIQR = ", round(r2_iqr[1], 3), "–", round(r2_iqr[2], 3),
      "\nRuns = ", n_runs
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "R² across refits",
    x = NULL,
    y = "R²"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_stab3 <- ggplot(stability_df, aes(x = topk_freq)) +
  geom_histogram(bins = 25, color = "black", fill = "grey80") +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "grey35") +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "red") +
  annotate(
    "text",
    x = 0.25,
    y = Inf,
    label = "25%",
    hjust = -0.1,
    vjust = 1.4,
    size = 3.5,
    fontface = "bold",
    colour = "red"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.25,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Selected ≥1×: ", n_selected_once,
      "\n≥10% runs: ", n_selected_10,
      "\n≥25% runs: ", n_selected_25
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = sprintf("Selection frequency distribution (top-%d)", top_k),
    x = "Selection frequency",
    y = "Number of features"
  )

type_composition <- stability_df %>%
  mutate(
    stability_group = case_when(
      topk_freq >= 0.25 ~ "≥25% of refits",
      topk_freq >= 0.10 ~ "10–24% of refits",
      topk_freq > 0 ~ "Selected at least once",
      TRUE ~ "Never selected"
    ),
    stability_group = factor(
      stability_group,
      levels = c("≥25% of refits", "10–24% of refits", "Selected at least once", "Never selected")
    )
  ) %>%
  filter(stability_group != "Never selected") %>%
  dplyr::count(stability_group, type, name = "n") %>%
  group_by(stability_group) %>%
  mutate(
    total_in_group = sum(n),
    prop = n / total_in_group,
    label = ifelse(n > 0, paste0(n, "\n", round(100 * prop), "%"), "")
  ) %>%
  ungroup()

p_stab4 <- ggplot(
  type_composition,
  aes(x = stability_group, y = n, fill = type)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    title = "Feature classes among selected predictors",
    x = NULL,
    y = "Number of features"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

stability_4panel <- (p_stab1 | p_stab2) / (p_stab3 | p_stab4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_stability_4panel,
  plot = stability_4panel,
  width = 13,
  height = 10.5,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_stability_4panel_pdf,
  plot = stability_4panel,
  width = 13,
  height = 10.5,
  bg = "white",
  limitsize = FALSE
)

safe_print(stability_4panel)

# ============================================================
# 23. ADDITIONAL SUPPLEMENTARY FIGURE — TOP-STABLE FEATURES
# ============================================================

extra_top_plot_df <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    feature_wrapped = str_wrap(humanise_feature_label(feature_label, type = type, max_chars = Inf), width = 34),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  group_by(type) %>%
  arrange(topk_freq, .by_group = TRUE) %>%
  mutate(feature_wrapped = factor(feature_wrapped, levels = unique(feature_wrapped))) %>%
  ungroup()

p_extra_top <- ggplot(
  extra_top_plot_df,
  aes(x = topk_freq, y = feature_wrapped, fill = type)
) +
  geom_col(width = 0.78, color = "black", linewidth = 0.25) +
  facet_wrap(~ type, scales = "free_y", ncol = 2, drop = TRUE) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    title = "Top stable RF predictors by class",
    x = sprintf("Selection frequency (top-%d across %d RF refits)", top_k, n_runs),
    y = "Predictive feature"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

p_extra_heat <- p_rf_heatmap_main +
  labs(title = "Spearman correlations among top stable RF predictors") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

extra_rf_figure <- p_extra_top / p_extra_heat +
  plot_layout(heights = c(1.05, 1.00)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_rf_extra,
  plot = extra_rf_figure,
  width = 14.5,
  height = 13,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_rf_extra_pdf,
  plot = extra_rf_figure,
  width = 14.5,
  height = 13,
  bg = "white",
  limitsize = FALSE
)

safe_print(extra_rf_figure)

# ============================================================
# 24. FINAL MESSAGES
# ============================================================

message("\n============================================================")
message("DONE")
message("Outputs written to:")
message(out_dir)
message("============================================================")

message("\nMain Figure 6:")
message(" - ", basename(out_fig6_png))
message(" - ", basename(out_fig6_pdf))
message(" - ", basename(out_csv_all))
message(" - ", basename(out_csv_sig))
message(" - ", basename(out_csv_nodes))
message(" - ", basename(out_csv_edges))

message("\nMetabolite translation diagnostics:")
message(" - ", basename(csv_metabolite_untranslated))
message(" - ", basename(txt_metabolite_template))

message("\nRF supplementary figures:")
message(" - ", basename(fig_qc_4panel))
message(" - ", basename(fig_qc_4panel_pdf))
message(" - ", basename(fig_stability_4panel))
message(" - ", basename(fig_stability_4panel_pdf))
message(" - ", basename(fig_rf_extra))
message(" - ", basename(fig_rf_extra_pdf))

message("\nRF CSVs:")
message(" - ", basename(csv_perf))
message(" - ", basename(csv_resid_stats))
message(" - ", basename(csv_importance_all))
message(" - ", basename(csv_selected_topk))
message(" - ", basename(csv_stability))
message(" - ", basename(csv_topstable_type))
message(" - ", basename(csv_topdisplay))
message(" - ", basename(csv_stablecand))
message(" - ", basename(csv_stab_overall))
message(" - ", basename(csv_stab_bytype))
message(" - ", basename(csv_rf_filtering))
message(" - ", basename(csv_rf_feature_heatmap))

message("\nRF refits completed: ", n_runs)
