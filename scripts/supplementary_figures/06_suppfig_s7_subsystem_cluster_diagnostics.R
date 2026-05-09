############################################################
# 06_suppfig_s7_subsystem_cluster_diagnostics.R
#
# Supplementary cluster diagnostics for exposed-sample subsystem heatmap structure.
#
# Paths have been converted to repo-relative locations:
#   data/processed/
#   data/metadata/
#   data/model_outputs/
#   results/supplementary_figures/
############################################################

# in Figure 4A subsystem heatmap
#
# Fixes included:
# - Duplicate taxonomy row names handled with make.unique()
# - metadata_exp already contains Sample_Name, so rownames_to_column()
#   is NOT used when joining metadata to cluster assignments
# - Duplicate metadata column names cleaned with make.unique()
#
# This script:
# 1. Recreates exposed-sample subsystem clustering
# 2. Assigns exposed samples to two subsystem-derived clusters
# 3. Tests whether clusters associate with metadata variables
# 4. Identifies subsystem/reaction features driving the split
# 5. Checks whether the same split is visible at the species level
# 6. Saves diagnostic plots and CSV outputs
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(phyloseq)
  library(vegan)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# 1. Paths
# ============================================================

subsystem_path <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
reaction_path  <- file.path("data", "model_outputs", "ReactionAbundance.csv")
metadata_path  <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

otu_file <- file.path("data", "processed", "otu_table.tsv")
tax_file <- file.path("data", "processed", "taxonomy_table.tsv")

output_dir <- file.path("results", "supplementary_figures", "cluster_diagnostics")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Settings
# ============================================================

condition_of_interest <- "PFAS_Exposed"
k_clusters <- 2
top_n_subsystems <- 35
top_n_species <- 100
presence_threshold <- 1e-9

candidate_technical_vars <- c(
  "Sequencing_Run",
  "Run",
  "Lane",
  "Flowcell",
  "Library_Batch",
  "Extraction_Batch",
  "DNA_Extraction_Batch",
  "Sequencing_Batch",
  "Plate",
  "Date",
  "Read_Depth",
  "Total_Reads",
  "Mapped_Reads",
  "Host_Reads",
  "N50",
  "Assembly_N50",
  "Number_of_Contigs",
  "Contigs"
)

candidate_biological_vars <- c(
  "Age",
  "Sex",
  "BMI",
  "weight",
  "PFOS",
  "PFHxS",
  "PFOA",
  "PFNA",
  "sum_PFAS",
  "average_k",
  "Average_k",
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

# ============================================================
# 3. Helper functions
# ============================================================

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c(
    "", "NA", "NaN", "nd", "ND", "na", "n/a", "N/A",
    "<LOD", "<LOQ", "< LOQ", "LOD", "LOQ"
  )] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

read_clean_delim <- function(path, sep = "\t") {
  df <- data.table::fread(
    path,
    sep = sep,
    data.table = FALSE,
    check.names = FALSE
  )
  colnames(df) <- make.unique(trimws(colnames(df)))
  df
}

find_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z[is.nan(mat_z)] <- 0
  mat_z[is.infinite(mat_z)] <- 0
  mat_z
}

matrix_from_feature_csv <- function(path, feature_col_name = NULL) {
  df <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  colnames(df) <- make.unique(trimws(colnames(df)))
  
  if (is.null(feature_col_name)) {
    feature_col_name <- colnames(df)[1]
  }
  
  feature_col_index <- which(colnames(df) == feature_col_name)[1]
  if (is.na(feature_col_index)) {
    feature_col_index <- 1
  }
  
  colnames(df)[feature_col_index] <- "Feature"
  
  df$Feature <- make.unique(as.character(df$Feature))
  
  numeric_cols <- setdiff(colnames(df), "Feature")
  df[numeric_cols] <- lapply(df[numeric_cols], to_numeric_safe)
  
  mat <- df %>%
    column_to_rownames("Feature") %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  mat
}

safe_wilcox <- function(df, value_col, group_col = "Cluster") {
  dat <- df %>%
    select(all_of(c(value_col, group_col))) %>%
    mutate(
      value = to_numeric_safe(.data[[value_col]]),
      group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (n_distinct(dat$group) != 2 || nrow(dat) < 4) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "wilcox"
    ))
  }
  
  res <- tryCatch(
    wilcox.test(value ~ group, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "wilcox"
    ))
  }
  
  tibble(
    variable = value_col,
    p_value = res$p.value,
    test = "wilcox"
  )
}

safe_fisher_or_chisq <- function(df, value_col, group_col = "Cluster") {
  dat <- df %>%
    select(all_of(c(value_col, group_col))) %>%
    mutate(
      value = as.factor(.data[[value_col]]),
      group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (n_distinct(dat$group) < 2 || n_distinct(dat$value) < 2) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "categorical"
    ))
  }
  
  tab <- table(dat$group, dat$value)
  
  res <- tryCatch(
    {
      if (any(tab < 5)) fisher.test(tab) else chisq.test(tab)
    },
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "categorical"
    ))
  }
  
  tibble(
    variable = value_col,
    p_value = res$p.value,
    test = ifelse(any(tab < 5), "fisher", "chisq")
  )
}

# ============================================================
# 4. Read metadata
# ============================================================

metadata <- read_clean_delim(metadata_path, sep = "\t")

sample_col <- find_first_existing_col(
  metadata,
  c("Sample_Name", "SampleID", "Sample_ID", "Sample")
)

condition_col <- find_first_existing_col(
  metadata,
  c("Condition", "Group")
)

if (is.null(sample_col) || is.null(condition_col)) {
  stop("Could not find sample or condition column in metadata.")
}

metadata <- metadata %>%
  mutate(
    Sample_Name = trimws(as.character(.data[[sample_col]])),
    Condition = trimws(as.character(.data[[condition_col]]))
  )

# Remove duplicate Sample_Name columns if make.unique created Sample_Name.1 etc.
# Keep the standardized Sample_Name created above.
metadata <- metadata %>%
  select(-matches("^Sample_Name\\.\\d+$"))

k_cols <- c(
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

available_k_cols <- intersect(k_cols, colnames(metadata))

if (length(available_k_cols) > 0) {
  metadata$Average_k <- rowMeans(
    as.data.frame(
      lapply(metadata[, available_k_cols, drop = FALSE], to_numeric_safe)
    ),
    na.rm = TRUE
  )
} else {
  metadata$Average_k <- NA_real_
}

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata sample column used:\n")
print(sample_col)

cat("\nMetadata condition column used:\n")
print(condition_col)

# ============================================================
# 5. Read subsystem and reaction matrices
# ============================================================

subsystem_mat <- matrix_from_feature_csv(subsystem_path)
reaction_mat <- matrix_from_feature_csv(reaction_path)

cat("\nSubsystem matrix dimensions:\n")
print(dim(subsystem_mat))

cat("\nReaction matrix dimensions:\n")
print(dim(reaction_mat))

# ============================================================
# 6. Build species matrix from OTU + taxonomy files
# ============================================================

otu_df <- read.delim(
  otu_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

otu_mat <- as.matrix(otu_df)
storage.mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
rownames(otu_mat) <- trimws(rownames(otu_mat))

cat("\nOTU matrix dimensions:\n")
print(dim(otu_mat))

tax_df_raw <- read.delim(
  tax_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(tax_df_raw) <- make.unique(trimws(colnames(tax_df_raw)))

tax_id_col <- colnames(tax_df_raw)[1]
tax_ids <- trimws(as.character(tax_df_raw[[tax_id_col]]))
tax_ids_unique <- make.unique(tax_ids)

tax_df <- tax_df_raw
rownames(tax_df) <- tax_ids_unique
tax_df[[tax_id_col]] <- tax_ids_unique

cat("\nTaxonomy table dimensions:\n")
print(dim(tax_df))

cat("\nNumber of duplicate raw taxonomy identifiers detected:\n")
print(sum(duplicated(tax_ids)))

otu_taxa_raw <- rownames(otu_mat)
otu_taxa_unique <- make.unique(otu_taxa_raw)
rownames(otu_mat) <- otu_taxa_unique

common_taxa <- intersect(rownames(otu_mat), rownames(tax_df))

if (length(common_taxa) == 0) {
  stop(
    "No matching taxa between OTU table and taxonomy table after make.unique(). ",
    "Check whether the first taxonomy column matches OTU rownames."
  )
}

otu_mat <- otu_mat[common_taxa, , drop = FALSE]
tax_df <- tax_df[common_taxa, , drop = FALSE]

cat("\nCommon taxa between OTU and taxonomy:\n")
print(length(common_taxa))

species_col <- find_first_existing_col(
  tax_df,
  c("Species", "species", "s", "Taxon", "taxon")
)

if (is.null(species_col)) {
  warning("No Species column found in taxonomy table. Using matched taxon IDs.")
  species_labels <- rownames(otu_mat)
} else {
  species_labels <- as.character(tax_df[[species_col]])
  species_labels <- trimws(species_labels)
  species_labels[is.na(species_labels) | species_labels == ""] <- rownames(otu_mat)[
    is.na(species_labels) | species_labels == ""
  ]
}

species_df <- as.data.frame(otu_mat) %>%
  rownames_to_column("Original_Taxon") %>%
  mutate(Species = species_labels) %>%
  select(-Original_Taxon) %>%
  group_by(Species) %>%
  summarise(
    across(everything(), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

species_mat <- species_df %>%
  column_to_rownames("Species") %>%
  as.matrix()

storage.mode(species_mat) <- "numeric"
species_mat[is.na(species_mat)] <- 0

sample_totals <- colSums(species_mat, na.rm = TRUE)
sample_totals[sample_totals == 0] <- NA

species_rel <- sweep(species_mat, 2, sample_totals, FUN = "/")
species_rel[is.na(species_rel)] <- 0

cat("\nSpecies relative-abundance matrix dimensions:\n")
print(dim(species_rel))

# ============================================================
# 7. Align samples across subsystem, reaction, species, metadata
# ============================================================

common_samples <- Reduce(
  intersect,
  list(
    colnames(subsystem_mat),
    colnames(reaction_mat),
    colnames(species_rel),
    metadata$Sample_Name
  )
)

cat("\nCommon samples across subsystem, reaction, species, metadata:\n")
print(length(common_samples))

if (length(common_samples) < 10) {
  stop("Too few common samples across subsystem, reaction, species, and metadata.")
}

metadata_aligned <- metadata %>%
  filter(Sample_Name %in% common_samples) %>%
  distinct(Sample_Name, .keep_all = TRUE) %>%
  arrange(match(Sample_Name, common_samples))

rownames(metadata_aligned) <- metadata_aligned$Sample_Name

subsystem_mat <- subsystem_mat[, metadata_aligned$Sample_Name, drop = FALSE]
reaction_mat  <- reaction_mat[, metadata_aligned$Sample_Name, drop = FALSE]
species_rel   <- species_rel[, metadata_aligned$Sample_Name, drop = FALSE]

exposed_samples <- metadata_aligned %>%
  filter(Condition == condition_of_interest) %>%
  pull(Sample_Name)

subsystem_exp <- subsystem_mat[, exposed_samples, drop = FALSE]
reaction_exp  <- reaction_mat[, exposed_samples, drop = FALSE]
species_exp   <- species_rel[, exposed_samples, drop = FALSE]

metadata_exp <- metadata_aligned %>%
  filter(Sample_Name %in% exposed_samples) %>%
  arrange(match(Sample_Name, exposed_samples))

rownames(metadata_exp) <- metadata_exp$Sample_Name

cat("\nNumber of PFAS-exposed samples used:\n")
print(length(exposed_samples))

# ============================================================
# 8. Filter variable features and select top subsystem rows
# ============================================================

subsystem_exp <- subsystem_exp[
  apply(subsystem_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

reaction_exp <- reaction_exp[
  apply(reaction_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

species_exp <- species_exp[
  apply(species_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

top_subsystems <- subsystem_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  mutate(total_abundance = rowSums(across(-Feature), na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = min(top_n_subsystems, nrow(.))) %>%
  pull(Feature)

subsystem_top <- subsystem_exp[top_subsystems, , drop = FALSE]
subsystem_top_z <- zscore_rows(subsystem_top)

# ============================================================
# 9. Derive two clusters from exposed subsystem heatmap
# ============================================================

sample_dist_subsystem <- dist(t(subsystem_top_z), method = "euclidean")
sample_hclust_subsystem <- hclust(sample_dist_subsystem, method = "ward.D2")

cluster_assignments <- cutree(sample_hclust_subsystem, k = k_clusters)

cluster_df <- tibble(
  Sample_Name = names(cluster_assignments),
  Cluster = paste0("Cluster_", cluster_assignments)
) %>%
  left_join(
    metadata_exp,
    by = "Sample_Name"
  )

write.csv(
  cluster_df,
  file.path(output_dir, "Exposed_SubsystemCluster_Assignments.csv"),
  row.names = FALSE
)

cat("\nCluster sizes:\n")
print(table(cluster_df$Cluster))

# ============================================================
# 10. Heatmap with subsystem-derived clusters
# ============================================================

annotation_col <- cluster_df %>%
  select(
    Sample_Name,
    Cluster,
    any_of(c("Average_k", "average_k", "Sex", "BMI", "Age"))
  ) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_Name
annotation_col$Sample_Name <- NULL

pheatmap(
  subsystem_top_z,
  cluster_rows = TRUE,
  cluster_cols = sample_hclust_subsystem,
  annotation_col = annotation_col,
  fontsize_row = 7,
  fontsize_col = 6,
  main = "PFAS-exposed samples: subsystem-derived clusters",
  filename = file.path(output_dir, "Diagnostic_Exposed_SubsystemHeatmap_WithClusters.png"),
  width = 12,
  height = 8
)

# ============================================================
# 11. Test cluster association with metadata variables
# ============================================================

available_technical_vars <- intersect(candidate_technical_vars, colnames(cluster_df))
available_biological_vars <- intersect(candidate_biological_vars, colnames(cluster_df))

all_candidate_vars <- unique(c(available_technical_vars, available_biological_vars))

metadata_tests <- purrr::map_dfr(
  all_candidate_vars,
  function(v) {
    x <- cluster_df[[v]]
    x_num <- to_numeric_safe(x)
    is_numeric_like <- sum(!is.na(x_num)) >= max(5, floor(0.5 * length(x)))
    
    if (is_numeric_like) {
      cluster_df[[v]] <- x_num
      safe_wilcox(cluster_df, v, group_col = "Cluster")
    } else {
      safe_fisher_or_chisq(cluster_df, v, group_col = "Cluster")
    }
  }
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    variable_type = case_when(
      variable %in% available_technical_vars ~ "technical",
      variable %in% available_biological_vars ~ "biological",
      TRUE ~ "other"
    )
  ) %>%
  arrange(p_value)

write.csv(
  metadata_tests,
  file.path(output_dir, "Cluster_Association_MetadataTests.csv"),
  row.names = FALSE
)

cat("\nTop metadata associations with subsystem-derived clusters:\n")
print(metadata_tests %>% slice_head(n = 20))

# ============================================================
# 12. Plot top numeric metadata associations
# ============================================================

top_numeric_vars <- metadata_tests %>%
  filter(test == "wilcox", !is.na(p_value)) %>%
  slice_head(n = 8) %>%
  pull(variable)

if (length(top_numeric_vars) > 0) {
  
  numeric_plot_df <- cluster_df %>%
    select(Sample_Name, Cluster, all_of(top_numeric_vars)) %>%
    pivot_longer(
      cols = all_of(top_numeric_vars),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Value = to_numeric_safe(Value)) %>%
    filter(!is.na(Value))
  
  p_numeric <- ggplot(
    numeric_plot_df,
    aes(x = Cluster, y = Value, fill = Cluster)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Value")
  
  ggsave(
    file.path(output_dir, "Cluster_Metadata_NumericVariables.png"),
    p_numeric,
    width = 12,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(output_dir, "Cluster_Metadata_NumericVariables.pdf"),
    p_numeric,
    width = 12,
    height = 7,
    bg = "white"
  )
}

# ============================================================
# 13. Identify subsystem features driving clusters
# ============================================================

subsystem_cluster_long <- subsystem_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

subsystem_cluster_tests <- subsystem_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    median_cluster_1 = median(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    median_cluster_2 = median(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1)
  ) %>%
  arrange(p_value)

write.csv(
  subsystem_cluster_tests,
  file.path(output_dir, "Subsystem_Features_Driving_ExposedClusters.csv"),
  row.names = FALSE
)

cat("\nTop subsystem features driving exposed clusters:\n")
print(subsystem_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 14. Identify reaction features driving clusters
# ============================================================

reaction_cluster_long <- reaction_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

reaction_cluster_tests <- reaction_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_cluster_1 = mean(Value[Cluster == "Cluster_1"] > presence_threshold, na.rm = TRUE),
    prevalence_cluster_2 = mean(Value[Cluster == "Cluster_2"] > presence_threshold, na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1),
    prevalence_difference = prevalence_cluster_2 - prevalence_cluster_1,
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(p_value)

write.csv(
  reaction_cluster_tests,
  file.path(output_dir, "Reaction_Features_Driving_ExposedClusters.csv"),
  row.names = FALSE
)

cat("\nTop reaction features driving exposed clusters:\n")
print(reaction_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 15. Identify species associated with subsystem-derived clusters
# ============================================================

species_cluster_long <- species_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

species_cluster_tests <- species_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    prevalence_cluster_1 = mean(Value[Cluster == "Cluster_1"] > 0, na.rm = TRUE),
    prevalence_cluster_2 = mean(Value[Cluster == "Cluster_2"] > 0, na.rm = TRUE),
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1),
    prevalence_difference = prevalence_cluster_2 - prevalence_cluster_1,
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(p_value)

write.csv(
  species_cluster_tests,
  file.path(output_dir, "Species_Features_Associated_With_SubsystemClusters.csv"),
  row.names = FALSE
)

cat("\nTop species associated with subsystem-derived clusters:\n")
print(species_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 16. Compare subsystem and species sample structure
# ============================================================

top_species <- species_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  mutate(row_sd = apply(species_exp, 1, sd, na.rm = TRUE)) %>%
  arrange(desc(row_sd)) %>%
  slice_head(n = min(top_n_species, nrow(.))) %>%
  pull(Feature)

species_top <- species_exp[top_species, , drop = FALSE]
species_top_z <- zscore_rows(species_top)

sample_dist_species <- dist(t(species_top_z), method = "euclidean")

mantel_res <- vegan::mantel(
  as.dist(sample_dist_subsystem),
  as.dist(sample_dist_species),
  method = "spearman",
  permutations = 999
)

capture.output(
  mantel_res,
  file = file.path(output_dir, "Mantel_Subsystem_vs_Species_Distances.txt")
)

cat("\nMantel test: subsystem vs species sample distances\n")
print(mantel_res)

# ============================================================
# 17. PCA comparison: subsystem vs species
# ============================================================

subsystem_pca <- prcomp(t(subsystem_top_z), scale. = FALSE)
species_pca <- prcomp(t(species_top_z), scale. = FALSE)

subsystem_pca_df <- as.data.frame(subsystem_pca$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

species_pca_df <- as.data.frame(species_pca$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

subsystem_var <- summary(subsystem_pca)$importance[2, 1:2] * 100
species_var <- summary(species_pca)$importance[2, 1:2] * 100

p_subsystem_pca <- ggplot(
  subsystem_pca_df,
  aes(x = PC1, y = PC2, color = Cluster, label = Sample_Name)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Subsystem-level PCA",
    x = paste0("PC1 (", round(subsystem_var[1], 1), "%)"),
    y = paste0("PC2 (", round(subsystem_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

p_species_pca <- ggplot(
  species_pca_df,
  aes(x = PC1, y = PC2, color = Cluster, label = Sample_Name)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Species-level PCA colored by subsystem cluster",
    x = paste0("PC1 (", round(species_var[1], 1), "%)"),
    y = paste0("PC2 (", round(species_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

p_pca_combined <- p_subsystem_pca | p_species_pca

ggsave(
  file.path(output_dir, "Subsystem_vs_Species_PCA_ClusterComparison.png"),
  p_pca_combined,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Subsystem_vs_Species_PCA_ClusterComparison.pdf"),
  p_pca_combined,
  width = 14,
  height = 6,
  bg = "white"
)

# ============================================================
# 18. Presence/absence diagnostics for possible stochastic features
# ============================================================

subsystem_presence_tests <- subsystem_cluster_long %>%
  mutate(Present = Value > presence_threshold) %>%
  group_by(Feature) %>%
  summarise(
    present_cluster_1 = mean(Present[Cluster == "Cluster_1"], na.rm = TRUE),
    present_cluster_2 = mean(Present[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_difference = present_cluster_2 - present_cluster_1,
    p_value = tryCatch(
      fisher.test(table(Cluster, Present))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(desc(abs_prevalence_difference), p_value)

write.csv(
  subsystem_presence_tests,
  file.path(output_dir, "Subsystem_PresenceAbsence_Driving_Clusters.csv"),
  row.names = FALSE
)

reaction_presence_tests <- reaction_cluster_long %>%
  mutate(Present = Value > presence_threshold) %>%
  group_by(Feature) %>%
  summarise(
    present_cluster_1 = mean(Present[Cluster == "Cluster_1"], na.rm = TRUE),
    present_cluster_2 = mean(Present[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_difference = present_cluster_2 - present_cluster_1,
    p_value = tryCatch(
      fisher.test(table(Cluster, Present))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(desc(abs_prevalence_difference), p_value)

write.csv(
  reaction_presence_tests,
  file.path(output_dir, "Reaction_PresenceAbsence_Driving_Clusters.csv"),
  row.names = FALSE
)

cat("\nTop subsystem presence/absence differences:\n")
print(subsystem_presence_tests %>% slice_head(n = 25))

cat("\nTop reaction presence/absence differences:\n")
print(reaction_presence_tests %>% slice_head(n = 25))

# ============================================================
# 19. Compact text report
# ============================================================

sink(file.path(output_dir, "Cluster_Diagnostic_Report.txt"))

cat("Diagnostic report for PFAS-exposed subsystem clustering\n")
cat("=====================================================\n\n")

cat("Paths used:\n")
cat("Subsystem: ", subsystem_path, "\n")
cat("Reaction:  ", reaction_path, "\n")
cat("Metadata:  ", metadata_path, "\n")
cat("OTU:       ", otu_file, "\n")
cat("Taxonomy:  ", tax_file, "\n\n")

cat("Number of exposed samples:\n")
print(length(exposed_samples))

cat("\nNumber of duplicate raw taxonomy identifiers detected:\n")
print(sum(duplicated(tax_ids)))

cat("\nNumber of common taxa after make.unique alignment:\n")
print(length(common_taxa))

cat("\nCluster sizes:\n")
print(table(cluster_df$Cluster))

cat("\nTop metadata associations:\n")
print(metadata_tests %>% slice_head(n = 20))

cat("\nTop subsystem feature differences:\n")
print(subsystem_cluster_tests %>% slice_head(n = 25))

cat("\nTop reaction feature differences:\n")
print(reaction_cluster_tests %>% slice_head(n = 25))

cat("\nTop subsystem presence/absence differences:\n")
print(subsystem_presence_tests %>% slice_head(n = 25))

cat("\nTop reaction presence/absence differences:\n")
print(reaction_presence_tests %>% slice_head(n = 25))

cat("\nTop species differences using subsystem-derived cluster labels:\n")
print(species_cluster_tests %>% slice_head(n = 25))

cat("\nMantel test: subsystem vs species sample distances:\n")
print(mantel_res)

sink()

cat("\nDone.\n")
cat("Diagnostic outputs saved to:\n", output_dir, "\n")
