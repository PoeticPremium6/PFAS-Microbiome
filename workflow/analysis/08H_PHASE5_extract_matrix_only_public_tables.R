options(stringsAsFactors = FALSE)

ph5 <- "/path/to/PFAS_mSystems_revision/04_public_dataset_harmonization"
raw_dir <- file.path(ph5, "01_curatedMetagenomicData")
meta_dir <- file.path(ph5, "02_metadata_filtered")
tax_dir <- file.path(ph5, "04_harmonized_taxonomy")
path_dir <- file.path(ph5, "05_humann_pathways")
ctrl_dir <- file.path(ph5, "00_control")

dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tax_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ctrl_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_mat <- function(mat, path) {
  df <- data.frame(feature_id = rownames(mat), as.data.frame(mat, check.names = FALSE), check.names = FALSE)
  write_tsv(df, path)
}

read_mat <- function(path) {
  x <- readRDS(path)
  if (!is.matrix(x)) x <- as.matrix(x)
  x
}

make_metadata_stub <- function(mat, cohort) {
  data.frame(
    sample_id = colnames(mat),
    cohort = cohort,
    body_site = "stool_or_public_dataset_default",
    adult_status = "assumed_or_not_filtered",
    health_status = "not_filtered_from_matrix_only_resource",
    inclusion_phase5_context = TRUE,
    metadata_status = "matrix_only_ExperimentHub_resource_no_colData",
    stringsAsFactors = FALSE
  )
}

make_genus_matrix <- function(mat) {
  feats <- rownames(mat)
  keep <- grepl("g__", feats) & !grepl("s__", feats)
  gm <- mat[keep, , drop = FALSE]

  genus <- rownames(gm)
  genus <- sub("^.*g__", "g__", genus)
  genus <- sub("\\|.*$", "", genus)
  rownames(gm) <- genus

  if (nrow(gm) > 0) {
    idx <- split(seq_len(nrow(gm)), rownames(gm))
    gm <- do.call(rbind, lapply(idx, function(i) colSums(gm[i, , drop = FALSE], na.rm = TRUE)))
    rownames(gm) <- names(idx)
  }

  gm
}

resources <- list(
  Nielsen_relative = file.path(raw_dir, "NielsenHB_2014__relative_abundance__EH5722.rds"),
  Nielsen_pathways = file.path(raw_dir, "NielsenHB_2014__pathway_abundance__EH5720.rds"),
  HMP_relative = file.path(raw_dir, "HMP_2012__relative_abundance__EH5584.rds"),
  HMP_pathways = file.path(raw_dir, "HMP_2012__pathway_abundance__EH5582.rds")
)

missing <- names(resources)[!file.exists(unlist(resources))]
if (length(missing) > 0) {
  stop(paste("Missing RDS resources:", paste(missing, collapse = ", ")))
}

nielsen_rel <- read_mat(resources$Nielsen_relative)
hmp_rel <- read_mat(resources$HMP_relative)
nielsen_path <- read_mat(resources$Nielsen_pathways)
hmp_path <- read_mat(resources$HMP_pathways)

nielsen_genus <- make_genus_matrix(nielsen_rel)
hmp_genus <- make_genus_matrix(hmp_rel)

write_tsv(make_metadata_stub(nielsen_rel, "NielsenHB_2014"), file.path(meta_dir, "Nielsen_metadata_filtered.tsv"))
write_tsv(make_metadata_stub(hmp_rel, "HMP_2012"), file.path(meta_dir, "HMP_metadata_filtered.tsv"))

write_mat(nielsen_rel, file.path(tax_dir, "Nielsen_relative_abundance_all_taxa_matrix_only.tsv"))
write_mat(hmp_rel, file.path(tax_dir, "HMP_relative_abundance_all_taxa_matrix_only.tsv"))
write_mat(nielsen_genus, file.path(tax_dir, "Nielsen_relative_abundance_genus_filtered.tsv"))
write_mat(hmp_genus, file.path(tax_dir, "HMP_relative_abundance_genus_filtered.tsv"))

write_mat(nielsen_path, file.path(path_dir, "Nielsen_pathway_abundance_matrix_only.tsv"))
write_mat(hmp_path, file.path(path_dir, "HMP_pathway_abundance_matrix_only.tsv"))

summary <- data.frame(
  Field = c(
    "status",
    "checked_at",
    "decision",
    "Nielsen_relative_matrix_dim",
    "HMP_relative_matrix_dim",
    "Nielsen_genus_matrix_dim",
    "HMP_genus_matrix_dim",
    "Nielsen_pathway_matrix_dim",
    "HMP_pathway_matrix_dim",
    "metadata_limitation",
    "raw_reprocessing_decision"
  ),
  Value = c(
    "COMPLETE",
    as.character(Sys.time()),
    "Use matrix-only ExperimentHub cMD resources for public contextual comparison unless metadata limitation becomes unacceptable.",
    paste(dim(nielsen_rel), collapse = "x"),
    paste(dim(hmp_rel), collapse = "x"),
    paste(dim(nielsen_genus), collapse = "x"),
    paste(dim(hmp_genus), collapse = "x"),
    paste(dim(nielsen_path), collapse = "x"),
    paste(dim(hmp_path), collapse = "x"),
    "ExperimentHub resources loaded as matrices without colData; sample-level metadata filtering not available from these objects. Public data used as broad contextual adult-stool reference, not direct PFAS-unexposed controls.",
    "Do not switch to raw FASTQs unless reviewers require stricter public metadata matching or unless sample metadata cannot be recovered elsewhere."
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary, file.path(ctrl_dir, "PHASE5_matrix_only_public_extract_summary.tsv"))

cat("Wrote matrix-only public extract summary\n")
