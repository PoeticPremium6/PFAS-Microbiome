options(stringsAsFactors = FALSE)

rev_dir <- "/path/to/PFAS_mSystems_revision"
ph5 <- file.path(rev_dir, "04_public_dataset_harmonization")
ph4 <- file.path(rev_dir, "03_taxonomy_metaphlan")

ctrl_dir <- file.path(ph5, "00_control")
meta_dir <- file.path(ph5, "02_metadata_filtered")
tax_dir <- file.path(ph5, "04_harmonized_taxonomy")
path_dir <- file.path(ph5, "05_humann_pathways")

dir.create(ctrl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tax_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_matrix_tsv <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"feature_id" %in% names(x)) {
    stop(paste("Missing feature_id column:", path))
  }
  feats <- x$feature_id
  x$feature_id <- NULL
  mat <- as.matrix(x)
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- feats
  mat[is.na(mat)] <- 0
  mat
}

write_matrix_tsv <- function(mat, path) {
  df <- data.frame(feature_id = rownames(mat), as.data.frame(mat, check.names = FALSE), check.names = FALSE)
  write_tsv(df, path)
}

prefix_cols <- function(mat, prefix) {
  colnames(mat) <- paste(prefix, colnames(mat), sep = "__")
  mat
}

collapse_duplicate_features <- function(mat) {
  if (nrow(mat) == 0) return(mat)
  idx <- split(seq_len(nrow(mat)), rownames(mat))
  out <- do.call(rbind, lapply(idx, function(i) colSums(mat[i, , drop = FALSE], na.rm = TRUE)))
  rownames(out) <- names(idx)
  out
}

harmonize_union <- function(mats) {
  all_features <- sort(unique(unlist(lapply(mats, rownames))))
  out <- lapply(mats, function(m) {
    missing <- setdiff(all_features, rownames(m))
    if (length(missing) > 0) {
      pad <- matrix(0, nrow = length(missing), ncol = ncol(m))
      rownames(pad) <- missing
      colnames(pad) <- colnames(m)
      m <- rbind(m, pad)
    }
    m[all_features, , drop = FALSE]
  })
  do.call(cbind, out)
}

harmonize_intersection <- function(mats) {
  shared <- Reduce(intersect, lapply(mats, rownames))
  shared <- sort(shared)
  out <- lapply(mats, function(m) m[shared, , drop = FALSE])
  do.call(cbind, out)
}

matrix_summary <- function(mat, cohort) {
  data.frame(
    cohort = cohort,
    n_features = nrow(mat),
    n_samples = ncol(mat),
    zero_features = sum(rowSums(mat, na.rm = TRUE) == 0),
    zero_samples = sum(colSums(mat, na.rm = TRUE) == 0),
    median_sample_sum = median(colSums(mat, na.rm = TRUE)),
    min_sample_sum = min(colSums(mat, na.rm = TRUE)),
    max_sample_sum = max(colSums(mat, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

make_public_sample_metadata <- function() {
  files <- c(
    NielsenHB_2014 = file.path(meta_dir, "Nielsen_metadata_filtered.tsv"),
    HMP_2012 = file.path(meta_dir, "HMP_metadata_filtered.tsv")
  )

  out <- list()
  for (cohort in names(files)) {
    if (!file.exists(files[[cohort]])) next
    md <- read.delim(files[[cohort]], check.names = FALSE, stringsAsFactors = FALSE)
    if (!"sample_id" %in% names(md)) next
    md$sample_id_original <- md$sample_id
    md$sample_id <- paste(cohort, md$sample_id_original, sep = "__")
    md$phase5_public_status <- "matrix_only_public_context"
    out[[cohort]] <- md
  }

  if (length(out) == 0) return(data.frame())
  cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- ""
    x[, cols, drop = FALSE]
  })
  do.call(rbind, out)
}

nielsen_genus_path <- file.path(tax_dir, "Nielsen_relative_abundance_genus_filtered.tsv")
hmp_genus_path <- file.path(tax_dir, "HMP_relative_abundance_genus_filtered.tsv")

if (!file.exists(nielsen_genus_path)) stop(paste("Missing:", nielsen_genus_path))
if (!file.exists(hmp_genus_path)) stop(paste("Missing:", hmp_genus_path))

nielsen_genus <- collapse_duplicate_features(read_matrix_tsv(nielsen_genus_path))
hmp_genus <- collapse_duplicate_features(read_matrix_tsv(hmp_genus_path))

nielsen_genus <- prefix_cols(nielsen_genus, "NielsenHB_2014")
hmp_genus <- prefix_cols(hmp_genus, "HMP_2012")

public_genus_union <- harmonize_union(list(NielsenHB_2014 = nielsen_genus, HMP_2012 = hmp_genus))
public_genus_shared <- harmonize_intersection(list(NielsenHB_2014 = nielsen_genus, HMP_2012 = hmp_genus))

write_matrix_tsv(public_genus_union, file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_union.tsv"))
write_matrix_tsv(public_genus_shared, file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_shared_features.tsv"))

sample_md <- make_public_sample_metadata()
if (nrow(sample_md) > 0) {
  write_tsv(sample_md, file.path(meta_dir, "Nielsen_HMP_public_sample_metadata_matrix_only.tsv"))
}

ronneby_candidates <- c(
  file.path(ph4, "03_merged", "Ronneby_MetaPhlAn_genus.tsv"),
  file.path(rev_dir, "12_tables", "Table_S6_MetaPhlAn_genus.tsv")
)
ronneby_genus_path <- ronneby_candidates[file.exists(ronneby_candidates)][1]

ronneby_status <- "WAITING_FOR_PHASE4_METAPHLAN_GENUS_TABLE"
ronneby_full_union_file <- ""
ronneby_full_shared_file <- ""

if (!is.na(ronneby_genus_path) && length(ronneby_genus_path) == 1) {
  ronneby_genus <- collapse_duplicate_features(read_matrix_tsv(ronneby_genus_path))
  ronneby_genus <- prefix_cols(ronneby_genus, "Ronneby")

  full_union <- harmonize_union(list(Ronneby = ronneby_genus, NielsenHB_2014 = nielsen_genus, HMP_2012 = hmp_genus))
  full_shared <- harmonize_intersection(list(Ronneby = ronneby_genus, NielsenHB_2014 = nielsen_genus, HMP_2012 = hmp_genus))

  ronneby_full_union_file <- file.path(tax_dir, "Ronneby_Nielsen_HMP_genus_harmonized_union.tsv")
  ronneby_full_shared_file <- file.path(tax_dir, "Ronneby_Nielsen_HMP_genus_harmonized_shared_features.tsv")

  write_matrix_tsv(full_union, ronneby_full_union_file)
  write_matrix_tsv(full_shared, ronneby_full_shared_file)

  ronneby_status <- "COMPLETE_FULL_HARMONIZATION"
} else {
  wait <- data.frame(
    Field = c("status", "needed_file", "public_union_ready", "public_shared_ready"),
    Value = c(
      "WAITING_FOR_PHASE4_METAPHLAN_GENUS_TABLE",
      paste(ronneby_candidates, collapse = ";"),
      file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_union.tsv"),
      file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_shared_features.tsv")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(wait, file.path(tax_dir, "Ronneby_Nielsen_HMP_genus_harmonization_WAITING_FOR_RONNEBY.tsv"))
}

nielsen_pathway_path <- file.path(path_dir, "Nielsen_pathway_abundance_matrix_only.tsv")
hmp_pathway_path <- file.path(path_dir, "HMP_pathway_abundance_matrix_only.tsv")

pathway_status <- "NOT_RUN_MISSING_PUBLIC_PATHWAY_INPUTS"
public_pathway_union_file <- ""
public_pathway_shared_file <- ""

if (file.exists(nielsen_pathway_path) && file.exists(hmp_pathway_path)) {
  nielsen_pathway <- collapse_duplicate_features(read_matrix_tsv(nielsen_pathway_path))
  hmp_pathway <- collapse_duplicate_features(read_matrix_tsv(hmp_pathway_path))

  nielsen_pathway <- prefix_cols(nielsen_pathway, "NielsenHB_2014")
  hmp_pathway <- prefix_cols(hmp_pathway, "HMP_2012")

  public_pathway_union <- harmonize_union(list(NielsenHB_2014 = nielsen_pathway, HMP_2012 = hmp_pathway))
  public_pathway_shared <- harmonize_intersection(list(NielsenHB_2014 = nielsen_pathway, HMP_2012 = hmp_pathway))

  public_pathway_union_file <- file.path(path_dir, "Nielsen_HMP_HUMAnN_pathways_public_harmonized_union.tsv")
  public_pathway_shared_file <- file.path(path_dir, "Nielsen_HMP_HUMAnN_pathways_public_harmonized_shared_features.tsv")

  write_matrix_tsv(public_pathway_union, public_pathway_union_file)
  write_matrix_tsv(public_pathway_shared, public_pathway_shared_file)

  pathway_status <- "PUBLIC_ONLY_COMPLETE_WAITING_FOR_RONNEBY_HUMANN_DECISION"
}

feature_summary <- rbind(
  matrix_summary(nielsen_genus, "NielsenHB_2014_genus"),
  matrix_summary(hmp_genus, "HMP_2012_genus"),
  matrix_summary(public_genus_union, "Nielsen_HMP_genus_union"),
  matrix_summary(public_genus_shared, "Nielsen_HMP_genus_shared")
)
write_tsv(feature_summary, file.path(tax_dir, "PHASE5_public_genus_harmonization_feature_summary.tsv"))

decision <- data.frame(
  Field = c(
    "status",
    "checked_at",
    "phase5_3_raw_public_fastq_decision",
    "phase5_3_reason",
    "phase5_4_public_genus_status",
    "phase5_4_full_ronneby_public_status",
    "phase5_5_public_pathway_status",
    "Nielsen_HMP_genus_union",
    "Nielsen_HMP_genus_shared_features",
    "Ronneby_Nielsen_HMP_genus_union",
    "Ronneby_Nielsen_HMP_genus_shared_features",
    "Nielsen_HMP_pathway_union",
    "Nielsen_HMP_pathway_shared_features",
    "sample_metadata",
    "metadata_limitation"
  ),
  Value = c(
    "COMPLETE_PUBLIC_PREHARMONIZATION",
    as.character(Sys.time()),
    "DEFER_RAW_FASTQ_REPROCESSING",
    "Public cMD/ExperimentHub matrices are usable as contextual references. Raw reprocessing remains fallback if reviewers require matched metadata or same-pipeline processing.",
    "COMPLETE_PUBLIC_ONLY",
    ronneby_status,
    pathway_status,
    file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_union.tsv"),
    file.path(tax_dir, "Nielsen_HMP_genus_public_harmonized_shared_features.tsv"),
    ronneby_full_union_file,
    ronneby_full_shared_file,
    public_pathway_union_file,
    public_pathway_shared_file,
    file.path(meta_dir, "Nielsen_HMP_public_sample_metadata_matrix_only.tsv"),
    "ExperimentHub matrix resources have no colData in this direct-loading route; metadata files are stub/sample inclusion records. Treat public cohorts as context, not matched controls."
  ),
  stringsAsFactors = FALSE
)

write_tsv(decision, file.path(ctrl_dir, "PHASE5_public_preharmonization_and_raw_decision_summary.tsv"))

log_txt <- file.path(ph5, "03_raw_reprocessing_decision", "public_processing_decision_log.txt")
dir.create(dirname(log_txt), recursive = TRUE, showWarnings = FALSE)
cat(
  "Phase 5 public data processing decision log\n\n",
  "Current decision:\n",
  "Use direct ExperimentHub matrix resources for NielsenHB_2014 and HMP_2012 as contextual public reference cohorts.\n\n",
  "Raw public FASTQ reprocessing decision:\n",
  "Deferred. Raw reprocessing is retained as a fallback, not the current main route.\n\n",
  "Reason:\n",
  "The cMD package itself could not be loaded reproducibly because of rbiom dependency conflicts, but the ExperimentHub matrix resources are directly usable. These matrices allow genus-level and pathway-level public-context harmonization while Ronneby MetaPhlAn continues. Since the direct resources lack sample-level colData in this route, the public cohorts should be described as contextual reference profiles, not matched unexposed controls.\n\n",
  "Trigger to switch to raw FASTQs:\n",
  "Switch only if reviewers require sample-level public metadata matching, exact same-pipeline public processing, or if the contextual-matrix limitation is judged unacceptable.\n\n",
  "Next dependency:\n",
  "Full Ronneby-Nielsen-HMP genus harmonization waits for the Phase 4 Ronneby MetaPhlAn genus table.\n",
  file = log_txt
)

cat("Wrote Phase 5 public preharmonization and raw decision summary\n")
