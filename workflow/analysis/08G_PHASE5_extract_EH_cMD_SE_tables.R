options(stringsAsFactors = FALSE)

rev_dir <- "/path/to/PFAS_mSystems_revision"
ph5 <- file.path(rev_dir, "04_public_dataset_harmonization")

out_raw <- file.path(ph5, "01_curatedMetagenomicData")
out_meta <- file.path(ph5, "02_metadata_filtered")
out_tax <- file.path(ph5, "04_harmonized_taxonomy")
out_path <- file.path(ph5, "05_humann_pathways")
out_ctrl <- file.path(ph5, "00_control")
cache_dir <- file.path(out_raw, "ExperimentHubCache")

dir.create(out_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out_meta, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tax, recursive = TRUE, showWarnings = FALSE)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ctrl, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(EXPERIMENT_HUB_CACHE = cache_dir)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_matrix_tsv <- function(mat, path) {
  df <- data.frame(feature_id = rownames(mat), as.data.frame(mat, check.names = FALSE), check.names = FALSE)
  write_tsv(df, path)
}

choose_col <- function(df, patterns) {
  nms <- names(df)
  hits <- unique(unlist(lapply(patterns, function(p) grep(p, nms, ignore.case = TRUE, value = TRUE))))
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

safe_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

filter_metadata <- function(md, cohort) {
  md$sample_id <- rownames(md)

  body_col <- choose_col(md, c("^body_site$", "body.*site", "bodysite", "sample_body_site"))
  disease_col <- choose_col(md, c("^disease$", "diagnosis", "condition", "study_condition", "health_status"))
  age_col <- choose_col(md, c("^age$", "age_years", "age_category"))
  subject_col <- choose_col(md, c("^subject_id$", "subject", "participant", "host", "person"))

  keep <- rep(TRUE, nrow(md))

  if (!is.na(body_col)) {
    body <- safe_string(md[[body_col]])
    keep <- keep & grepl("stool|feces|faeces|fecal|faecal", body, ignore.case = TRUE)
  }

  if (!is.na(disease_col)) {
    dz <- safe_string(md[[disease_col]])
    keep <- keep & (dz == "" | grepl("healthy|control|none|no disease|^n$", dz, ignore.case = TRUE))
  }

  if (!is.na(age_col)) {
    age_num <- suppressWarnings(as.numeric(md[[age_col]]))
    keep <- keep & (is.na(age_num) | age_num >= 18)
  }

  md_f <- md[keep, , drop = FALSE]

  if (!is.na(subject_col) && nrow(md_f) > 0) {
    sid <- safe_string(md_f[[subject_col]])
    missing <- sid == ""
    sid[missing] <- paste0("missing_subject_", seq_len(sum(missing)))
    md_f$.subject_for_dedup <- sid
    md_f <- md_f[order(md_f$.subject_for_dedup, md_f$sample_id), , drop = FALSE]
    md_f <- md_f[!duplicated(md_f$.subject_for_dedup), , drop = FALSE]
    md_f$.subject_for_dedup <- NULL
  }

  md_f$cohort_for_revision <- cohort
  md_f
}

extract_assay_matrix <- function(obj) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment not available")
  }

  an <- SummarizedExperiment::assayNames(obj)
  if (length(an) > 0) {
    mat <- SummarizedExperiment::assay(obj, an[1])
  } else {
    mat <- SummarizedExperiment::assay(obj)
  }

  mat <- as.matrix(mat)
  mat
}

extract_metadata <- function(obj) {
  md <- as.data.frame(SummarizedExperiment::colData(obj))
  md[] <- lapply(md, function(z) {
    if (is.list(z)) vapply(z, function(y) paste(as.character(y), collapse = ";"), character(1)) else as.character(z)
  })
  md
}

make_genus_matrix <- function(mat) {
  feats <- rownames(mat)
  genus_keep <- grepl("g__", feats) & !grepl("s__", feats)
  gm <- mat[genus_keep, , drop = FALSE]

  genus <- rownames(gm)
  genus <- sub("^.*g__", "g__", genus)
  genus <- sub("\\|.*$", "", genus)

  rownames(gm) <- genus

  if (nrow(gm) > 0) {
    split_idx <- split(seq_len(nrow(gm)), rownames(gm))
    collapsed <- do.call(rbind, lapply(split_idx, function(i) colSums(gm[i, , drop = FALSE], na.rm = TRUE)))
    rownames(collapsed) <- names(split_idx)
    gm <- collapsed
  }

  gm
}

resources <- data.frame(
  cohort = c("NielsenHB_2014", "NielsenHB_2014", "NielsenHB_2014", "HMP_2012", "HMP_2012", "HMP_2012"),
  feature_type = c("relative_abundance", "marker_abundance", "pathway_abundance", "relative_abundance", "marker_abundance", "pathway_abundance"),
  EH_ID = c("EH5722", "EH5718", "EH5720", "EH5584", "EH5580", "EH5582"),
  stringsAsFactors = FALSE
)

if (!requireNamespace("ExperimentHub", quietly = TRUE)) stop("ExperimentHub missing")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment missing")

suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(SummarizedExperiment))

eh <- ExperimentHub::ExperimentHub()

summary_rows <- list()
metadata_by_cohort <- list()
relative_by_cohort <- list()

for (i in seq_len(nrow(resources))) {
  cohort <- resources$cohort[i]
  ftype <- resources$feature_type[i]
  ehid <- resources$EH_ID[i]

  prefix <- paste(cohort, ftype, ehid, sep = "__")

  message("Downloading/extracting ", prefix)

  obj <- tryCatch(eh[[ehid]], error = function(e) e)

  if (inherits(obj, "error")) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      cohort = cohort, feature_type = ftype, EH_ID = ehid,
      status = "FAILED_DOWNLOAD", error = conditionMessage(obj),
      n_features = NA, n_samples = NA, stringsAsFactors = FALSE
    )
    next
  }

  saveRDS(obj, file.path(out_raw, paste0(prefix, ".rds")))

  mat <- tryCatch(extract_assay_matrix(obj), error = function(e) e)
  md <- tryCatch(extract_metadata(obj), error = function(e) e)

  if (inherits(mat, "error") || inherits(md, "error")) {
    err <- paste(
      if (inherits(mat, "error")) conditionMessage(mat) else "",
      if (inherits(md, "error")) conditionMessage(md) else ""
    )
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      cohort = cohort, feature_type = ftype, EH_ID = ehid,
      status = "FAILED_EXTRACT", error = err,
      n_features = NA, n_samples = NA, stringsAsFactors = FALSE
    )
    next
  }

  write_matrix_tsv(mat, file.path(out_raw, paste0(prefix, "__matrix.tsv")))
  md$sample_id <- rownames(md)
  write_tsv(md, file.path(out_raw, paste0(prefix, "__metadata.tsv")))

  if (ftype == "relative_abundance") {
    metadata_by_cohort[[cohort]] <- md
    relative_by_cohort[[cohort]] <- mat
  }

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    cohort = cohort, feature_type = ftype, EH_ID = ehid,
    status = "COMPLETE", error = "",
    n_features = nrow(mat), n_samples = ncol(mat), stringsAsFactors = FALSE
  )
}

for (cohort in names(metadata_by_cohort)) {
  md <- metadata_by_cohort[[cohort]]
  mat <- relative_by_cohort[[cohort]]

  md_f <- filter_metadata(md, cohort)
  out_md <- if (cohort == "NielsenHB_2014") {
    file.path(out_meta, "Nielsen_metadata_filtered.tsv")
  } else {
    file.path(out_meta, "HMP_metadata_filtered.tsv")
  }

  write_tsv(md_f, out_md)

  sample_keep <- intersect(md_f$sample_id, colnames(mat))
  genus_mat <- make_genus_matrix(mat[, sample_keep, drop = FALSE])

  out_genus <- if (cohort == "NielsenHB_2014") {
    file.path(out_tax, "Nielsen_relative_abundance_genus_filtered.tsv")
  } else {
    file.path(out_tax, "HMP_relative_abundance_genus_filtered.tsv")
  }

  write_matrix_tsv(genus_mat, out_genus)
}

summary <- do.call(rbind, summary_rows)
write_tsv(summary, file.path(out_ctrl, "PHASE5_EH_cMD_extract_resource_summary.tsv"))

decision <- data.frame(
  Field = c(
    "status",
    "checked_at",
    "Nielsen_relative_abundance_EH",
    "HMP_relative_abundance_EH",
    "Nielsen_metadata_filtered",
    "HMP_metadata_filtered",
    "Nielsen_genus_matrix",
    "HMP_genus_matrix"
  ),
  Value = c(
    "COMPLETE",
    as.character(Sys.time()),
    "EH5722",
    "EH5584",
    file.path(out_meta, "Nielsen_metadata_filtered.tsv"),
    file.path(out_meta, "HMP_metadata_filtered.tsv"),
    file.path(out_tax, "Nielsen_relative_abundance_genus_filtered.tsv"),
    file.path(out_tax, "HMP_relative_abundance_genus_filtered.tsv")
  ),
  stringsAsFactors = FALSE
)

write_tsv(decision, file.path(out_ctrl, "PHASE5_EH_cMD_extract_summary.tsv"))

cat("Wrote:", file.path(out_ctrl, "PHASE5_EH_cMD_extract_summary.tsv"), "\n")
