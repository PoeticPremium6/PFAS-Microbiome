options(stringsAsFactors = FALSE)

rev_dir <- "/path/to/PFAS_mSystems_revision"
ph5 <- file.path(rev_dir, "04_public_dataset_harmonization")

dir.create(file.path(ph5, "00_control"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ph5, "01_curatedMetagenomicData"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ph5, "02_metadata_filtered"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ph5, "03_raw_reprocessing_decision"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ph5, "06_logs"), recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(ph5, "00_control", "PHASE5_cMD_inventory_summary.tsv")

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
  out <- data.frame(
    Field = c("status", "reason", "needed_package"),
    Value = c("FAILED", "curatedMetagenomicData not available in current R library", "curatedMetagenomicData")
  )
  write_tsv(out, summary_file)
  stop("curatedMetagenomicData not available")
}

suppressPackageStartupMessages(library(curatedMetagenomicData))

data("sampleMetadata", package = "curatedMetagenomicData")
md <- as.data.frame(sampleMetadata)
md$cMD_row_id <- rownames(md)

cols <- names(md)

find_col <- function(patterns) {
  hits <- unique(unlist(lapply(patterns, function(p) grep(p, cols, ignore.case = TRUE, value = TRUE))))
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

study_col <- find_col(c("^study_name$", "^study$", "study"))
body_col <- find_col(c("^body_site$", "body.*site", "bodysite"))
disease_col <- find_col(c("^disease$", "diagnosis", "condition", "study_condition"))
age_col <- find_col(c("^age$", "age_years"))
subject_col <- find_col(c("^subject_id$", "subject", "participant", "host"))
sample_col <- find_col(c("^sample_id$", "sample"))

accession_cols <- unique(grep("accession|sra|ena|run|ncbi|read|fastq|ftp|url", cols, ignore.case = TRUE, value = TRUE))

metadata_cols <- data.frame(
  column_name = cols,
  is_accession_like = cols %in% accession_cols,
  stringsAsFactors = FALSE
)
write_tsv(metadata_cols, file.path(ph5, "01_curatedMetagenomicData", "cMD_sampleMetadata_columns.tsv"))

if (is.na(study_col)) stop("Could not identify study column in curatedMetagenomicData metadata")

study_counts <- as.data.frame(table(md[[study_col]], useNA = "ifany"))
names(study_counts) <- c("study_name", "n_samples")
study_counts <- study_counts[order(study_counts$study_name), ]
write_tsv(study_counts, file.path(ph5, "01_curatedMetagenomicData", "cMD_study_counts.tsv"))

nielsen <- md[md[[study_col]] == "NielsenHB_2014", , drop = FALSE]

hmp_study_names <- unique(md[[study_col]][grepl("^HMP|HMP_", md[[study_col]], ignore.case = TRUE)])
hmp <- md[md[[study_col]] %in% hmp_study_names, , drop = FALSE]

is_stool <- function(x) {
  if (is.na(body_col)) return(rep(TRUE, nrow(x)))
  grepl("stool|feces|faeces|fecal|faecal", x[[body_col]], ignore.case = TRUE)
}

is_healthy_or_control <- function(x) {
  if (is.na(disease_col)) return(rep(TRUE, nrow(x)))
  z <- x[[disease_col]]
  is.na(z) | z == "" | grepl("healthy|control|none|no disease|n$", z, ignore.case = TRUE)
}

is_adult_or_unknown <- function(x) {
  if (is.na(age_col)) return(rep(TRUE, nrow(x)))
  age_num <- suppressWarnings(as.numeric(x[[age_col]]))
  is.na(age_num) | age_num >= 18
}

nielsen_filt <- nielsen[is_stool(nielsen) & is_healthy_or_control(nielsen) & is_adult_or_unknown(nielsen), , drop = FALSE]
hmp_filt <- hmp[is_stool(hmp) & is_healthy_or_control(hmp) & is_adult_or_unknown(hmp), , drop = FALSE]

dedup_one_per_subject <- function(x) {
  if (nrow(x) == 0) return(x)
  if (is.na(subject_col)) return(x)
  sid <- x[[subject_col]]
  sid[is.na(sid) | sid == ""] <- paste0("missing_subject_", seq_len(sum(is.na(sid) | sid == "")))
  x$.subject_for_dedup <- sid
  x <- x[order(x$.subject_for_dedup, x$cMD_row_id), , drop = FALSE]
  x <- x[!duplicated(x$.subject_for_dedup), , drop = FALSE]
  x$.subject_for_dedup <- NULL
  x
}

hmp_filt_one_per_subject <- dedup_one_per_subject(hmp_filt)

write_tsv(nielsen, file.path(ph5, "01_curatedMetagenomicData", "NielsenHB_2014_metadata_all.tsv"))
write_tsv(hmp, file.path(ph5, "01_curatedMetagenomicData", "HMP_metadata_all.tsv"))
write_tsv(nielsen_filt, file.path(ph5, "02_metadata_filtered", "Nielsen_metadata_filtered.tsv"))
write_tsv(hmp_filt_one_per_subject, file.path(ph5, "02_metadata_filtered", "HMP_metadata_filtered.tsv"))

make_accession_long <- function(x, cohort) {
  if (nrow(x) == 0 || length(accession_cols) == 0) {
    return(data.frame())
  }

  sample_id <- if (!is.na(sample_col)) x[[sample_col]] else x$cMD_row_id
  subject_id <- if (!is.na(subject_col)) x[[subject_col]] else NA

  out <- list()
  k <- 1

  for (ac in accession_cols) {
    vals <- x[[ac]]
    keep <- !is.na(vals) & vals != ""
    if (any(keep)) {
      out[[k]] <- data.frame(
        cohort = cohort,
        cMD_row_id = x$cMD_row_id[keep],
        sample_id = sample_id[keep],
        subject_id = subject_id[keep],
        accession_column = ac,
        accession_value = vals[keep],
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

acc_nielsen <- make_accession_long(nielsen_filt, "NielsenHB_2014")
acc_hmp <- make_accession_long(hmp_filt_one_per_subject, "HMP")

acc_all <- rbind(acc_nielsen, acc_hmp)
write_tsv(acc_all, file.path(ph5, "03_raw_reprocessing_decision", "public_raw_accession_inventory.tsv"))

if (nrow(acc_all) > 0) {
  acc_summary <- aggregate(accession_value ~ cohort + accession_column, acc_all, function(z) length(unique(z)))
  names(acc_summary)[3] <- "unique_accession_values"
} else {
  acc_summary <- data.frame(cohort = character(), accession_column = character(), unique_accession_values = integer())
}
write_tsv(acc_summary, file.path(ph5, "03_raw_reprocessing_decision", "public_raw_accession_summary.tsv"))

decision_log <- file.path(ph5, "03_raw_reprocessing_decision", "public_processing_decision_log.txt")
cat(
  "Phase 5 public processing decision log\n\n",
  "Current decision:\n",
  "Use curatedMetagenomicData metadata/profiles for Phase 5.1 and 5.2 immediately.\n",
  "Assess raw public FASTQ reprocessing feasibility using accession inventory before launching downloads.\n\n",
  "Raw reprocessing rule:\n",
  "Only start public raw FASTQ reprocessing after Ronneby Phase 4 MetaPhlAn settings are stable and after Nielsen/HMP accession IDs are confirmed.\n\n",
  "Generated files:\n",
  "- 02_metadata_filtered/Nielsen_metadata_filtered.tsv\n",
  "- 02_metadata_filtered/HMP_metadata_filtered.tsv\n",
  "- 03_raw_reprocessing_decision/public_raw_accession_inventory.tsv\n",
  "- 03_raw_reprocessing_decision/public_raw_accession_summary.tsv\n",
  sep = "",
  file = decision_log
)

summary <- data.frame(
  Field = c(
    "status",
    "checked_at",
    "study_column",
    "body_site_column",
    "disease_column",
    "age_column",
    "subject_column",
    "sample_column",
    "accession_like_columns",
    "Nielsen_all_samples",
    "Nielsen_filtered_samples",
    "HMP_study_names",
    "HMP_all_samples",
    "HMP_filtered_one_per_subject_samples",
    "raw_accession_inventory_rows",
    "raw_accession_summary_rows",
    "Nielsen_metadata_filtered",
    "HMP_metadata_filtered",
    "raw_accession_inventory",
    "raw_accession_summary"
  ),
  Value = c(
    "COMPLETE",
    as.character(Sys.time()),
    study_col,
    body_col,
    disease_col,
    age_col,
    subject_col,
    sample_col,
    paste(accession_cols, collapse = ","),
    nrow(nielsen),
    nrow(nielsen_filt),
    paste(hmp_study_names, collapse = ","),
    nrow(hmp),
    nrow(hmp_filt_one_per_subject),
    nrow(acc_all),
    nrow(acc_summary),
    file.path(ph5, "02_metadata_filtered", "Nielsen_metadata_filtered.tsv"),
    file.path(ph5, "02_metadata_filtered", "HMP_metadata_filtered.tsv"),
    file.path(ph5, "03_raw_reprocessing_decision", "public_raw_accession_inventory.tsv"),
    file.path(ph5, "03_raw_reprocessing_decision", "public_raw_accession_summary.tsv")
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary, summary_file)
cat("Wrote summary:", summary_file, "\n")
