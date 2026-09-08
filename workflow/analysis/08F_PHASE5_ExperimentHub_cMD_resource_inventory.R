options(stringsAsFactors = FALSE)

rev_dir <- "/path/to/PFAS_mSystems_revision"
ph5 <- file.path(rev_dir, "04_public_dataset_harmonization")

out_dir <- file.path(ph5, "01_curatedMetagenomicData")
log_dir <- file.path(ph5, "06_logs")
ctrl_dir <- file.path(ph5, "00_control")
cache_dir <- file.path(out_dir, "ExperimentHubCache")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ctrl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(EXPERIMENT_HUB_CACHE = cache_dir)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

flatten_df <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(out)) {
    if (is.list(out[[nm]])) {
      out[[nm]] <- vapply(out[[nm]], function(z) paste(as.character(z), collapse = ";"), character(1))
    } else {
      out[[nm]] <- as.character(out[[nm]])
    }
  }
  out
}

summary_file <- file.path(ctrl_dir, "PHASE5_ExperimentHub_cMD_resource_inventory_summary.tsv")

if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
  write_tsv(
    data.frame(
      Field = c("status", "reason", "needed_package"),
      Value = c("FAILED", "ExperimentHub not available in current R environment", "ExperimentHub")
    ),
    summary_file
  )
  stop("ExperimentHub not available")
}

suppressPackageStartupMessages(library(ExperimentHub))

eh <- ExperimentHub::ExperimentHub()
eh_md <- flatten_df(S4Vectors::mcols(eh))
eh_md$EH_ID <- names(eh)

if (!"title" %in% names(eh_md)) eh_md$title <- ""
if (!"description" %in% names(eh_md)) eh_md$description <- ""

write_tsv(eh_md, file.path(out_dir, "ExperimentHub_all_resources_metadata.tsv"))

row_text <- apply(eh_md, 1, paste, collapse = " | ")

cmd_hits <- eh_md[grepl("curatedMetagenomicData|curated metagenomic", row_text, ignore.case = TRUE), , drop = FALSE]
nielsen_hits <- eh_md[grepl("NielsenHB_2014|Nielsen", row_text, ignore.case = TRUE), , drop = FALSE]
hmp_hits <- eh_md[grepl("HMP|Human Microbiome|Huttenhower|Lloyd-Price", row_text, ignore.case = TRUE), , drop = FALSE]

write_tsv(cmd_hits, file.path(out_dir, "ExperimentHub_curatedMetagenomicData_resource_hits.tsv"))
write_tsv(nielsen_hits, file.path(out_dir, "ExperimentHub_Nielsen_resource_hits.tsv"))
write_tsv(hmp_hits, file.path(out_dir, "ExperimentHub_HMP_resource_hits.tsv"))

query_to_df <- function(query_terms, label) {
  q <- tryCatch(ExperimentHub::query(eh, query_terms), error = function(e) e)
  if (inherits(q, "error")) {
    return(data.frame(query_label = label, error = conditionMessage(q), stringsAsFactors = FALSE))
  }

  if (length(q) == 0) {
    return(data.frame(query_label = label, error = "no_results", stringsAsFactors = FALSE))
  }

  md <- flatten_df(S4Vectors::mcols(q))
  md$EH_ID <- names(q)
  md$query_label <- label
  md
}

q1 <- query_to_df(c("curatedMetagenomicData", "NielsenHB_2014"), "cMD_NielsenHB_2014")
q2 <- query_to_df(c("curatedMetagenomicData", "HMP"), "cMD_HMP")
q3 <- query_to_df(c("curatedMetagenomicData", "relative_abundance"), "cMD_relative_abundance")
q4 <- query_to_df(c("curatedMetagenomicData", "marker_abundance"), "cMD_marker_abundance")

queries <- rbind(
  q1[, union(names(q1), character()), drop = FALSE],
  q2[, union(names(q2), character()), drop = FALSE],
  q3[, union(names(q3), character()), drop = FALSE],
  q4[, union(names(q4), character()), drop = FALSE]
)

all_names <- unique(unlist(list(names(q1), names(q2), names(q3), names(q4))))
pad <- function(x) {
  missing <- setdiff(all_names, names(x))
  for (m in missing) x[[m]] <- ""
  x[, all_names, drop = FALSE]
}

queries <- rbind(pad(q1), pad(q2), pad(q3), pad(q4))
write_tsv(queries, file.path(out_dir, "ExperimentHub_cMD_target_query_results.tsv"))

summary <- data.frame(
  Field = c(
    "status",
    "checked_at",
    "R_version",
    "ExperimentHub_available",
    "ExperimentHub_cache",
    "all_EH_resources",
    "curatedMetagenomicData_hit_rows",
    "Nielsen_hit_rows",
    "HMP_hit_rows",
    "target_query_rows",
    "all_resources_file",
    "cMD_hits_file",
    "Nielsen_hits_file",
    "HMP_hits_file",
    "target_query_results_file"
  ),
  Value = c(
    "COMPLETE",
    as.character(Sys.time()),
    R.version.string,
    "TRUE",
    cache_dir,
    nrow(eh_md),
    nrow(cmd_hits),
    nrow(nielsen_hits),
    nrow(hmp_hits),
    nrow(queries),
    file.path(out_dir, "ExperimentHub_all_resources_metadata.tsv"),
    file.path(out_dir, "ExperimentHub_curatedMetagenomicData_resource_hits.tsv"),
    file.path(out_dir, "ExperimentHub_Nielsen_resource_hits.tsv"),
    file.path(out_dir, "ExperimentHub_HMP_resource_hits.tsv"),
    file.path(out_dir, "ExperimentHub_cMD_target_query_results.tsv")
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary, summary_file)

cat("Wrote summary:", summary_file, "\n")
