rev <- "/path/to/PFAS_mSystems_revision"
ph6 <- file.path(rev, "05_humann_functional_profiles")
amat <- file.path(ph6, "04_analysis_matrices")
outdir <- file.path(ph6, "07_association_tests")
tabdir <- file.path(rev, "12_tables")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(tabdir, recursive=TRUE, showWarnings=FALSE)

meta_master_file <- file.path(rev, "01_metadata_freeze", "sample_metadata_revision_master.tsv")
pfas_file <- file.path(rev, "PFAS_Metadata.csv")
elim_file <- file.path(rev, "01_metadata_freeze", "elimination_rate_strata.tsv")

pathway_file <- file.path(amat, "Ronneby_HUMAnN_pathabundance_cpm_unstratified_no_unmapped_prevalence10.tsv")
genefam_file <- file.path(amat, "Ronneby_HUMAnN_genefamilies_cpm_top_variable_features.tsv")
genefam_matrix_file <- file.path(amat, "Ronneby_HUMAnN_genefamilies_cpm_unstratified_no_unmapped_prevalence10.tsv")

read_tsv <- function(x) read.table(x, header=TRUE, sep="\t", check.names=FALSE, quote="", comment.char="")
clean_sample <- function(x) {
  x <- as.character(x)
  vapply(x, function(z) {
    m <- regexpr("[A-H][0-9]{2}", z)
    if (m[1] > 0) {
      regmatches(z, m)
    } else {
      z
    }
  }, character(1))
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

cat("Reading metadata...\n")
meta <- read_tsv(meta_master_file)
pfas <- read.csv(pfas_file, check.names=FALSE)
elim <- read_tsv(elim_file)

if (!"SampleID" %in% names(meta)) stop("SampleID missing from sample_metadata_revision_master.tsv")
if (!"SampleID" %in% names(elim)) stop("SampleID missing from elimination_rate_strata.tsv")

pfas_sample_col <- intersect(c("SampleID", "Sample_Name", "sample_id", "sample", "Sample"), names(pfas))[1]
if (is.na(pfas_sample_col)) stop("No sample column detected in PFAS_Metadata.csv")
names(pfas)[names(pfas) == pfas_sample_col] <- "SampleID"

meta$SampleID <- as.character(meta$SampleID)
pfas$SampleID <- as.character(pfas$SampleID)
elim$SampleID <- as.character(elim$SampleID)

md <- merge(meta, pfas, by="SampleID", all.x=TRUE, suffixes=c("", "_PFAS"))
md <- merge(md, elim, by="SampleID", all.x=TRUE, suffixes=c("", "_ELIM"))

# Coerce common covariates.
for (cc in intersect(c("Age", "BMI"), names(md))) md[[cc]] <- to_num(md[[cc]])
for (cc in intersect(c("Sex", "ExposureGroup", "SequencingBatch"), names(md))) md[[cc]] <- as.factor(md[[cc]])

# Candidate phenotypes.
serum_pfas_candidates <- c(
  "PFOA","PFNA","PFDA","PFUnDA","PFDoA","PFTrDA",
  "PFBS","PFPeS","PFHxS","PFHpS","PFOS",
  "PFOS_branched_MP11","PFOS_branchedMP3_4_51","PFOS_branchedMP2_61",
  "PFNS","PFDS"
)
serum_pfas <- intersect(serum_pfas_candidates, names(md))

elim_numeric <- grep("^(k_|mean_k_all_available)", names(md), value=TRUE)
elim_numeric <- elim_numeric[sapply(md[elim_numeric], function(x) sum(!is.na(to_num(x))) >= 10)]

categorical <- intersect(c("ExposureGroup", "primary_elimination_stratum"), names(md))

phenotypes <- unique(c(categorical, serum_pfas, elim_numeric))

# Remove phenotypes with insufficient data.
ok_pheno <- c()
for (p in phenotypes) {
  if (p %in% categorical) {
    z <- md[[p]]
    z <- z[!is.na(z) & z != ""]
    if (length(unique(z)) >= 2 && length(z) >= 10) ok_pheno <- c(ok_pheno, p)
  } else {
    z <- to_num(md[[p]])
    if (sum(!is.na(z)) >= 10 && length(unique(z[!is.na(z)])) >= 4) ok_pheno <- c(ok_pheno, p)
  }
}
phenotypes <- ok_pheno

choose_covariates <- function(d, pheno) {
  covars <- c()
  if ("Age" %in% names(d) && pheno != "Age" && sum(!is.na(d$Age)) >= 20) covars <- c(covars, "Age")
  if ("Sex" %in% names(d) && pheno != "Sex" && length(unique(na.omit(d$Sex))) >= 2) covars <- c(covars, "Sex")
  if ("BMI" %in% names(d) && pheno != "BMI" && sum(!is.na(d$BMI)) >= 20) covars <- c(covars, "BMI")
  if ("SequencingBatch" %in% names(d) && pheno != "SequencingBatch" && length(unique(na.omit(d$SequencingBatch))) >= 2) covars <- c(covars, "SequencingBatch")
  covars
}

read_feature_matrix <- function(path, feature_col_name) {
  x <- read_tsv(path)
  feats <- x[[1]]
  mat <- as.matrix(x[, -1, drop=FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feats
  colnames(mat) <- clean_sample(colnames(mat))
  list(feature_col=feature_col_name, features=feats, mat=mat)
}

subset_genefam_top <- function(matrix_file, top_file, n_top=1000) {
  top <- read_tsv(top_file)
  top_features <- head(top[[1]], n_top)

  con <- file(matrix_file, "r")
  header <- strsplit(readLines(con, n=1), "\t", fixed=TRUE)[[1]]
  samples <- clean_sample(header[-1])

  vals <- list()
  feats <- c()
  repeat {
    line <- readLines(con, n=1)
    if (length(line) == 0) break
    parts <- strsplit(line, "\t", fixed=TRUE)[[1]]
    feat <- parts[1]
    if (feat %in% top_features) {
      nums <- suppressWarnings(as.numeric(parts[-1]))
      nums[is.na(nums)] <- 0
      vals[[feat]] <- nums
      feats <- c(feats, feat)
    }
  }
  close(con)

  mat <- do.call(rbind, vals)
  rownames(mat) <- names(vals)
  colnames(mat) <- samples
  list(feature_col="GeneFamily", features=rownames(mat), mat=mat)
}

run_tests <- function(obj, md, target_name) {
  mat <- obj$mat
  common <- intersect(colnames(mat), md$SampleID)
  mat <- mat[, common, drop=FALSE]
  d0 <- md[match(common, md$SampleID), , drop=FALSE]

  results <- list()
  idx <- 1

  for (pheno in phenotypes) {
    d <- d0

    if (!(pheno %in% names(d))) next

    if (pheno %in% categorical) {
      d[[pheno]] <- as.factor(d[[pheno]])
      keep <- !is.na(d[[pheno]]) & d[[pheno]] != ""
      if (sum(keep) < 10 || length(unique(d[[pheno]][keep])) < 2) next
      d[[pheno]] <- droplevels(d[[pheno]])
      pheno_for_model <- pheno
    } else {
      d[[pheno]] <- log10(to_num(d[[pheno]]) + 1)
      keep <- !is.na(d[[pheno]])
      if (sum(keep) < 10 || length(unique(d[[pheno]][keep])) < 4) next
      pheno_for_model <- pheno
    }

    covars <- choose_covariates(d, pheno)
    model_terms <- c(pheno_for_model, covars)
    form <- as.formula(paste("y ~", paste(model_terms, collapse=" + ")))

    pvals <- rep(NA_real_, nrow(mat))
    betas <- rep(NA_real_, nrow(mat))
    ns <- rep(NA_integer_, nrow(mat))

    for (i in seq_len(nrow(mat))) {
      y <- log1p(as.numeric(mat[i, ]))
      dd <- d
      dd$y <- y

      complete_vars <- c("y", model_terms)
      cc <- complete.cases(dd[, complete_vars, drop=FALSE])
      if (sum(cc) < 10) next

      fit <- tryCatch(lm(form, data=dd[cc, , drop=FALSE]), error=function(e) NULL)
      if (is.null(fit)) next

      a <- tryCatch(anova(fit), error=function(e) NULL)
      if (!is.null(a) && pheno_for_model %in% rownames(a)) {
        pvals[i] <- a[pheno_for_model, "Pr(>F)"]
      }

      co <- tryCatch(coef(summary(fit)), error=function(e) NULL)
      if (!is.null(co) && pheno_for_model %in% rownames(co)) {
        betas[i] <- co[pheno_for_model, "Estimate"]
      } else {
        betas[i] <- NA_real_
      }
      ns[i] <- sum(cc)
    }

    qvals <- p.adjust(pvals, method="BH")

    res <- data.frame(
      target=target_name,
      feature_type=obj$feature_col,
      feature=rownames(mat),
      phenotype=pheno,
      model=paste("log1p_feature_cpm ~", paste(model_terms, collapse=" + ")),
      n=ns,
      beta_or_group_effect=betas,
      p_value=pvals,
      q_value=qvals,
      interpretation=ifelse(!is.na(qvals) & qvals < 0.05, "FDR_significant",
                       ifelse(!is.na(pvals) & pvals < 0.05, "nominal_exploratory", "not_significant")),
      stringsAsFactors=FALSE
    )

    results[[idx]] <- res
    idx <- idx + 1
  }

  if (length(results) == 0) {
    return(data.frame())
  }
  do.call(rbind, results)
}

cat("Reading pathway matrix...\n")
path_obj <- read_feature_matrix(pathway_file, "Pathway")

cat("Reading top-variable gene-family subset...\n")
gf_obj <- subset_genefam_top(genefam_matrix_file, genefam_file, n_top=1000)

cat("Running pathway tests...\n")
path_res <- run_tests(path_obj, md, "HUMAnN_pathways_prevalence10")

cat("Running gene-family tests...\n")
gf_res <- run_tests(gf_obj, md, "HUMAnN_gene_families_top1000_variable")


result_cols <- c("target","feature_type","feature","phenotype","model","n","beta_or_group_effect","p_value","q_value","interpretation")

empty_results <- function() {
  x <- as.data.frame(setNames(replicate(length(result_cols), character(), simplify=FALSE), result_cols))
  x$n <- integer()
  x$beta_or_group_effect <- numeric()
  x$p_value <- numeric()
  x$q_value <- numeric()
  x
}

ensure_cols <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0 || !("q_value" %in% names(x)) || !("p_value" %in% names(x))) {
    return(empty_results())
  }
  missing <- setdiff(result_cols, names(x))
  for (m in missing) x[[m]] <- NA
  x[, result_cols, drop=FALSE]
}

path_res <- ensure_cols(path_res)
gf_res <- ensure_cols(gf_res)


all_res <- rbind(path_res, gf_res)
all_res <- all_res[order(all_res$q_value, all_res$p_value), ]

out_all <- file.path(outdir, "Ronneby_HUMAnN_pathway_gene_family_exposure_elimination_associations.tsv")
write.table(all_res, out_all, sep="\t", quote=FALSE, row.names=FALSE)

out_path <- file.path(outdir, "Ronneby_HUMAnN_pathway_exposure_elimination_associations.tsv")
write.table(path_res[order(path_res$q_value, path_res$p_value), ], out_path, sep="\t", quote=FALSE, row.names=FALSE)

out_gf <- file.path(outdir, "Ronneby_HUMAnN_gene_family_top1000_exposure_elimination_associations.tsv")
write.table(gf_res[order(gf_res$q_value, gf_res$p_value), ], out_gf, sep="\t", quote=FALSE, row.names=FALSE)

table_s12 <- file.path(tabdir, "Table_S12_HUMAnN_pathway_gene_family_association_tests.tsv")
write.table(all_res, table_s12, sep="\t", quote=FALSE, row.names=FALSE)

sig_summary <- aggregate(
  cbind(
    n_tests=!is.na(all_res$p_value),
    n_nominal=all_res$p_value < 0.05,
    n_fdr=all_res$q_value < 0.05
  ) ~ target + phenotype,
  data=all_res,
  FUN=sum,
  na.rm=TRUE
)
sumfile <- file.path(outdir, "Ronneby_HUMAnN_association_test_summary.tsv")
write.table(sig_summary, sumfile, sep="\t", quote=FALSE, row.names=FALSE)

phenofile <- file.path(outdir, "Ronneby_HUMAnN_association_phenotypes_used.tsv")
write.table(
  data.frame(
    phenotype=phenotypes,
    class=ifelse(phenotypes %in% categorical, "categorical", "numeric_log10_plus1"),
    stringsAsFactors=FALSE
  ),
  phenofile, sep="\t", quote=FALSE, row.names=FALSE
)

cat("DONE\n")
cat(out_all, "\n")
cat(sumfile, "\n")
cat(table_s12, "\n")
