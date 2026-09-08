args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Usage: 85A_rebuild_item3_toxicokinetic_summary.R ROOT OUTDIR METADATA")
ROOT <- normalizePath(args[[1]], mustWork=TRUE)
OUT <- args[[2]]
meta_path <- normalizePath(args[[3]], mustWork=TRUE)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"provenance"), recursive=TRUE, showWarnings=FALSE)

sep <- if (grepl("\\.csv$", meta_path, ignore.case=TRUE)) "," else "\t"
d <- read.table(meta_path, sep=sep, header=TRUE, check.names=FALSE,
                stringsAsFactors=FALSE, quote="\"", comment.char="")

kcols <- c("k_PFOA","k_PFPeS","k_PFHxS","k_PFHpS",
           "k_LPFOS","k_PFOS_MP1","k_PFOS_MP345","k_PFOS_MP26")
missing <- setdiff(kcols, names(d))
if (length(missing)) stop("Missing k columns: ", paste(missing, collapse=", "))

truthy <- function(x) tolower(trimws(as.character(x))) %in%
  c("1","true","yes","y","include","included")

if ("Include_elimination_rate_analysis" %in% names(d)) {
  keep <- truthy(d$Include_elimination_rate_analysis)
  eligibility_rule <- "Include_elimination_rate_analysis"
} else if ("ExposureGroup" %in% names(d)) {
  keep <- grepl("expos", d$ExposureGroup, ignore.case=TRUE)
  eligibility_rule <- "ExposureGroup contains expos"
} else {
  keep <- rowSums(!is.na(d[,kcols,drop=FALSE])) == length(kcols)
  eligibility_rule <- "all eight k parameters complete"
}

x <- d[keep,,drop=FALSE]
if (nrow(x) != 47) stop("Expected 47 eligible participants; found ", nrow(x))

for (k in kcols) {
  x[[k]] <- suppressWarnings(as.numeric(x[[k]]))
  if (sum(is.finite(x[[k]])) != 47) stop(k, " is not complete 47/47")
}

qfun <- function(v,p) as.numeric(quantile(v,p,na.rm=TRUE,names=FALSE,type=7))

dist_rows <- do.call(rbind,lapply(kcols,function(k){
  v <- x[[k]]
  data.frame(
    parameter=k,
    parameter_class="compound_specific_elimination_rate_constant",
    n=length(v), mean=mean(v), sd=sd(v), median=median(v),
    q1=qfun(v,0.25), q3=qfun(v,0.75), min=min(v), max=max(v),
    analysis_scope="47 PFAS-exposed participants with longitudinal elimination estimates",
    stringsAsFactors=FALSE
  )
}))

pairs <- combn(kcols,2,simplify=FALSE)
cor_rows <- do.call(rbind,lapply(pairs,function(pr){
  a <- x[[pr[[1]]]]
  b <- x[[pr[[2]]]]
  z <- suppressWarnings(cor.test(a,b,method="spearman",exact=FALSE))
  data.frame(
    parameter_x=pr[[1]], parameter_y=pr[[2]],
    n=sum(is.finite(a)&is.finite(b)),
    spearman_rho=unname(z$estimate), p_value=z$p.value,
    stringsAsFactors=FALSE
  )
}))
cor_rows$q_value_bh <- p.adjust(cor_rows$p_value,method="BH")
cor_rows$analysis_scope <-
  "47 PFAS-exposed participants with longitudinal elimination estimates"
cor_rows$multiplicity_note <-
  "BH correction across the 28 pairwise elimination-rate correlations"

fig1_candidates <- c(
  file.path(ROOT,"21_item3_toxicokinetic_audit",
            "item3_final_freeze_20260731_054949",
            "provenance","Figure_1D_toxicokinetic_correlation_summary.tsv"),
  file.path(ROOT,"submission","source_data","main","Figure_1",
            "Figure_1D_toxicokinetic_correlation_summary.tsv"),
  file.path(ROOT,"submission","source_data","main","Figure_01",
            "Figure_1D_toxicokinetic_correlation_summary.tsv")
)
fig1_path <- fig1_candidates[file.exists(fig1_candidates)][1]
if (is.na(fig1_path)) stop("Could not locate final Figure 1D toxicokinetic provenance")
fig1 <- read.delim(fig1_path,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
if (!all(c("compound","relationship","rho","n") %in% names(fig1)))
  stop("Figure 1D provenance schema mismatch")
if (nrow(fig1) < 15) stop("Figure 1D provenance unexpectedly sparse")

write.table(dist_rows,
  file.path(OUT,"tables","Table_S02a_PFAS_elimination_rate_distributions.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE,na="")
write.table(cor_rows,
  file.path(OUT,"tables","Table_S02b_PFAS_elimination_rate_correlations.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE,na="")
write.table(fig1,
  file.path(OUT,"tables","Table_S02c_Figure_1D_toxicokinetic_relationships.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE,na="")

combined_dist <- data.frame(
  table_section="elimination_rate_distribution",
  parameter=dist_rows$parameter, parameter_x="", parameter_y="",
  n=dist_rows$n, mean=dist_rows$mean, sd=dist_rows$sd,
  median=dist_rows$median, q1=dist_rows$q1, q3=dist_rows$q3,
  min=dist_rows$min, max=dist_rows$max,
  spearman_rho="", p_value="", q_value_bh="",
  analysis_scope=dist_rows$analysis_scope,
  stringsAsFactors=FALSE
)
combined_cor <- data.frame(
  table_section="elimination_rate_pairwise_correlation",
  parameter="", parameter_x=cor_rows$parameter_x,
  parameter_y=cor_rows$parameter_y, n=cor_rows$n,
  mean="",sd="",median="",q1="",q3="",min="",max="",
  spearman_rho=cor_rows$spearman_rho,
  p_value=cor_rows$p_value, q_value_bh=cor_rows$q_value_bh,
  analysis_scope=cor_rows$analysis_scope,
  stringsAsFactors=FALSE
)
write.table(rbind(combined_dist,combined_cor),
  file.path(OUT,"tables","Table_S02_final_PFAS_toxicokinetic_phenotype_summary.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE,na="")

coverage <- data.frame(
  parameter=kcols,
  complete_n=sapply(kcols,function(k) sum(is.finite(x[[k]]))),
  expected_n=47
)
write.table(coverage,
  file.path(OUT,"provenance","participant_parameter_coverage.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE)
file.copy(fig1_path,
  file.path(OUT,"provenance","Figure_1D_toxicokinetic_correlation_summary.tsv"),
  overwrite=TRUE)

writeLines(c(
  paste0("SELECTED_METADATA=",meta_path),
  paste0("ELIGIBILITY_RULE=",eligibility_rule),
  "ELIGIBLE_EXPOSED_PARTICIPANTS=47",
  "ELIMINATION_RATE_PARAMETERS=8",
  "NEW_TOXICOKINETIC_MODELS_FITTED=0"
), file.path(OUT,"provenance","T2K2_metadata_selection.txt"))

closeout <- data.frame(
  metric=c(
    "eligible_exposed_participants","continuous_rate_parameters",
    "parameters_complete_47_of_47","distribution_rows",
    "pairwise_correlation_rows","expected_pairwise_rows",
    "figure1_relationship_rows","new_toxicokinetic_models_fitted",
    "ITEM3_CLOSEOUT_STATUS"
  ),
  value=c(
    47,length(kcols),sum(coverage$complete_n==47),nrow(dist_rows),
    nrow(cor_rows),choose(length(kcols),2),nrow(fig1),0,
    ifelse(nrow(dist_rows)==8 && nrow(cor_rows)==28 &&
           all(coverage$complete_n==47),"PASS","FAIL")
  ),
  stringsAsFactors=FALSE
)
write.table(closeout,file.path(OUT,"item3_final_closeout_summary.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

cat("PHASE_T2K2_ITEM3_REBUILD=PASS\n")
cat("SELECTED_METADATA=",meta_path,"\n",sep="")
cat("ELIGIBLE_EXPOSED=47\n")
cat("DISTRIBUTION_ROWS=",nrow(dist_rows),"\n",sep="")
cat("PAIRWISE_CORRELATION_ROWS=",nrow(cor_rows),"\n",sep="")
cat("FIGURE1_RELATIONSHIP_ROWS=",nrow(fig1),"\n",sep="")
cat("NEW_TOXICOKINETIC_MODELS_FITTED=0\n")
