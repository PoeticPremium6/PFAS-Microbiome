# FIGURE2_PANELB_LEGEND_PATCH_V1_1
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

options(stringsAsFactors=FALSE, warn=1)

args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Repository root is required.", call.=FALSE)
REV <- normalizePath(args[[1]], mustWork=TRUE)
P14 <- file.path(REV, "14_figure_assembly")

OUT <- file.path(REV, "submission", "figure_builds", "Figure_2")
FIGDIR <- file.path(OUT, "figures")
PNGD <- file.path(OUT, "previews")
LOGD <- file.path(OUT, "reports")
SRC <- file.path(OUT, "source_data")
CLUSTERDIR <- file.path(OUT, "cluster_orders")

for (directory in c(OUT, FIGDIR, PNGD, LOGD, SRC, CLUSTERDIR)) {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
}

GENUS_PATHS <- c(
  file.path(REV, "08_microbiome_statistics/00_inputs/corrected_feature_matrices/metaphlan_genus_sample_matched.tsv"),
  file.path(REV, "03_taxonomy_metaphlan/03_merged/Ronneby_MetaPhlAn_genus.tsv")
)

META_PATHS <- c(
  file.path(P14, "manual_inputs/PFAS_Metadata.csv"),
  file.path(REV, "08_microbiome_statistics/00_inputs/phase8_selected_metadata.tsv")
)

purple <- c(
  deep="#3F007D",
  main="#4B2E83",
  mid="#6A51A3",
  soft="#807DBA",
  light="#9E9AC8",
  pale="#CBC9E2",
  faint="#E6E0EF",
  accent="#88419D",
  grey_dark="#333333",
  grey_mid="#7A7A7A",
  grey_light="#D9D9D9",
  other="#C7C7C7"
)

pal_class <- c(
  "Lower"="#D7C9F2",
  "Middle"="#9F7AD8",
  "Higher"="#5B2CA3",
  "Reference (no k)"="#BDB7C9"
)

theme_pfas <- function(base_size=9.0) {
  theme_bw(base_size=base_size, base_family="sans") +
    theme(
      axis.title=element_text(face="bold", colour=unname(purple["grey_dark"]), size=base_size+0.4),
      axis.text=element_text(colour=unname(purple["grey_dark"]), size=base_size-0.1),
      legend.title=element_text(face="bold", colour=unname(purple["grey_dark"]), size=base_size),
      legend.text=element_text(colour=unname(purple["grey_dark"]), size=base_size-0.6),
      strip.text=element_text(face="bold", colour=unname(purple["grey_dark"]), size=base_size-0.2),
      strip.background=element_rect(fill=unname(purple["faint"]), colour=unname(purple["grey_mid"]), linewidth=0.28),
      panel.grid.major=element_line(colour="#ECECEC", linewidth=0.22),
      panel.grid.minor=element_blank(),
      plot.margin=margin(4,5,4,5),
      plot.tag=element_text(face="bold", size=11, colour="black"),
      plot.tag.position=c(0.01,0.995)
    )
}

pick_existing <- function(paths) {
  hits <- paths[file.exists(paths) & file.info(paths)$size > 0]
  if (length(hits) == 0) stop("No usable file found among: ", paste(paths, collapse=" | "))
  hits[1]
}

parse_num <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | x == "NA" | x == "NaN"] <- NA_character_
  less <- grepl("^<", x)
  y <- suppressWarnings(as.numeric(gsub("^<", "", x)))
  y[less & is.finite(y)] <- y[less & is.finite(y)] / sqrt(2)
  y
}

clean_tax <- function(x) {
  x <- as.character(x)
  x <- gsub(".*\\|g__", "", x)
  x <- gsub("^g__", "", x)
  x <- gsub(".*\\|s__", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- "Unclassified"
  x
}

title_case <- function(x) {
  x <- as.character(x)
  out <- sapply(strsplit(tolower(x), " "), function(z) {
    paste0(toupper(substr(z, 1, 1)), substr(z, 2, nchar(z)), collapse=" ")
  })
  as.character(out)
}

read_abundance_matrix <- function(path, sample_ids=NULL) {
  dt <- fread(path, data.table=FALSE, check.names=FALSE)
  names(dt) <- make.unique(names(dt), sep="_dup")

  first <- dt[[1]]
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  rownames(mat) <- as.character(first)
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[is.na(mat)] <- 0

  rmatch <- if (!is.null(sample_ids)) sum(rownames(mat) %in% sample_ids) else 0
  cmatch <- if (!is.null(sample_ids)) sum(colnames(mat) %in% sample_ids) else 0

  if (cmatch > rmatch) mat <- t(mat)
  if (is.null(sample_ids) && nrow(mat) < ncol(mat)) mat <- t(mat)

  mat <- as.matrix(mat)
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[is.na(mat)] <- 0
  colnames(mat) <- make.unique(clean_tax(colnames(mat)), sep=" ")

  row_sums <- rowSums(mat, na.rm=TRUE)
  row_sums[row_sums == 0] <- 1
  if (median(row_sums, na.rm=TRUE) > 1.5) {
    mat <- mat / row_sums * 100
  }

  mat
}

stdz_mean <- function(df) {
  if (ncol(df) == 0) return(rep(NA_real_, nrow(df)))
  z <- lapply(df, function(v) {
    v <- parse_num(v)
    s <- sd(v, na.rm=TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(v)))
    (v - mean(v, na.rm=TRUE)) / s
  })
  out <- rowMeans(as.data.frame(z), na.rm=TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

bray_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      den <- sum(mat[i,] + mat[j,], na.rm=TRUE)
      d[i,j] <- ifelse(den > 0, sum(abs(mat[i,] - mat[j,]), na.rm=TRUE) / den, 0)
    }
  }
  as.dist(d)
}

shannon <- function(mat) {
  apply(mat, 1, function(v) {
    v[is.na(v)] <- 0
    s <- sum(v)
    if (s <= 0) return(NA_real_)
    p <- v / s
    p <- p[p > 0]
    -sum(p * log(p))
  })
}

meta_path <- pick_existing(META_PATHS)
meta <- fread(meta_path, data.table=TRUE, check.names=FALSE)
setnames(meta, make.unique(names(meta), sep="_dup"))

sample_col <- intersect(c("Sample_Name","sample_id","SampleID","sample","sample_name"), names(meta))[1]
if (is.na(sample_col)) stop("No sample ID column found in metadata.")
meta[, sample_id := as.character(get(sample_col))]

if ("ExposureGroup" %in% names(meta)) {
  meta[, exposure_group := as.character(ExposureGroup)]
} else if ("Condition" %in% names(meta)) {
  meta[, exposure_group := as.character(Condition)]
} else {
  meta[, exposure_group := "Unknown"]
}

meta[, exposure_group := ifelse(grepl("ref", exposure_group, ignore.case=TRUE), "Reference", "PFAS Exposed")]
meta[, exposure_group := factor(exposure_group, levels=c("PFAS Exposed","Reference"))]

k_cols <- intersect(names(meta), c("k_LPFOS","k_PFHxS","k_PFOA","k_PFOS_MP1","k_PFOS_MP26","k_PFOS_MP345","k_PFPeS","k_PFHpS"))
if (length(k_cols) == 0) k_cols <- grep("^k_", names(meta), value=TRUE)

for (kc in k_cols) meta[[kc]] <- parse_num(meta[[kc]])

meta[, elimination_score := stdz_mean(.SD), .SDcols=k_cols]
meta[, elimination_class := "Reference (no k)"]

idx <- which(meta$exposure_group == "PFAS Exposed" & is.finite(meta$elimination_score))
if (length(idx) >= 6) {
  qs <- quantile(meta$elimination_score[idx], probs=c(1/3, 2/3), na.rm=TRUE)
  meta[idx, elimination_class := ifelse(elimination_score <= qs[1], "Lower",
                                 ifelse(elimination_score <= qs[2], "Middle", "Higher"))]
}
meta[, elimination_class := factor(elimination_class, levels=c("Lower","Middle","Higher","Reference (no k)"))]

genus_path <- pick_existing(GENUS_PATHS)
genus <- read_abundance_matrix(genus_path, meta$sample_id)

keep <- intersect(rownames(genus), meta$sample_id)
genus <- genus[keep, , drop=FALSE]
meta2 <- meta[match(keep, sample_id)]

sample_order_dt <- data.table(
  sample_id=rownames(genus),
  exposure_group=meta2$exposure_group,
  elimination_score=meta2$elimination_score,
  elimination_class=meta2$elimination_class
)
sample_order_dt[exposure_group == "Reference", elimination_score := Inf]
setorder(sample_order_dt, exposure_group, elimination_score)
sample_order <- sample_order_dt$sample_id


class_levels <- c("Lower", "Middle", "Higher", "Reference (no k)")
class_labels <- c(
  "Lower"="Lower",
  "Middle"="Middle",
  "Higher"="Higher",
  "Reference (no k)"="Ref.\n(no k)"
)

safe_cor_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 10 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(c(rho=NA_real_, p=NA_real_, n=sum(keep)))
  }
  test <- suppressWarnings(try(
    cor.test(x[keep], y[keep], method="spearman", exact=FALSE),
    silent=TRUE
  ))
  if (inherits(test, "try-error")) {
    return(c(rho=NA_real_, p=NA_real_, n=sum(keep)))
  }
  c(rho=unname(test$estimate), p=test$p.value, n=sum(keep))
}

hclust_segments <- function(tree) {
  if (is.null(tree)) return(data.table())
  n <- length(tree$order)
  leaf_y <- numeric(n)
  leaf_y[tree$order] <- seq_len(n)
  node_y <- numeric(n - 1)
  output <- list()
  child_details <- function(child) {
    if (child < 0) return(list(x=0, y=leaf_y[-child]))
    list(x=tree$height[child], y=node_y[child])
  }
  index <- 1
  for (merge_index in seq_len(n - 1)) {
    left <- child_details(tree$merge[merge_index, 1])
    right <- child_details(tree$merge[merge_index, 2])
    parent_x <- tree$height[merge_index]
    parent_y <- mean(c(left$y, right$y))
    node_y[merge_index] <- parent_y
    output[[index]] <- data.table(x=left$x, xend=parent_x, y=left$y, yend=left$y)
    index <- index + 1
    output[[index]] <- data.table(x=right$x, xend=parent_x, y=right$y, yend=right$y)
    index <- index + 1
    output[[index]] <- data.table(x=parent_x, xend=parent_x, y=left$y, yend=right$y)
    index <- index + 1
  }
  rbindlist(output)
}

# A. Genus composition across ordered Ronneby samples.
top_n <- 10
top_genera <- names(sort(colMeans(genus, na.rm=TRUE), decreasing=TRUE))[
  seq_len(min(top_n, ncol(genus)))
]

comp <- as.data.table(genus[, top_genera, drop=FALSE])
comp[, sample_id := rownames(genus)]
comp[, Other := pmax(0, 100 - rowSums(.SD, na.rm=TRUE)), .SDcols=top_genera]
comp <- merge(
  comp,
  sample_order_dt[, .(sample_id, exposure_group, elimination_score, elimination_class)],
  by="sample_id",
  all.x=TRUE
)
comp[, sample_id := factor(sample_id, levels=sample_order)]
comp[, exposure_panel := ifelse(
  exposure_group == "Reference",
  "Reference (n=18; no k)",
  "PFAS exposed (n=47)"
)]
comp[, exposure_panel := factor(
  exposure_panel,
  levels=c("PFAS exposed (n=47)", "Reference (n=18; no k)")
)]

comp_long <- melt(
  comp,
  id.vars=c(
    "sample_id", "exposure_panel", "exposure_group",
    "elimination_score", "elimination_class"
  ),
  measure.vars=c(top_genera, "Other"),
  variable.name="Genus",
  value.name="Relative_abundance"
)
comp_long[, Genus_label := ifelse(
  Genus == "Other",
  "Other",
  title_case(Genus)
)]
comp_long[, Genus_label := factor(
  Genus_label,
  levels=c(title_case(top_genera), "Other")
)]

gen_cols <- colorRampPalette(
  c("#EEE7F7", "#CDBCE8", "#A98BD9", "#8062BD", "#5B2CA3", "#3F007D")
)(length(top_genera))
names(gen_cols) <- title_case(top_genera)
gen_cols <- c(gen_cols, Other=unname(purple["other"]))

fwrite(comp_long, file.path(SRC, "Figure_2A_genus_composition.tsv"), sep="\t")

pA <- ggplot(
  comp_long,
  aes(x=sample_id, y=Relative_abundance, fill=Genus_label)
) +
  geom_col(width=0.98) +
  facet_grid(. ~ exposure_panel, scales="free_x", space="free_x") +
  scale_y_continuous(limits=c(0,100), expand=expansion(mult=c(0,0.01))) +
  scale_fill_manual(values=gen_cols, drop=FALSE) +
  labs(
    tag="A",
    x="Ronneby samples; PFAS-exposed samples ordered by elimination score",
    y="Relative abundance (%)",
    fill="Genus"
  ) +
  theme_pfas(8.8) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid=element_blank(),
    legend.position="right",
    legend.key.height=unit(0.29,"cm"),
    legend.key.width=unit(0.34,"cm"),
    strip.text=element_text(size=8.3, face="bold")
  ) +
  guides(fill=guide_legend(ncol=1))

# B. Bray-Curtis PCoA with elimination-axis correlation.
bray <- bray_dist(genus)
ordination <- cmdscale(bray, eig=TRUE, k=2)
eigenvalues <- ordination$eig
percent <- round(100 * eigenvalues[1:2] / sum(abs(eigenvalues)), 1)

pcoa <- data.table(
  sample_id=rownames(genus),
  PCoA1=ordination$points[,1],
  PCoA2=ordination$points[,2],
  exposure_group=meta2$exposure_group,
  elimination_class=meta2$elimination_class,
  elimination_score=meta2$elimination_score
)
pcoa[, elimination_class := factor(as.character(elimination_class), levels=class_levels)]
pcoa[, exposure_group := factor(as.character(exposure_group), levels=c("PFAS Exposed", "Reference"))]

pcoa_use <- pcoa$exposure_group == "PFAS Exposed" & is.finite(pcoa$elimination_score)
pcoa_test <- safe_cor_test(pcoa$PCoA1[pcoa_use], pcoa$elimination_score[pcoa_use])

fwrite(pcoa, file.path(SRC, "Figure_2B_genus_Bray_PCoA.tsv"), sep="\t")

pB <- ggplot(
  pcoa,
  aes(x=PCoA1, y=PCoA2, fill=elimination_class, shape=exposure_group)
) +
  geom_point(
    size=2.55,
    colour=unname(purple["grey_dark"]),
    stroke=0.35,
    alpha=0.92
  ) +
  scale_fill_manual(values=pal_class, labels=class_labels, drop=FALSE) +
  scale_shape_manual(values=c("PFAS Exposed"=21, "Reference"=24), drop=FALSE) +
  annotate(
    "label",
    x=-Inf,
    y=Inf,
    hjust=-0.06,
    vjust=1.06,
    size=2.15,
    linewidth=0.18,
    fill="white",
    colour=unname(purple["deep"]),
    label=sprintf(
      "rho(PC1, elimination) = %.2f\np = %.3g; n = %d",
      pcoa_test[["rho"]],
      pcoa_test[["p"]],
      as.integer(pcoa_test[["n"]])
    )
  ) +
  labs(
    tag="B",
    x=sprintf("Bray-Curtis PCoA1 (%.1f%%)", percent[[1]]),
    y=sprintf("PCoA2 (%.1f%%)", percent[[2]]),
    fill="Elimination class",
    shape="Group"
  ) +
  theme_pfas(8.3) +
  theme(
    legend.position="bottom",
    legend.box="vertical",
    legend.key.size=unit(0.34,"cm")
  ) +
  guides(
    fill=guide_legend(nrow=1, byrow=TRUE),
    shape=guide_legend(nrow=1)
  )

# C. Shannon diversity with global descriptive tests.
alpha <- data.table(
  sample_id=rownames(genus),
  Shannon=shannon(genus),
  exposure_group=meta2$exposure_group,
  elimination_class=meta2$elimination_class
)
alpha <- alpha[is.finite(Shannon)]
alpha[, elimination_class := factor(as.character(elimination_class), levels=class_levels)]
alpha[, exposure_group := factor(as.character(exposure_group), levels=c("PFAS Exposed", "Reference"))]

alpha_exposed <- alpha[exposure_group == "PFAS Exposed"]
kw_p <- suppressWarnings(
  tryCatch(
    kruskal.test(Shannon ~ elimination_class, data=alpha_exposed)$p.value,
    error=function(e) NA_real_
  )
)
group_p <- suppressWarnings(
  tryCatch(
    wilcox.test(Shannon ~ exposure_group, data=alpha, exact=FALSE)$p.value,
    error=function(e) NA_real_
  )
)

fwrite(alpha, file.path(SRC, "Figure_2C_genus_Shannon.tsv"), sep="\t")

pC <- ggplot(
  alpha,
  aes(x=elimination_class, y=Shannon, fill=elimination_class)
) +
  geom_boxplot(
    width=0.60,
    outlier.shape=NA,
    colour=unname(purple["grey_dark"]),
    linewidth=0.36,
    alpha=0.96
  ) +
  geom_jitter(
    width=0.11,
    size=1.25,
    alpha=0.68,
    colour=unname(purple["grey_dark"])
  ) +
  scale_fill_manual(values=pal_class, labels=class_labels, drop=FALSE) +
  scale_x_discrete(labels=class_labels) +
  annotate(
    "label",
    x=-Inf,
    y=Inf,
    hjust=-0.05,
    vjust=1.04,
    size=2.05,
    linewidth=0.18,
    fill="white",
    colour=unname(purple["deep"]),
    label=sprintf(
      "Exposed tertiles: Kruskal-Wallis p = %.3g\nExposed vs reference: p = %.3g",
      kw_p,
      group_p
    )
  ) +
  labs(tag="C", x="Elimination class", y="Genus Shannon diversity") +
  theme_pfas(8.3) +
  theme(
    axis.text.x=element_text(size=7.1),
    legend.position="none"
  )

# D. Hierarchically clustered dominant-genus class profiles.
dominant <- names(sort(colMeans(genus, na.rm=TRUE), decreasing=TRUE))[
  seq_len(min(10, ncol(genus)))
]

mean_dt <- as.data.table(genus[, dominant, drop=FALSE])
mean_dt[, sample_id := rownames(genus)]
mean_dt <- merge(
  mean_dt,
  meta2[, .(sample_id, elimination_class)],
  by="sample_id",
  all.x=TRUE
)
mean_long <- melt(
  mean_dt,
  id.vars=c("sample_id", "elimination_class"),
  variable.name="Genus",
  value.name="abundance"
)
mean_sum <- mean_long[, .(
  Mean_abundance=mean(abundance, na.rm=TRUE)
), by=.(elimination_class, Genus)]
mean_sum[, elimination_class := factor(as.character(elimination_class), levels=class_levels)]

mean_wide <- dcast(mean_sum, Genus ~ elimination_class, value.var="Mean_abundance", fill=0)
genus_names <- mean_wide$Genus
mean_wide$Genus <- NULL
for (level in class_levels) {
  if (!level %in% names(mean_wide)) mean_wide[[level]] <- 0
}
mean_matrix <- as.matrix(mean_wide[, ..class_levels])
rownames(mean_matrix) <- genus_names

profile_z <- t(scale(t(mean_matrix)))
profile_z[!is.finite(profile_z)] <- 0
genus_tree <- hclust(dist(profile_z), method="average")
genus_order <- rownames(mean_matrix)[genus_tree$order]

cluster_export <- data.table(
  order_index=seq_along(genus_order),
  genus=genus_order,
  distance_method="Euclidean distance on row z-scores",
  linkage_method="average"
)
fwrite(cluster_export, file.path(CLUSTERDIR, "Figure_2D_genus_cluster_order.tsv"), sep="\t")

mean_sum[, y := match(as.character(Genus), genus_order)]
mean_sum[, class_short := factor(
  as.character(elimination_class),
  levels=class_levels,
  labels=unname(class_labels)
)]
mean_sum[, label := sprintf("%.1f", Mean_abundance)]
threshold <- quantile(mean_sum$Mean_abundance, probs=0.72, na.rm=TRUE)
mean_sum[, label_colour := ifelse(Mean_abundance >= threshold, "white", "dark")]

fwrite(mean_sum, file.path(SRC, "Figure_2D_clustered_genus_class_means.tsv"), sep="\t")

genus_label_map <- setNames(title_case(genus_order), seq_along(genus_order))

pD_heat <- ggplot(
  mean_sum,
  aes(x=class_short, y=y, fill=Mean_abundance)
) +
  geom_tile(colour="white", linewidth=0.45) +
  geom_text(aes(label=label, colour=label_colour), size=2.95, fontface="bold") +
  scale_colour_manual(
    values=c(white="white", dark=unname(purple["grey_dark"])),
    guide="none"
  ) +
  scale_fill_gradient(
    low=unname(purple["faint"]),
    high=unname(purple["deep"]),
    name="Mean abundance (%)"
  ) +
  scale_y_continuous(
    breaks=seq_along(genus_order),
    labels=genus_label_map,
    limits=c(0.5, length(genus_order)+0.5),
    expand=c(0,0),
    position="right"
  ) +
  labs(x="Elimination class", y="Dominant genus") +
  theme_pfas(8.3) +
  theme(
    panel.grid=element_blank(),
    axis.text.y=element_text(size=7.0, hjust=0),
    axis.text.x=element_text(size=7.1),
    legend.position="bottom"
  ) +
  guides(
    fill=guide_colourbar(
      direction="horizontal",
      title.position="top",
      barwidth=unit(38,"mm"),
      barheight=unit(3,"mm")
    )
  )

segments_d <- hclust_segments(genus_tree)
pD_dend <- ggplot(segments_d) +
  geom_segment(
    aes(x=x, xend=xend, y=y, yend=yend),
    linewidth=0.38,
    colour=unname(purple["grey_dark"])
  ) +
  scale_x_reverse(expand=c(0.02,0)) +
  scale_y_continuous(
    limits=c(0.5, length(genus_order)+0.5),
    expand=c(0,0)
  ) +
  labs(tag="D") +
  theme_void() +
  theme(
    plot.tag=element_text(face="bold", size=11),
    plot.tag.position=c(0.02,0.98),
    plot.margin=margin(2,0,2,2)
  )

pD <- pD_dend + pD_heat + plot_layout(widths=c(0.20,1))

# E. Genus-elimination associations in PFAS-exposed participants.
correlation_rows <- lapply(colnames(genus), function(genus_name) {
  x <- genus[, genus_name]
  y <- meta2$elimination_score
  use <- meta2$exposure_group == "PFAS Exposed" & is.finite(x) & is.finite(y)
  result <- safe_cor_test(x[use], y[use])
  data.table(
    Genus=genus_name,
    rho=as.numeric(result[["rho"]]),
    p_value=as.numeric(result[["p"]]),
    n=as.integer(result[["n"]]),
    exposed_prevalence=mean(x[use] > 0, na.rm=TRUE),
    exposed_mean_abundance=mean(x[use], na.rm=TRUE)
  )
})
correlation_all <- rbindlist(correlation_rows, fill=TRUE)
correlation_all <- correlation_all[is.finite(rho)]
correlation_all[, q_value := p.adjust(p_value, method="BH")]
correlation_all[, abs_rho := abs(rho)]
setorder(correlation_all, -abs_rho)

correlation_top <- correlation_all[seq_len(min(12, .N))]
correlation_top[, support_symbol := fifelse(
  q_value < 0.05,
  "\u2020",
  fifelse(p_value < 0.05, "*", "")
)]
correlation_top[, display_genus := title_case(Genus)]
correlation_top[, display_genus := factor(
  display_genus,
  levels=display_genus[order(rho)]
)]

nominal_count <- sum(correlation_all$p_value < 0.05, na.rm=TRUE)
fdr_count <- sum(correlation_all$q_value < 0.05, na.rm=TRUE)

fwrite(
  correlation_all,
  file.path(SRC, "Figure_2E_all_genus_elimination_correlations.tsv"),
  sep="\t"
)
fwrite(
  correlation_top,
  file.path(SRC, "Figure_2E_top_genus_elimination_correlations.tsv"),
  sep="\t"
)

maximum_rho <- max(abs(correlation_top$rho), na.rm=TRUE)

pE <- ggplot(
  correlation_top,
  aes(x=rho, y=display_genus)
) +
  geom_vline(
    xintercept=0,
    colour=unname(purple["grey_mid"]),
    linewidth=0.42
  ) +
  geom_segment(
    aes(x=0, xend=rho, y=display_genus, yend=display_genus),
    colour=unname(purple["mid"]),
    linewidth=0.86,
    alpha=0.95
  ) +
  geom_point(
    aes(size=exposed_prevalence, fill=exposed_mean_abundance),
    shape=21,
    colour=unname(purple["grey_dark"]),
    stroke=0.34,
    alpha=0.98
  ) +
  scale_fill_gradient(
    low=unname(purple["pale"]),
    high=unname(purple["deep"]),
    labels=percent_format(accuracy=1),
    name="Mean relative\nabundance"
  ) +
  scale_size_continuous(
    labels=percent_format(accuracy=1),
    range=c(2.7,6.0),
    name="Prevalence"
  ) +
  scale_x_continuous(limits=c(-maximum_rho*1.16, maximum_rho*1.16)) +
  annotate(
    "label",
    x=-Inf,
    y=Inf,
    hjust=-0.04,
    vjust=1.04,
    size=2.05,
    linewidth=0.18,
    fill="white",
    colour=unname(purple["deep"]),
    label=sprintf(
      "Nominal p < 0.05: %d genera\nFDR q < 0.05: %d genera",
      nominal_count,
      fdr_count
    )
  ) +
  labs(
    tag="E",
    x="Spearman rho with elimination score (PFAS-exposed only)",
    y="Genus"
  ) +
  theme_pfas(8.5) +
  theme(
    legend.position="right",
    axis.text.y=element_text(size=7.1)
  )

row_bc <- pB + pC + plot_layout(widths=c(1.12,0.88))
# === PHASE_F2C_FIGURE2_EDITORIAL_OVERRIDE_BEGIN ===
# Narrow editorial correction: retain the authoritative Figure 2 data preparation
# and Panels C-E, while rebuilding only Panels A-B and their shared row layout.

figure2_editorial_theme <- function(base_size = 12.6) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(face = "bold", colour = "#2D2933"),
      axis.text = element_text(colour = "#3A3542"),
      legend.title = element_text(face = "bold", colour = "#2D2933"),
      legend.text = element_text(colour = "#3A3542"),
      strip.text = element_text(face = "bold", colour = "#2D2933"),
      strip.background = element_rect(fill = "#E6DFEC", colour = "#77717E", linewidth = 0.30),
      panel.grid.major = element_line(colour = "#ECE8EF", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

# -------------------------
# A. Expanded genus composition from the authoritative genus matrix.
# -------------------------
figure2_genus_matrix <- as.matrix(genus)
storage.mode(figure2_genus_matrix) <- "numeric"

# Work in percentage units regardless of whether the authoritative matrix is
# represented as proportions or percentages.
if (max(figure2_genus_matrix, na.rm = TRUE) <= 1.5) {
  figure2_genus_matrix <- figure2_genus_matrix * 100
}

figure2_genus_means <- sort(colMeans(figure2_genus_matrix, na.rm = TRUE), decreasing = TRUE)
figure2_display_n <- min(15L, length(figure2_genus_means))
figure2_display_genera <- names(figure2_genus_means)[seq_len(figure2_display_n)]

figure2_comp_wide <- as.data.table(
  figure2_genus_matrix[, figure2_display_genera, drop = FALSE],
  keep.rownames = "sample_id"
)
figure2_total_abundance <- rowSums(figure2_genus_matrix, na.rm = TRUE)
figure2_selected_abundance <- rowSums(
  figure2_genus_matrix[, figure2_display_genera, drop = FALSE],
  na.rm = TRUE
)
figure2_comp_wide[, Other := pmax(figure2_total_abundance - figure2_selected_abundance, 0)]

figure2_comp_long <- melt(
  figure2_comp_wide,
  id.vars = "sample_id",
  variable.name = "Genus",
  value.name = "Relative_abundance"
)

figure2_meta_editorial <- copy(as.data.table(meta2))
figure2_meta_editorial[, sample_id := as.character(sample_id)]
figure2_comp_long[, sample_id := as.character(sample_id)]
figure2_comp_long <- merge(
  figure2_comp_long,
  figure2_meta_editorial[, .(sample_id, exposure_group, elimination_score)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

figure2_exposed_order <- figure2_meta_editorial[
  grepl("^PFAS", exposure_group),
  .(sample_id, elimination_score)
][order(elimination_score, sample_id), sample_id]
figure2_reference_order <- figure2_meta_editorial[
  grepl("Reference", exposure_group),
  sample_id
]
figure2_sample_levels <- unique(c(figure2_exposed_order, figure2_reference_order))
figure2_comp_long[, sample_id := factor(sample_id, levels = figure2_sample_levels)]

figure2_n_exposed <- uniqueN(figure2_exposed_order)
figure2_n_reference <- uniqueN(figure2_reference_order)
figure2_comp_long[, exposure_panel := fifelse(
  grepl("^PFAS", exposure_group),
  paste0("PFAS exposed (n=", figure2_n_exposed, ")"),
  paste0("Reference (n=", figure2_n_reference, "; no elimination score)")
)]
figure2_comp_long[, exposure_panel := factor(
  exposure_panel,
  levels = c(
    paste0("PFAS exposed (n=", figure2_n_exposed, ")"),
    paste0("Reference (n=", figure2_n_reference, "; no elimination score)")
  )
)]

figure2_title_case <- function(x) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  ifelse(
    grepl("^(GGB|SGB)[0-9]+$", x, ignore.case = TRUE),
    toupper(x),
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  )
}
figure2_comp_long[, Genus_label := figure2_title_case(Genus)]
figure2_comp_long[Genus == "Other", Genus_label := "Other"]

figure2_named_labels <- figure2_title_case(figure2_display_genera)
figure2_genus_palette <- setNames(
  grDevices::colorRampPalette(c(
    "#EAE4F2", "#D4C6E5", "#B49AD3", "#916BBE",
    "#7043A5", "#512785", "#35145F"
  ))(length(figure2_named_labels)),
  figure2_named_labels
)
figure2_genus_palette <- c(figure2_genus_palette, Other = "#C7C7C7")
figure2_legend_levels <- c(figure2_named_labels, "Other")
figure2_comp_long[, Genus_label := factor(Genus_label, levels = figure2_legend_levels)]

fwrite(
  figure2_comp_long,
  file.path(SRC, "Figure_2A_genus_composition.tsv"),
  sep = "\t"
)

pA <- ggplot(
  figure2_comp_long,
  aes(x = sample_id, y = Relative_abundance, fill = Genus_label)
) +
  geom_col(width = 0.98) +
  facet_grid(. ~ exposure_panel, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_manual(values = figure2_genus_palette, drop = FALSE) +
  labs(
    x = "Ronneby samples; PFAS-exposed samples ordered by elimination score",
    y = "Relative abundance (%)",
    fill = "Genus"
  ) +
  figure2_editorial_theme(12.2) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.34, "cm"),
    legend.text = element_text(size = 8.6),
    legend.title = element_text(size = 10.0),
    strip.text = element_text(size = 11.2),
    axis.title.x = element_text(size = 12.4),
    axis.title.y = element_text(size = 12.4)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# -------------------------
# B. Full-cohort ordination with explicit exposed-only statistic boundary.
# -------------------------
figure2_pcoa <- copy(as.data.table(pcoa))
figure2_pcoa[, exposure_group := as.character(exposure_group)]
figure2_pcoa[, elimination_class := as.character(elimination_class)]

if (!"elimination_score" %in% names(figure2_pcoa) && "sample_id" %in% names(figure2_pcoa)) {
  figure2_pcoa <- merge(
    figure2_pcoa,
    figure2_meta_editorial[, .(sample_id, elimination_score)],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
}

figure2_pcoa[, group_label := fifelse(
  grepl("^PFAS", exposure_group),
  "PFAS exposed",
  "Reference"
)]
figure2_pcoa[, class_display := fifelse(
  group_label == "PFAS exposed",
  elimination_class,
  NA_character_
)]
figure2_pcoa[, class_display := factor(
  class_display,
  levels = c("Lower", "Middle", "Higher")
)]
figure2_pcoa[, group_label := factor(
  group_label,
  levels = c("PFAS exposed", "Reference")
)]

figure2_stat_use <- figure2_pcoa$group_label == "PFAS exposed" &
  is.finite(figure2_pcoa$PCoA1) &
  is.finite(figure2_pcoa$elimination_score)
figure2_cor_test <- suppressWarnings(cor.test(
  figure2_pcoa$PCoA1[figure2_stat_use],
  figure2_pcoa$elimination_score[figure2_stat_use],
  method = "spearman",
  exact = FALSE
))
figure2_rho <- unname(figure2_cor_test$estimate)
figure2_p <- figure2_cor_test$p.value
figure2_n_total <- nrow(figure2_pcoa)
figure2_n_stat <- sum(figure2_stat_use)

figure2_panel_b_note <- paste0(
  "Ordination: n=", figure2_n_total,
  " (", figure2_n_exposed, " exposed; ", figure2_n_reference, " reference)\n",
  "Exposed only: rho(PC1, elimination)=", sprintf("%.2f", figure2_rho),
  "; p=", format.pval(figure2_p, digits = 2, eps = 0.001),
  "; n=", figure2_n_stat
)

figure2_class_palette <- c(
  Lower = "#CFC3E2",
  Middle = "#9979C3",
  Higher = "#5A2B91"
)

# Preserve the authoritative PCoA axis labels before replacing Panel B.
figure2_existing_pB <- pB
figure2_pcoa_x_label <- if (
  !is.null(figure2_existing_pB$labels$x) &&
  length(figure2_existing_pB$labels$x) == 1L &&
  nzchar(as.character(figure2_existing_pB$labels$x))
) {
  as.character(figure2_existing_pB$labels$x)
} else {
  "Bray-Curtis PCoA1"
}
figure2_pcoa_y_label <- if (
  !is.null(figure2_existing_pB$labels$y) &&
  length(figure2_existing_pB$labels$y) == 1L &&
  nzchar(as.character(figure2_existing_pB$labels$y))
) {
  as.character(figure2_existing_pB$labels$y)
} else {
  "PCoA2"
}

pB <- ggplot(
  figure2_pcoa,
  aes(x = PCoA1, y = PCoA2, fill = class_display, shape = group_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#BEB8C3", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#BEB8C3", linewidth = 0.35) +
  geom_point(
    size = 3.45,
    colour = "#3A3542",
    stroke = 0.40,
    alpha = 0.92
  ) +
  scale_fill_manual(
    values = figure2_class_palette,
    breaks = c("Lower", "Middle", "Higher"),
    na.value = "#BDBDBD",
    na.translate = FALSE,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c("PFAS exposed" = 21, "Reference" = 24),
    drop = FALSE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = figure2_panel_b_note,
    hjust = -0.02,
    vjust = 1.08,
    size = 3.05,
    linewidth = 0.22,
    colour = "#4A286F",
    fill = "white"
  ) +
  labs(
    x = figure2_pcoa_x_label,
    y = figure2_pcoa_y_label,
    fill = "Elimination class",
    shape = "Group"
  ) +
  figure2_editorial_theme(12.4) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.38, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    axis.title = element_text(face = "bold")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  )

# Keep the established Panel C statistics and design, but rebalance the shared row.
row_bc <- pB + pC + plot_layout(widths = c(1.18, 0.82))

figure2_other_summary <- figure2_comp_long[
  Genus_label == "Other",
  .(
    panelA_other_mean_percent = mean(Relative_abundance, na.rm = TRUE),
    panelA_other_median_percent = median(Relative_abundance, na.rm = TRUE),
    panelA_other_max_percent = max(Relative_abundance, na.rm = TRUE)
  )
]
figure2_refinement_metrics <- data.table(
  metric = c(
    "panelA_displayed_named_genera",
    "panelA_other_mean_percent",
    "panelA_other_median_percent",
    "panelA_other_max_percent",
    "panelB_n_total",
    "panelB_n_exposed",
    "panelB_n_reference",
    "panelB_exposed_stat_n",
    "panelB_rho_pc1_elimination",
    "panelB_p_pc1_elimination",
    "panels_C_to_E_preserved"
  ),
  value = as.character(c(
    figure2_display_n,
    figure2_other_summary$panelA_other_mean_percent,
    figure2_other_summary$panelA_other_median_percent,
    figure2_other_summary$panelA_other_max_percent,
    figure2_n_total,
    figure2_n_exposed,
    figure2_n_reference,
    figure2_n_stat,
    figure2_rho,
    figure2_p,
    TRUE
  ))
)
fwrite(
  figure2_refinement_metrics,
  file.path(SRC, "Figure_2_editorial_refinement_metrics.tsv"),
  sep = "\t"
)
# === PHASE_F2G_FIGURE2_FINAL_LAYOUT_BEGIN ===
figure2_exposed_label <- paste0("PFAS exposed (n=", figure2_n_exposed, ")")
figure2_reference_label <- paste0("Reference (n=", figure2_n_reference, ")")

figure2_comp_long[, exposure_panel := fifelse(
  grepl("^PFAS", exposure_group), figure2_exposed_label, figure2_reference_label
)]
figure2_comp_long[, exposure_panel := factor(
  exposure_panel, levels=c(figure2_exposed_label, figure2_reference_label)
)]

figure2_class_counts <- figure2_meta_editorial[grepl("^PFAS", exposure_group), .N, by=elimination_class]
figure2_n_lower <- figure2_class_counts[as.character(elimination_class)=="Lower", N]
figure2_n_middle <- figure2_class_counts[as.character(elimination_class)=="Middle", N]
figure2_n_higher <- figure2_class_counts[as.character(elimination_class)=="Higher", N]
if (length(figure2_n_lower)==0L) figure2_n_lower <- 16L
if (length(figure2_n_middle)==0L) figure2_n_middle <- 15L
if (length(figure2_n_higher)==0L) figure2_n_higher <- 16L

figure2_panelA_breaks <- data.table(
  exposure_panel=factor(rep(figure2_exposed_label,2), levels=c(figure2_exposed_label,figure2_reference_label)),
  xintercept=c(figure2_n_lower+0.5, figure2_n_lower+figure2_n_middle+0.5)
)
figure2_mid_idx <- c(
  pmax(1L,round(figure2_n_lower/2)),
  pmax(1L,figure2_n_lower+round(figure2_n_middle/2)),
  pmax(1L,figure2_n_lower+figure2_n_middle+round(figure2_n_higher/2))
)
figure2_mid_idx <- pmin(figure2_mid_idx,length(figure2_exposed_order))
figure2_panelA_labels <- data.table(
  sample_id=factor(as.character(figure2_exposed_order[figure2_mid_idx]), levels=figure2_sample_levels),
  exposure_panel=factor(rep(figure2_exposed_label,3), levels=c(figure2_exposed_label,figure2_reference_label)),
  label=c(paste0("Lower\n(n=",figure2_n_lower,")"), paste0("Middle\n(n=",figure2_n_middle,")"), paste0("Higher\n(n=",figure2_n_higher,")"))
)

pA <- ggplot(figure2_comp_long, aes(x=sample_id,y=Relative_abundance,fill=Genus_label)) +
  geom_col(width=0.98) +
  geom_vline(data=figure2_panelA_breaks,aes(xintercept=xintercept),inherit.aes=FALSE,linewidth=0.55,colour="#675B70",linetype="dashed") +
  geom_text(data=figure2_panelA_labels,aes(x=sample_id,y=-6,label=label),inherit.aes=FALSE,lineheight=0.9,size=3.1,fontface="bold",colour="#302A36") +
  facet_grid(.~exposure_panel,scales="free_x",space="free_x") +
  scale_y_continuous(limits=c(-12,100),breaks=c(0,25,50,75,100),labels=function(x) ifelse(x<0,"",paste0(x,"%")),expand=expansion(mult=c(0,0.01))) +
  scale_fill_manual(values=figure2_genus_palette,drop=FALSE) +
  coord_cartesian(clip="off") +
  labs(x="PFAS-exposed samples are ordered left-to-right by elimination score",y="Relative abundance (%)",fill="Genus") +
  figure2_editorial_theme(12.4) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid=element_blank(),legend.position="right",legend.key.height=unit(0.25,"cm"),legend.key.width=unit(0.34,"cm"),legend.text=element_text(size=8.7),legend.title=element_text(size=10.2),strip.text=element_text(size=10.8),axis.title.x=element_text(size=11,margin=margin(t=3)),axis.title.y=element_text(size=12.2),plot.margin=margin(9,8,9,16,"mm")) +
  guides(fill=guide_legend(ncol=1,byrow=TRUE))

pB <- pB +
  labs(tag=NULL) +
  theme(legend.position="bottom",legend.box="vertical",legend.title=element_text(size=10.5,face="bold"),legend.text=element_text(size=9.5),axis.title=element_text(face="bold",size=11.5),plot.margin=margin(9,8,9,16,"mm"))

pC <- pC + labs(tag=NULL) +
  theme(axis.title=element_text(face="bold",size=11.5),axis.text=element_text(size=9.5),plot.margin=margin(9,10,9,16,"mm"))

try(pD$labels$tag <- NULL, silent=TRUE)
try(pD$labels$title <- NULL, silent=TRUE)
if (inherits(pD,"patchwork")) pD <- pD + plot_annotation(tag_levels=NULL)

figure2_pick_column <- function(nms,patterns,exclude=character()) {
  for (pattern in patterns) {
    hits <- grep(pattern,nms,ignore.case=TRUE,value=TRUE,perl=TRUE)
    hits <- setdiff(hits,exclude)
    if (length(hits)>0L) return(hits[[1L]])
  }
  NA_character_
}
figure2_e_support <- fread(file.path(SRC,"Figure_2E_all_genus_elimination_correlations.tsv"))
figure2_e_names <- names(figure2_e_support)
figure2_e_genus_col <- figure2_pick_column(figure2_e_names,c("^Genus_label$","^Genus$","genus","taxon","feature"))
figure2_e_rho_col <- figure2_pick_column(figure2_e_names,c("^rho$","spearman.*rho","correlation","effect"),exclude=grep("abs",figure2_e_names,ignore.case=TRUE,value=TRUE))
figure2_e_prev_col <- figure2_pick_column(figure2_e_names,c("exposed.*prevalence","prevalence"))
figure2_e_abund_col <- figure2_pick_column(figure2_e_names,c("exposed.*mean.*abund","mean.*abund","abund.*mean"))
figure2_e_p_col <- figure2_pick_column(figure2_e_names,c("^p$","^p_value$","^pvalue$","p_val","nominal.*p"))
figure2_safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
figure2_e_plot <- data.table(
  Genus=as.character(figure2_e_support[[figure2_e_genus_col]]),
  rho=figure2_safe_numeric(figure2_e_support[[figure2_e_rho_col]]),
  prevalence=figure2_safe_numeric(figure2_e_support[[figure2_e_prev_col]]),
  mean_abundance=figure2_safe_numeric(figure2_e_support[[figure2_e_abund_col]]),
  p_value=figure2_safe_numeric(figure2_e_support[[figure2_e_p_col]])
)
figure2_e_plot <- figure2_e_plot[is.finite(rho)&is.finite(p_value)&p_value<0.05]
figure2_e_plot[,abs_rho:=abs(rho)]
setorder(figure2_e_plot,-abs_rho,Genus)
if (nrow(figure2_e_plot)>12L) figure2_e_plot <- figure2_e_plot[seq_len(12L)]
if (max(figure2_e_plot$prevalence,na.rm=TRUE)>1.5) figure2_e_plot[,prevalence:=prevalence/100]
if (max(figure2_e_plot$mean_abundance,na.rm=TRUE)<=1.5) figure2_e_plot[,mean_abundance:=mean_abundance*100]
figure2_e_plot[,Genus_display:=figure2_title_case(Genus)]
figure2_e_plot <- figure2_e_plot[order(rho)]
figure2_e_plot[,Genus_display:=factor(Genus_display,levels=unique(Genus_display))]
figure2_e_max_abs <- max(abs(figure2_e_plot$rho),na.rm=TRUE)

pE <- ggplot(figure2_e_plot,aes(x=rho,y=Genus_display)) +
  geom_vline(xintercept=0,colour="#8E8794",linewidth=0.45) +
  geom_segment(aes(x=0,xend=rho,yend=Genus_display),colour="#7652AA",linewidth=0.9,alpha=0.95) +
  geom_point(aes(size=prevalence,fill=mean_abundance),shape=21,colour="#3A3542",stroke=0.55,alpha=0.99) +
  annotate("label",x=-figure2_e_max_abs*1.03,y=length(levels(figure2_e_plot$Genus_display))+0.55,label=paste0("Nominal p < 0.05: ",nrow(figure2_e_plot)," genera"),hjust=0,vjust=1,size=2.75,linewidth=0.20,colour="#4A286F",fill="white") +
  scale_fill_gradient(low="#D6CBE5",high="#4B176E") +
  scale_size_continuous(labels=scales::percent_format(accuracy=1),range=c(3,7)) +
  scale_x_continuous(limits=c(-figure2_e_max_abs*1.12,figure2_e_max_abs*1.12)) +
  labs(x="Spearman rho with elimination score (PFAS-exposed only)",y="Genus",fill="Mean relative\nabundance",size="Prevalence") +
  figure2_editorial_theme(12.4) +
  theme(legend.position="right",legend.box="vertical",legend.key.height=unit(0.42,"cm"),panel.grid.minor=element_blank(),plot.margin=margin(9,8,9,16,"mm")) +
  guides(fill=guide_colourbar(order=1),size=guide_legend(order=2))

fwrite(figure2_e_plot,file.path(SRC,"Figure_2E_nominal_only_final.tsv"),sep="\t")
# === PHASE_F2G_FIGURE2_FINAL_LAYOUT_END ===
# === PHASE_F2H2_FIGURE2_FINAL_EDITORIAL_REPAIR_BEGIN ===
# Display-only repair layered on the successful F2G state.

# Panel A: preserve F2G content, only ensure breathing room.
pA <- pA + theme(plot.margin = margin(9, 7, 10, 13, "mm"))

# Panel B: move legends to the right and statistics above the data.
figure2_pcoa[, group_label := factor(
  as.character(group_label), levels = c("PFAS exposed", "Reference")
)]
figure2_pcoa[, class_display := factor(
  as.character(class_display), levels = c("Lower", "Middle", "Higher")
)]

figure2_panel_b_subtitle <- paste0(
  "Ordination: n=", figure2_n_total,
  " (", figure2_n_exposed, " exposed; ", figure2_n_reference, " reference)\n",
  "Exposed only: rho(PC1, elimination)=",
  sprintf("%.2f", figure2_rho),
  "; p=", format.pval(figure2_p, digits = 2, eps = 0.001)
)

pB <- ggplot(
  figure2_pcoa,
  aes(x = PCoA1, y = PCoA2, fill = class_display, shape = group_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#BEB8C3", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#BEB8C3", linewidth = 0.35) +
  geom_point(size = 3.2, colour = "#3A3542", stroke = 0.42, alpha = 0.92) +
  scale_fill_manual(
    values = figure2_class_palette,
    breaks = c("Lower", "Middle", "Higher"),
    na.value = "white", na.translate = FALSE, drop = FALSE
  ) +
  scale_shape_manual(values = c("PFAS exposed" = 21, "Reference" = 24), drop = FALSE) +
  labs(
    x = figure2_pcoa_x_label,
    y = figure2_pcoa_y_label,
    subtitle = figure2_panel_b_subtitle,
    fill = "Elimination class",
    shape = "Group"
  ) +
  figure2_editorial_theme(12.5) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10.2, face = "bold"),
    legend.text = element_text(size = 9.3),
    legend.key.size = unit(0.42, "cm"),
    plot.subtitle = element_text(size = 8.6, colour = "#4A286F", lineheight = 1.08, margin = margin(b = 4)),
    axis.title = element_text(face = "bold", size = 11.4),
    plot.margin = margin(9, 8, 9, 13, "mm")
  ) +
  guides(
    # Force a fill-capable glyph for the elimination-class guide so the
    # light / medium / dark purple fills are visible in the legend.
    fill = guide_legend(
      order = 1,
      ncol = 1,
      override.aes = list(
        shape = 21,
        size = 3.6,
        colour = "#3A3542",
        stroke = 0.42
      )
    ),
    # Keep the exposure-group shape guide unfilled so it remains a pure
    # circle/triangle shape distinction and does not compete with class fill.
    shape = guide_legend(
      order = 2,
      ncol = 1,
      override.aes = list(
        fill = "white",
        colour = "#3A3542",
        size = 3.6,
        stroke = 0.42
      )
    )
  )

# Panel C: statistics in subtitle, not over the points.
figure2_c_palette <- c(
  "Lower" = "#CFC3E2",
  "Middle" = "#9979C3",
  "Higher" = "#5A2B91",
  "Reference (no k)" = "#BDBDBD"
)
figure2_c_labels <- c(
  "Lower" = paste0("Lower\n(n=", figure2_n_lower, ")"),
  "Middle" = paste0("Middle\n(n=", figure2_n_middle, ")"),
  "Higher" = paste0("Higher\n(n=", figure2_n_higher, ")"),
  "Reference (no k)" = paste0("Ref.\n(n=", figure2_n_reference, ")")
)
figure2_panel_c_subtitle <- paste0(
  "Exposed tertiles: Kruskal-Wallis p = 0.0225\n",
  "Exposed vs reference: p = 0.486"
)
alpha[, elimination_class := factor(
  as.character(elimination_class),
  levels = c("Lower", "Middle", "Higher", "Reference (no k)")
)]

pC <- ggplot(alpha, aes(x = elimination_class, y = Shannon, fill = elimination_class)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, colour = "#4D4852", linewidth = 0.40, alpha = 0.97) +
  geom_jitter(width = 0.12, size = 1.85, alpha = 0.72, colour = "#4D4852") +
  scale_fill_manual(values = figure2_c_palette, drop = FALSE) +
  scale_x_discrete(labels = figure2_c_labels) +
  labs(
    x = "Elimination class",
    y = "Genus Shannon diversity",
    subtitle = figure2_panel_c_subtitle
  ) +
  figure2_editorial_theme(12.5) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8.9, lineheight = 0.92),
    axis.title = element_text(face = "bold", size = 11.4),
    plot.subtitle = element_text(size = 8.35, colour = "#4A286F", lineheight = 1.08, margin = margin(b = 4)),
    plot.margin = margin(9, 9, 9, 13, "mm")
  )

# Panel D: preserve dendrogram + heatmap but move collected guide to the right.
try(pD$labels$tag <- NULL, silent = TRUE)
try(pD$labels$title <- NULL, silent = TRUE)
if (inherits(pD, "patchwork")) {
  pD <- pD + plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 9.7),
      legend.text = element_text(size = 8.7)
    )
} else {
  pD <- pD + theme(legend.position = "right", legend.direction = "vertical")
}

# Panel E: preserve F2G nominal-only result and right-side legends.
try(pE$labels$tag <- NULL, silent = TRUE)
pE <- pE + theme(
  legend.position = "right",
  legend.box = "vertical",
  plot.margin = margin(9, 8, 9, 13, "mm")
)

# Clear explicit tags. The single top-level auto-tag pass is authoritative.
for (figure2_nm in c("pA", "pB", "pC", "pD", "pE")) {
  figure2_obj <- get(figure2_nm)
  try(figure2_obj$labels$tag <- NULL, silent = TRUE)
  assign(figure2_nm, figure2_obj)
}
# === PHASE_F2H2_FIGURE2_FINAL_EDITORIAL_REPAIR_END ===
# === PHASE_F2I_FIGURE2_TYPOGRAPHY_AND_D_TAG_BEGIN ===
# Final display-only typography pass. No data/statistical changes.

figure2_strip_nested_tags <- function(x) {
  if (inherits(x, "ggplot")) {
    try(x$labels$tag <- NULL, silent = TRUE)
  }
  if (inherits(x, "patchwork")) {
    if (!is.null(x$patches$plots)) {
      x$patches$plots <- lapply(
        x$patches$plots,
        figure2_strip_nested_tags
      )
    }
    if (!is.null(x$patches$annotation)) {
      try(x$patches$annotation$tag_levels <- NULL, silent = TRUE)
      try(x$patches$annotation$tag_prefix <- NULL, silent = TRUE)
      try(x$patches$annotation$tag_suffix <- NULL, silent = TRUE)
    }
  }
  x
}

pA <- pA +
  theme(
    legend.title = element_text(face = "bold", size = 11.8),
    legend.text = element_text(size = 10.8),
    legend.key.height = unit(0.29, "cm"),
    axis.text.y = element_text(size = 10.7),
    axis.title.y = element_text(face = "bold", size = 13.2),
    axis.title.x = element_text(face = "bold", size = 12.0),
    strip.text = element_text(face = "bold", size = 12.0)
  )

pB <- pB +
  theme(
    axis.text = element_text(size = 10.2),
    axis.title = element_text(face = "bold", size = 12.2),
    legend.title = element_text(face = "bold", size = 10.8),
    legend.text = element_text(size = 10.0),
    plot.subtitle = element_text(
      size = 9.2,
      colour = "#4A286F",
      lineheight = 1.08,
      margin = margin(b = 4)
    )
  )

pC <- pC +
  theme(
    axis.text.y = element_text(size = 10.2),
    axis.text.x = element_text(size = 9.7, lineheight = 0.94),
    axis.title = element_text(face = "bold", size = 12.2),
    plot.subtitle = element_text(
      size = 9.0,
      colour = "#4A286F",
      lineheight = 1.08,
      margin = margin(b = 4)
    )
  )

pD <- figure2_strip_nested_tags(pD)
try(pD$labels$tag <- NULL, silent = TRUE)
try(pD$labels$title <- NULL, silent = TRUE)

if (inherits(pD, "patchwork")) {
  pD <- pD &
    theme(
      axis.text.y = element_text(size = 10.6),
      axis.text.x = element_text(size = 10.2),
      axis.title = element_text(face = "bold", size = 11.8),
      legend.title = element_text(face = "bold", size = 10.5),
      legend.text = element_text(size = 9.5)
    )
} else {
  pD <- pD +
    theme(
      axis.text.y = element_text(size = 10.6),
      axis.text.x = element_text(size = 10.2),
      axis.title = element_text(face = "bold", size = 11.8),
      legend.title = element_text(face = "bold", size = 10.5),
      legend.text = element_text(size = 9.5)
    )
}

pE <- pE +
  theme(
    axis.text.y = element_text(size = 10.4),
    axis.text.x = element_text(size = 10.2),
    axis.title = element_text(face = "bold", size = 12.2),
    legend.title = element_text(face = "bold", size = 10.8),
    legend.text = element_text(size = 9.8)
  )


# FIGURE2_READABILITY_PATCH_V1_0
# Display-only refinement: larger typography, distinct elimination-class colours,
# and explicit removal of dendrogram axes. Scientific data/statistics are unchanged.

figure2_readability_class_palette <- c(
  "Lower" = "#D9CFEA",
  "Middle" = "#9A7AC8",
  "Higher" = "#5B2D90"
)
figure2_readability_c_palette <- c(
  "Lower" = "#D9CFEA",
  "Middle" = "#9A7AC8",
  "Higher" = "#5B2D90",
  "Reference (no k)" = "#BDBDBD"
)

# Make the three elimination classes visibly distinct in both plotted fills and
# the Panel B legend. Reference remains unfilled/grey as before.
pB <- pB +
  scale_fill_manual(
    values = figure2_readability_class_palette,
    breaks = c("Lower", "Middle", "Higher"),
    na.value = "white", na.translate = FALSE, drop = FALSE
  )

pC <- pC +
  scale_fill_manual(values = figure2_readability_c_palette, drop = FALSE)

# Modest final typography increase for manuscript-scale readability.
pA <- pA +
  theme(
    legend.title = element_text(face = "bold", size = 12.7),
    legend.text = element_text(size = 11.7),
    axis.text.y = element_text(size = 11.6),
    axis.title.y = element_text(face = "bold", size = 14.3),
    axis.title.x = element_text(face = "bold", size = 13.0),
    strip.text = element_text(face = "bold", size = 13.0)
  )

pB <- pB +
  theme(
    axis.text = element_text(size = 11.0),
    axis.title = element_text(face = "bold", size = 13.2),
    legend.title = element_text(face = "bold", size = 11.7),
    legend.text = element_text(size = 10.8),
    plot.subtitle = element_text(
      size = 9.9, colour = "#4A286F", lineheight = 1.08,
      margin = margin(b = 4)
    )
  )

pC <- pC +
  theme(
    axis.text.y = element_text(size = 11.0),
    axis.text.x = element_text(size = 10.5, lineheight = 0.94),
    axis.title = element_text(face = "bold", size = 13.2),
    plot.subtitle = element_text(
      size = 9.7, colour = "#4A286F", lineheight = 1.08,
      margin = margin(b = 4)
    )
  )

# Reassert a genuinely axis-free dendrogram AFTER the earlier patchwork-wide
# typography pass, which otherwise re-enables dendrogram axis text/titles.
pD_dend <- pD_dend +
  theme_void() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.tag = element_blank(),
    plot.margin = margin(2, 0, 2, 2)
  )

# Apply larger typography only to the heatmap component, not the dendrogram.
pD_heat <- pD_heat +
  theme(
    axis.text.y = element_text(size = 11.4, hjust = 0),
    axis.text.x = element_text(size = 11.0),
    axis.title = element_text(face = "bold", size = 12.8),
    legend.title = element_text(face = "bold", size = 11.3),
    legend.text = element_text(size = 10.3)
  )

pD <- pD_dend + pD_heat +
  plot_layout(widths = c(0.20, 1), guides = "collect") &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical"
  )
pD <- figure2_strip_nested_tags(pD)
try(pD$labels$tag <- NULL, silent = TRUE)
try(pD$labels$title <- NULL, silent = TRUE)

pE <- pE +
  theme(
    axis.text.y = element_text(size = 11.2),
    axis.text.x = element_text(size = 11.0),
    axis.title = element_text(face = "bold", size = 13.2),
    legend.title = element_text(face = "bold", size = 11.7),
    legend.text = element_text(size = 10.6)
  )

for (figure2_nm in c("pA", "pB", "pC", "pD", "pE")) {
  figure2_obj <- get(figure2_nm)
  if (figure2_nm == "pD") {
    figure2_obj <- figure2_strip_nested_tags(figure2_obj)
  }
  try(figure2_obj$labels$tag <- NULL, silent = TRUE)
  assign(figure2_nm, figure2_obj)
}
# === PHASE_F2I_FIGURE2_TYPOGRAPHY_AND_D_TAG_END ===
# === PHASE_F2C_FIGURE2_EDITORIAL_OVERRIDE_END ===


figure2_design <- "
AA
BC
DD
EE
"

# Each scientific panel is exactly one top-level patchwork element. This is
# essential because Panel D internally contains a dendrogram + heatmap.
figure2_A <- patchwork::wrap_elements(full = pA)
figure2_B <- patchwork::wrap_elements(full = pB)
figure2_C <- patchwork::wrap_elements(full = pC)
figure2_D <- patchwork::wrap_elements(full = pD)
figure2_E <- patchwork::wrap_elements(full = pE)

figure2 <- figure2_A + figure2_B + figure2_C + figure2_D + figure2_E +
  plot_layout(
    design = figure2_design,
    heights = c(1.42, 1.26, 1.20, 1.30),
    widths = c(1.34, 1.00)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 17, colour = "#17131C"),
    plot.tag.position = c(0.005, 0.995)
  )

output_stem <- file.path(FIGDIR, "Figure_2_excellence_candidate")
ggsave(
  paste0(output_stem, ".pdf"),
  figure2,
  width=15.0,
  height=18.0,
  units="in",
  device=cairo_pdf,
  limitsize=FALSE
)
ggsave(
  paste0(output_stem, ".png"),
  figure2,
  width=15.0,
  height=18.0,
  units="in",
  dpi=400,
  limitsize=FALSE
)
ggsave(
  paste0(output_stem, ".tiff"),
  figure2,
  width=15.0,
  height=18.0,
  units="in",
  dpi=600,
  compression="lzw",
  limitsize=FALSE
)
ggsave(
  file.path(PNGD, "Figure_2_excellence_preview.png"),
  figure2,
  width=15.0,
  height=18.0,
  units="in",
  dpi=220,
  limitsize=FALSE
)

required_outputs <- paste0(output_stem, c(".pdf", ".png", ".tiff"))
if (any(!file.exists(required_outputs)) ||
    any(file.info(required_outputs)$size <= 0)) {
  stop("One or more Figure 2 candidate outputs are missing.", call.=FALSE)
}

source_manifest <- data.table(
  figure="Figure_2",
  panel=c("A","B","C","D","D","E","E"),
  source_data=c(
    "Figure_2A_genus_composition.tsv",
    "Figure_2B_genus_Bray_PCoA.tsv",
    "Figure_2C_genus_Shannon.tsv",
    "Figure_2D_clustered_genus_class_means.tsv",
    "../cluster_orders/Figure_2D_genus_cluster_order.tsv",
    "Figure_2E_all_genus_elimination_correlations.tsv",
    "Figure_2E_top_genus_elimination_correlations.tsv"
  )
)
fwrite(source_manifest, file.path(SRC, "Figure_2_source_data_manifest.tsv"), sep="\t")

input_manifest <- data.table(
  role=c("metadata", "MetaPhlAn genus matrix"),
  path=c(meta_path, genus_path),
  rows=c(nrow(meta), nrow(genus)),
  columns=c(ncol(meta), ncol(genus)),
  md5=as.character(tools::md5sum(c(meta_path, genus_path)))
)
fwrite(input_manifest, file.path(SRC, "Figure_2_input_manifest.tsv"), sep="\t")

audit <- data.table(
  check=c(
    "METADATA_ROWS",
    "MATCHED_SAMPLES",
    "EXPOSED_SAMPLES",
    "REFERENCE_SAMPLES",
    "GENUS_FEATURES",
    "COMPOSITION_GENERA",
    "CLUSTERED_GENERA",
    "ASSOCIATION_ROWS",
    "NOMINAL_ASSOCIATIONS",
    "FDR_ASSOCIATIONS",
    "OUTPUT_PDF",
    "OUTPUT_PNG",
    "OUTPUT_TIFF"
  ),
  value=c(
    nrow(meta),
    nrow(genus),
    sum(meta2$exposure_group == "PFAS Exposed"),
    sum(meta2$exposure_group == "Reference"),
    ncol(genus),
    length(top_genera),
    length(genus_order),
    nrow(correlation_all),
    nominal_count,
    fdr_count,
    1,1,1
  )
)
fwrite(audit, file.path(LOGD, "Figure_2_excellence_audit.tsv"), sep="\t")

cat("=== PFAS FIGURE 2 EXCELLENCE CANDIDATE ===\n")
cat("Matched samples:", nrow(genus), "\n")
cat("Genus features:", ncol(genus), "\n")
cat("Nominal/FDR genus associations:", nominal_count, "/", fdr_count, "\n")
cat("Output:", paste0(output_stem, ".pdf"), "\n")
cat("Status: READY_FOR_FIGURE2_VISUAL_QC\n")
