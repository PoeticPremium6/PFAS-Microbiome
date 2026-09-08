#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, warn=1)
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Repository root is required.", call.=FALSE)
ROOT <- normalizePath(args[[1]], mustWork=TRUE)

OUT <- file.path(ROOT, "30_figure_program", "figure3_excellence")
FIGDIR <- file.path(OUT, "figures")
SRCDIR <- file.path(OUT, "source_data")
CLUSTERDIR <- file.path(OUT, "cluster_orders")
REPORTDIR <- file.path(ROOT, "30_figure_program", "reports")
FINALDIR <- file.path(ROOT, "29_figure_excellence_rollout", "excellent_figures", "main_text")
for (d in c(OUT, FIGDIR, SRCDIR, CLUSTERDIR, REPORTDIR, FINALDIR)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

required_packages <- c("data.table", "ggplot2", "patchwork", "scales", "cluster")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse=", "), call.=FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

purple <- c(
  deep="#3F007D", dark="#54278F", main="#6A51A3", mid="#807DBA",
  light="#9E9AC8", pale="#DADAEB", palest="#F2F0F7",
  grey_dark="#3F3F3F", grey_mid="#8A8A8A", grey_light="#D9D9D9", missing="#BDBDBD"
)
class_palette <- c("Lower"="#D7C9F2", "Middle"="#9F7AD8", "Higher"="#5B2CA3", "Reference (no k)"="#BDB7C9")
class_labels <- c("Lower"="Lower", "Middle"="Middle", "Higher"="Higher", "Reference (no k)"="Ref.\n(no k)")
class_levels <- names(class_labels)

base_theme <- function(base_size=8.0) {
  theme_bw(base_size=base_size, base_family="sans") +
    theme(
      axis.title=element_text(face="bold", colour=purple["grey_dark"], size=base_size+0.3),
      axis.text=element_text(colour=purple["grey_dark"], size=base_size-0.2),
      legend.title=element_text(face="bold", size=base_size-0.1),
      legend.text=element_text(size=base_size-0.5),
      strip.text=element_text(face="bold", size=base_size-0.3),
      strip.background=element_rect(fill=purple["palest"], colour=purple["grey_mid"], linewidth=0.3),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour="#ECECEC", linewidth=0.22),
      plot.margin=margin(5,7,5,7),
      plot.tag=element_text(face="bold", size=13.0, colour="black"),
      plot.tag.position=c(0.01, 0.995)
    )
}

pick_existing <- function(paths, label) {
  found <- paths[file.exists(paths)]
  found <- found[file.info(found)$size > 0]
  if (!length(found)) stop("No usable ", label, " input among: ", paste(paths, collapse=" | "), call.=FALSE)
  found[[1]]
}

parse_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  less <- grepl("^<", x)
  y <- suppressWarnings(as.numeric(sub("^<", "", x)))
  y[less & is.finite(y)] <- y[less & is.finite(y)] / sqrt(2)
  y
}

wrap_text <- function(x, width=28) {
  vapply(as.character(x), function(value) paste(strwrap(value, width=width), collapse="\n"), character(1))
}

clean_feature_label <- function(x, width=26) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  wrap_text(x, width=width)
}

source_shapes <- c("HUMAnN"=21, "AGORA2"=22)

clean_humann <- function(x) {
  x <- sub("\\|.*$", "", as.character(x))
  x <- sub("^\\s*[A-Za-z0-9_.-]*PWY[-_0-9A-Za-z.]*[: ]+", "", x, ignore.case=TRUE)
  x <- sub("^\\s*PWY[-_0-9A-Za-z.]*[: ]+", "", x, ignore.case=TRUE)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- "Unclassified pathway"
  x
}

clean_subsystem <- function(x) {
  x <- gsub("[_\\.]", " ", as.character(x))
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- "Unclassified subsystem"
  x
}

read_feature_matrix <- function(path, sample_ids, cleaner, aggregate_base=FALSE) {
  dt <- fread(path, data.table=FALSE, check.names=FALSE)
  names(dt) <- make.unique(names(dt), sep="_dup")
  if (ncol(dt) < 2) stop("Matrix has fewer than two columns: ", path, call.=FALSE)
  identifiers <- as.character(dt[[1]])
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  rownames(mat) <- identifiers
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[!is.finite(mat)] <- 0
  row_matches <- sum(rownames(mat) %in% sample_ids)
  column_matches <- sum(colnames(mat) %in% sample_ids)
  if (column_matches > row_matches) mat <- t(mat)
  feature_names <- cleaner(colnames(mat))
  keep <- !grepl("UNMAPPED|UNINTEGRATED|UNGROUPED", feature_names, ignore.case=TRUE)
  if (sum(keep) >= 10) {
    mat <- mat[, keep, drop=FALSE]
    feature_names <- feature_names[keep]
  }
  if (aggregate_base) {
    feature_names <- cleaner(sub("\\|.*$", "", feature_names))
    aggregate_matrix <- rowsum(t(mat), group=feature_names, reorder=FALSE)
    mat <- t(as.matrix(aggregate_matrix))
    colnames(mat) <- rownames(aggregate_matrix)
  } else {
    colnames(mat) <- make.unique(feature_names, sep=" ")
  }
  row_total <- rowSums(mat, na.rm=TRUE)
  row_total[row_total <= 0] <- 1
  if (median(row_total, na.rm=TRUE) > 1.5) mat <- mat / row_total * 100
  mat
}

standardized_mean <- function(df) {
  if (!ncol(df)) return(rep(NA_real_, nrow(df)))
  z <- lapply(df, function(value) {
    value <- parse_num(value)
    s <- sd(value, na.rm=TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(value)))
    (value - mean(value, na.rm=TRUE)) / s
  })
  answer <- rowMeans(as.data.frame(z), na.rm=TRUE)
  answer[!is.finite(answer)] <- NA_real_
  answer
}

zscore_columns <- function(mat) {
  answer <- scale(mat)
  answer[!is.finite(answer)] <- 0
  as.matrix(answer)
}

safe_cor_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 10 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(c(rho=NA_real_, p=NA_real_, n=sum(keep)))
  test <- suppressWarnings(try(cor.test(x[keep], y[keep], method="spearman", exact=FALSE), silent=TRUE))
  if (inherits(test, "try-error")) return(c(rho=NA_real_, p=NA_real_, n=sum(keep)))
  c(rho=unname(test$estimate), p=test$p.value, n=sum(keep))
}

feature_statistics <- function(mat, meta, family) {
  rows <- lapply(colnames(mat), function(feature) {
    x <- mat[, feature]
    use <- meta$exposure_group == "PFAS Exposed" & is.finite(meta$elimination_score) & is.finite(x)
    result <- safe_cor_test(x[use], meta$elimination_score[use])
    lower <- mean(x[use & meta$elimination_class == "Lower"], na.rm=TRUE)
    higher <- mean(x[use & meta$elimination_class == "Higher"], na.rm=TRUE)
    data.table(
      feature=feature,
      family=family,
      rho=as.numeric(result[["rho"]]),
      p_value=as.numeric(result[["p"]]),
      n=as.integer(result[["n"]]),
      prevalence=mean(x[use] > 0, na.rm=TRUE),
      mean_value=mean(x[use], na.rm=TRUE),
      higher_minus_lower=higher-lower
    )
  })
  result <- rbindlist(rows, fill=TRUE)
  result <- result[is.finite(rho)]
  result[, q_value := p.adjust(p_value, method="BH")]
  result[, abs_rho := abs(rho)]
  result
}

class_means <- function(mat, meta, features, family) {
  z <- zscore_columns(mat[, features, drop=FALSE])
  zdt <- as.data.table(z)
  zdt[, sample_id := rownames(mat)]
  zdt <- merge(zdt, meta[, .(sample_id, elimination_class)], by="sample_id", all.x=TRUE)
  long <- melt(zdt, id.vars=c("sample_id", "elimination_class"), variable.name="feature", value.name="z")
  result <- long[, .(mean_z=mean(z, na.rm=TRUE)), by=.(feature, elimination_class)]
  result[, family := family]
  result[, elimination_class := factor(as.character(elimination_class), levels=class_levels)]
  result
}

cluster_order <- function(wide_matrix) {
  matrix_value <- as.matrix(wide_matrix)
  storage.mode(matrix_value) <- "double"
  matrix_value[!is.finite(matrix_value)] <- 0
  if (nrow(matrix_value) <= 1) return(list(order=rownames(matrix_value), tree=NULL))
  correlation <- suppressWarnings(cor(t(matrix_value), method="spearman", use="pairwise.complete.obs"))
  if (any(!is.finite(correlation))) {
    row_z <- t(scale(t(matrix_value)))
    row_z[!is.finite(row_z)] <- 0
    distance <- dist(row_z)
  } else {
    diag(correlation) <- 1
    distance_matrix <- 1 - correlation
    distance_matrix[distance_matrix < 0] <- 0
    diag(distance_matrix) <- 0
    distance <- as.dist(distance_matrix)
  }
  tree <- hclust(distance, method="average")
  list(order=rownames(matrix_value)[tree$order], tree=tree)
}

hclust_segments <- function(tree) {
  if (is.null(tree)) return(data.table())
  n <- length(tree$order)
  leaf_y <- numeric(n)
  leaf_y[tree$order] <- seq_len(n)
  node_y <- numeric(n-1)
  output <- list()
  child_details <- function(child) {
    if (child < 0) return(list(x=0, y=leaf_y[-child]))
    list(x=tree$height[child], y=node_y[child])
  }
  index <- 1
  for (merge_index in seq_len(n-1)) {
    left <- child_details(tree$merge[merge_index, 1])
    right <- child_details(tree$merge[merge_index, 2])
    parent_x <- tree$height[merge_index]
    parent_y <- mean(c(left$y, right$y))
    node_y[merge_index] <- parent_y
    output[[index]] <- data.table(x=left$x, xend=parent_x, y=left$y, yend=left$y); index <- index+1
    output[[index]] <- data.table(x=right$x, xend=parent_x, y=right$y, yend=right$y); index <- index+1
    output[[index]] <- data.table(x=parent_x, xend=parent_x, y=left$y, yend=right$y); index <- index+1
  }
  rbindlist(output)
}

heatmap_composite <- function(heat, family_label, panel_tag, cluster_file) {
  wide <- dcast(heat, feature ~ elimination_class, value.var="mean_z")
  feature_names <- wide$feature
  wide$feature <- NULL
  for (level in class_levels) if (!level %in% names(wide)) wide[[level]] <- 0
  wide <- wide[, ..class_levels]
  wide_matrix <- as.matrix(wide)
  rownames(wide_matrix) <- feature_names

  clustered <- cluster_order(wide_matrix)
  order_labels <- clustered$order

  order_table <- data.table(
    order_index=seq_along(order_labels),
    feature=order_labels,
    family=family_label,
    distance_method="1-Spearman; Euclidean fallback",
    linkage_method="average"
  )
  fwrite(order_table, cluster_file, sep="	")

  heat[, y := match(as.character(feature), order_labels)]
  heat[, class_short := factor(
    as.character(elimination_class),
    levels=class_levels,
    labels=unname(class_labels)
  )]

  label_map <- setNames(
    clean_feature_label(order_labels, width=28),
    seq_along(order_labels)
  )

  heat_plot <- ggplot(heat, aes(x=class_short, y=y, fill=mean_z)) +
    geom_tile(colour="white", linewidth=0.38) +
    scale_y_continuous(
      breaks=seq_along(order_labels),
      labels=label_map,
      limits=c(0.5, length(order_labels)+0.5),
      position="right",
      expand=c(0,0)
    ) +
    scale_fill_gradient2(
      low="#C5B1DA",
      mid="#F0EAF5",
      high=purple["deep"],
      midpoint=0,
      name="Mean z-score"
    ) +
    labs(x="Elimination class", y=family_label) +
    base_theme(9.3) +
    theme(
      panel.grid=element_blank(),
      axis.text.y=element_text(
        size=7.2,
        lineheight=0.92,
        hjust=0
      ),
      axis.text.x=element_text(size=7.8),
      axis.title.y=element_text(
        size=9.2,
        margin=margin(l=7)
      ),
      legend.position="bottom",
      plot.margin=margin(5,16,5,1)
    ) +
    guides(fill=guide_colourbar(
      direction="horizontal",
      title.position="top",
      barwidth=unit(34,"mm"),
      barheight=unit(3.2,"mm")
    ))

  segments <- hclust_segments(clustered$tree)
  dendrogram_plot <- ggplot(segments) +
    geom_segment(
      aes(x=x, xend=xend, y=y, yend=yend),
      linewidth=0.42,
      colour=purple["grey_dark"]
    ) +
    scale_x_reverse(expand=c(0.02,0)) +
    scale_y_continuous(
      limits=c(0.5, length(order_labels)+0.5),
      expand=c(0,0)
    ) +
    labs(tag=panel_tag) +
    theme_void() +
    theme(
      plot.tag=element_text(face="bold", size=13.0, colour="black"),
      plot.tag.position=c(0.02,0.985),
      plot.margin=margin(5,0,5,5)
    )

  wrap_plots(
    dendrogram_plot,
    heat_plot,
    nrow=1,
    widths=c(0.18,1)
  )
}

pcoa_data <- function(mat, meta, max_features, family) {
  variances <- apply(mat, 2, var, na.rm=TRUE)
  variances[!is.finite(variances)] <- 0
  selected <- names(sort(variances, decreasing=TRUE))[seq_len(min(max_features, length(variances)))]
  x <- zscore_columns(mat[, selected, drop=FALSE])
  ordination <- cmdscale(dist(x), eig=TRUE, k=2)
  percent <- round(100 * ordination$eig[1:2] / sum(abs(ordination$eig)), 1)
  result <- data.table(sample_id=rownames(mat), PCoA1=ordination$points[,1], PCoA2=ordination$points[,2])
  result <- merge(result, meta[, .(sample_id, exposure_group, elimination_class, elimination_score)], by="sample_id", all.x=TRUE)
  use <- result$exposure_group == "PFAS Exposed" & is.finite(result$elimination_score)
  axis_test <- safe_cor_test(result$PCoA1[use], result$elimination_score[use])
  list(data=result, percent=percent, rho=axis_test[["rho"]], p=axis_test[["p"]], n=axis_test[["n"]], family=family)
}

META_PATHS <- c(
  file.path(ROOT, "14_figure_assembly/manual_inputs/PFAS_Metadata.csv"),
  file.path(ROOT, "08_microbiome_statistics/00_inputs/phase8_selected_metadata.tsv")
)
HUMANN_PATHS <- c(
  file.path(ROOT, "08_microbiome_statistics/00_inputs/corrected_feature_matrices/humann_pathways_sample_matched.tsv"),
  file.path(ROOT, "05_functional_humann/merged_pathabundance_sample_matched.tsv"),
  file.path(ROOT, "05_functional_humann/humann_pathways_sample_matched.tsv")
)
AGORA_PATHS <- c(
  file.path(ROOT, "08_microbiome_statistics/00_inputs/corrected_feature_matrices/agora2_subsystem_sample_matched.tsv"),
  file.path(ROOT, "07_AGORA2_community_modelling/tables/agora2_subsystem_sample_matched.tsv")
)

meta_path <- pick_existing(META_PATHS, "metadata")
humann_path <- pick_existing(HUMANN_PATHS, "HUMAnN")
agora_path <- pick_existing(AGORA_PATHS, "AGORA2 subsystem")

meta <- fread(meta_path, data.table=TRUE, check.names=FALSE)
setnames(meta, make.unique(names(meta), sep="_dup"))
sample_column <- intersect(c("Sample_Name", "sample_id", "SampleID", "sample", "sample_name"), names(meta))[[1]]
if (is.na(sample_column)) stop("No sample ID column in metadata.", call.=FALSE)
meta[, sample_id := as.character(get(sample_column))]
if ("ExposureGroup" %in% names(meta)) meta[, exposure_group := as.character(ExposureGroup)] else if ("Condition" %in% names(meta)) meta[, exposure_group := as.character(Condition)] else meta[, exposure_group := "Unknown"]
meta[, exposure_group := ifelse(grepl("ref", exposure_group, ignore.case=TRUE), "Reference", "PFAS Exposed")]
meta[, exposure_group := factor(exposure_group, levels=c("PFAS Exposed", "Reference"))]
k_columns <- intersect(names(meta), c("k_LPFOS", "k_PFHxS", "k_PFOA", "k_PFOS_MP1", "k_PFOS_MP26", "k_PFOS_MP345", "k_PFPeS", "k_PFHpS"))
if (!length(k_columns)) k_columns <- grep("^k_", names(meta), value=TRUE)
for (column in k_columns) meta[[column]] <- parse_num(meta[[column]])
meta[, elimination_score := standardized_mean(.SD), .SDcols=k_columns]
meta[, elimination_class := "Reference (no k)"]
exposed_index <- which(meta$exposure_group == "PFAS Exposed" & is.finite(meta$elimination_score))
if (length(exposed_index) >= 6) {
  cut_points <- quantile(meta$elimination_score[exposed_index], probs=c(1/3,2/3), na.rm=TRUE)
  meta[exposed_index, elimination_class := ifelse(elimination_score <= cut_points[[1]], "Lower", ifelse(elimination_score <= cut_points[[2]], "Middle", "Higher"))]
}
meta[, elimination_class := factor(elimination_class, levels=class_levels)]

humann <- read_feature_matrix(humann_path, meta$sample_id, clean_humann, aggregate_base=TRUE)
agora <- read_feature_matrix(agora_path, meta$sample_id, clean_subsystem, aggregate_base=FALSE)
keep_h <- intersect(rownames(humann), meta$sample_id)
keep_a <- intersect(rownames(agora), meta$sample_id)
humann <- humann[keep_h,,drop=FALSE]
agora <- agora[keep_a,,drop=FALSE]
meta_h <- meta[match(keep_h, sample_id)]
meta_a <- meta[match(keep_a, sample_id)]

stats_h <- feature_statistics(humann, meta_h, "HUMAnN pathway")
stats_a <- feature_statistics(agora, meta_a, "AGORA2 subsystem")
all_stats <- rbindlist(list(stats_h, stats_a), fill=TRUE)

select_heat_features <- function(statistics, n=10) {
  statistics <- copy(statistics)
  statistics[, score := 0.65*rank(-abs_rho, ties.method="average")/max(.N,1) + 0.35*rank(-abs(higher_minus_lower), ties.method="average")/max(.N,1)]
  setorder(statistics, score, -abs_rho)
  unique(statistics$feature)[seq_len(min(n, uniqueN(statistics$feature)))]
}
features_h <- select_heat_features(stats_h, 10)
features_a <- select_heat_features(stats_a, 10)
heat_h <- class_means(humann, meta_h, features_h, "HUMAnN pathway")
heat_a <- class_means(agora, meta_a, features_a, "AGORA2 subsystem")
fwrite(heat_h, file.path(SRCDIR, "Figure_3A_HUMAnN_clustered_class_means.tsv"), sep="\t")
fwrite(heat_a, file.path(SRCDIR, "Figure_3B_AGORA2_clustered_class_means.tsv"), sep="\t")

pA <- heatmap_composite(heat_h, "HUMAnN pathway", "A", file.path(CLUSTERDIR, "Figure_3A_HUMAnN_cluster_order.tsv"))
pB <- heatmap_composite(heat_a, "AGORA2 subsystem", "B", file.path(CLUSTERDIR, "Figure_3B_AGORA2_cluster_order.tsv"))

pcoa_h <- pcoa_data(humann, meta_h, 700, "HUMAnN")
pcoa_a <- pcoa_data(agora, meta_a, 300, "AGORA2")
fwrite(pcoa_h$data, file.path(SRCDIR, "Figure_3C_HUMAnN_PCoA.tsv"), sep="\t")
fwrite(pcoa_a$data, file.path(SRCDIR, "Figure_3D_AGORA2_PCoA.tsv"), sep="\t")

# >>> PHASE_F3A2_FIGURE3_EDITORIAL_REFINEMENT_BEGIN
# Display-only helpers for Figure 3 editorial refinement.
figure3_theme_readable <- ggplot2::theme(
  text = ggplot2::element_text(size = 12),
  axis.title = ggplot2::element_text(size = 11.5, face = "bold"),
  axis.text = ggplot2::element_text(size = 9.5),
  legend.title = ggplot2::element_text(size = 10.5, face = "bold"),
  legend.text = ggplot2::element_text(size = 9.5),
  strip.text = ggplot2::element_text(size = 10.5, face = "bold"),
  plot.subtitle = ggplot2::element_text(
    size = 8.5,
    colour = "#4A286F",
    hjust = 0,
    lineheight = 1.05,
    margin = ggplot2::margin(b = 4)
  )
)

figure3_heatmap_editorial <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  p +
    figure3_theme_readable +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      axis.text.y = ggplot2::element_text(size = 9.0, lineheight = 0.90),
      axis.text.y.right = ggplot2::element_text(size = 9.0, lineheight = 0.90),
      axis.text.x = ggplot2::element_text(size = 9.2),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(52, "pt"),
        barwidth = grid::unit(9, "pt")
      ),
      colour = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(52, "pt"),
        barwidth = grid::unit(9, "pt")
      )
    )
}

figure3_move_stat_annotation <- function(p) {
  if (!inherits(p, "ggplot") || length(p$layers) == 0L) return(p)
  keep <- rep(TRUE, length(p$layers))
  labels <- character()
  for (i in seq_along(p$layers)) {
    layer <- p$layers[[i]]
    is_text <- inherits(layer$geom, "GeomLabel") || inherits(layer$geom, "GeomText")
    dat <- layer$data
    if (!is_text || !is.data.frame(dat) || !("label" %in% names(dat))) next
    lab <- as.character(dat$label)
    stat_like <- grepl("rho\\(PC1|elimination|p\\s*=", lab, ignore.case = TRUE)
    if (any(stat_like, na.rm = TRUE)) {
      labels <- c(labels, lab[stat_like])
      keep[i] <- FALSE
    }
  }
  if (length(labels) > 0L) {
    p$layers <- p$layers[keep]
    p <- p + ggplot2::labs(subtitle = paste(unique(labels), collapse = "\n"))
  }
  p
}

figure3_ordination_editorial <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  p <- figure3_move_stat_annotation(p)
  p +
    figure3_theme_readable +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 10.2, face = "bold"),
      legend.text = ggplot2::element_text(size = 9.2),
      axis.title = ggplot2::element_text(size = 11.8, face = "bold"),
      axis.text = ggplot2::element_text(size = 9.8),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
}

figure3_general_editorial <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  p + figure3_theme_readable +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 9.5),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}
# <<< PHASE_F3A2_FIGURE3_EDITORIAL_REFINEMENT_END


pC <- ggplot(pcoa_h$data, aes(PCoA1, PCoA2, fill=elimination_class, shape=exposure_group)) +
  geom_point(size=2.5, colour=purple["grey_dark"], stroke=0.32, alpha=0.9) +
  scale_fill_manual(values=class_palette, labels=class_labels, drop=FALSE) +
  scale_shape_manual(values=c("PFAS Exposed"=21, "Reference"=24), drop=FALSE) +
  annotate(
    "label",
    x=-Inf,
    y=Inf,
    hjust=-0.08,
    vjust=1.08,
    size=2.1,
    linewidth=0.18,
    fill="white",
    colour=purple["deep"],
    label=sprintf(
      "rho(PC1, elimination) = %.2f\np = %.3g; n = %d",
      pcoa_h$rho,
      pcoa_h$p,
      pcoa_h$n
    )
  ) +
  labs(
    tag="C",
    x=sprintf("HUMAnN PCoA1 (%.1f%%)", pcoa_h$percent[[1]]),
    y=sprintf("PCoA2 (%.1f%%)", pcoa_h$percent[[2]]),
    fill="Elimination class",
    shape="Group"
  ) +
  base_theme(9.3) +
  theme(
    legend.position="bottom",
    legend.box="vertical",
    legend.justification="left",
    legend.key.height=unit(3.2,"mm"),
    legend.key.width=unit(4.0,"mm"),
    legend.margin=margin(1,0,0,0),
    plot.margin=margin(5,5,3,5)
  ) +
  guides(
    fill=guide_legend(nrow=2, byrow=TRUE, title.position="top"),
    shape=guide_legend(nrow=1, title.position="top")
  )


pC <- figure3_ordination_editorial(pC)
# >>> PHASE_F3B2_FIGURE3_PANEL_C_SIDE_LEGEND_BEGIN
# Focused Panel C legend placement override.
# Existing scales, colours, shapes and data are preserved.
pC <- pC +
  ggplot2::theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.justification = c(0.5, 0.5),
    legend.title = ggplot2::element_text(size = 10.0, face = "bold"),
    legend.text = ggplot2::element_text(size = 9.0),
    legend.key.height = grid::unit(10, "pt"),
    legend.key.width = grid::unit(10, "pt"),
    legend.spacing.y = grid::unit(2, "pt"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 2, 0, 5),
    plot.margin = ggplot2::margin(6, 7, 6, 6)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(order = 1, ncol = 1),
    shape = ggplot2::guide_legend(order = 2, ncol = 1)
  )
# <<< PHASE_F3B2_FIGURE3_PANEL_C_SIDE_LEGEND_END
pD <- ggplot(pcoa_a$data, aes(PCoA1, PCoA2, fill=elimination_class, shape=exposure_group)) +
  geom_point(size=2.5, colour=purple["grey_dark"], stroke=0.32, alpha=0.9) +
  scale_fill_manual(values=class_palette, labels=class_labels, drop=FALSE) +
  scale_shape_manual(values=c("PFAS Exposed"=21, "Reference"=24), drop=FALSE) +
  annotate(
    "label",
    x=-Inf,
    y=Inf,
    hjust=-0.08,
    vjust=1.08,
    size=2.1,
    linewidth=0.18,
    fill="white",
    colour=purple["deep"],
    label=sprintf(
      "rho(PC1, elimination) = %.2f\np = %.3g; n = %d",
      pcoa_a$rho,
      pcoa_a$p,
      pcoa_a$n
    )
  ) +
  labs(
    tag="D",
    x=sprintf("AGORA2 PCoA1 (%.1f%%)", pcoa_a$percent[[1]]),
    y=sprintf("PCoA2 (%.1f%%)", pcoa_a$percent[[2]]),
    fill="Elimination class",
    shape="Group"
  ) +
  base_theme(9.3) +
  theme(
    legend.position="none",
    plot.margin=margin(5,5,3,5)
  )

correlation_distribution <- all_stats[, .(rho, p_value, q_value, family)]
fwrite(correlation_distribution, file.path(SRCDIR, "Figure_3E_correlation_distribution.tsv"), sep="\t")
family_summary <- all_stats[, .(
  feature_count=.N,
  median_rho=median(rho, na.rm=TRUE),
  nominal_count=sum(p_value < 0.05, na.rm=TRUE),
  fdr_count=sum(q_value < 0.05, na.rm=TRUE)
), by=family]
y_annotation <- max(abs(all_stats$rho), na.rm=TRUE) * 1.08

pD <- figure3_ordination_editorial(pD)
pE <- ggplot(all_stats, aes(x=family, y=rho, fill=family)) +
  geom_hline(yintercept=0, linewidth=0.4, colour=purple["grey_mid"]) +
  geom_violin(width=0.72, colour=purple["grey_dark"], linewidth=0.3, alpha=0.9) +
  geom_boxplot(width=0.15, outlier.shape=NA, fill="white", linewidth=0.35) +
  geom_text(data=family_summary, aes(x=family, y=y_annotation, label=sprintf("n=%d\nmedian=%.2f\nP<0.05: %d; FDR: %d", feature_count, median_rho, nominal_count, fdr_count)), inherit.aes=FALSE, size=2.0, vjust=0) +
  scale_fill_manual(values=c("HUMAnN pathway"=unname(purple["light"]), "AGORA2 subsystem"=unname(purple["dark"]))) +
  scale_y_continuous(limits=c(-y_annotation*1.05, y_annotation*1.34)) +
  labs(tag="E", x=NULL, y="Spearman rho") +
  base_theme(9.3) + theme(legend.position="none", axis.text.x=element_text(size=7.8))

# Clustered functional evidence matrix.

pE <- pE + theme(
  plot.tag.position=c(0.02, 0.985),
  plot.margin=margin(6, 10, 6, 14)
)

selected_h <- stats_h[order(-abs_rho)][seq_len(min(7,.N))]
selected_a <- stats_a[order(-abs_rho)][seq_len(min(7,.N))]
selected <- rbindlist(list(selected_h, selected_a), fill=TRUE)
selected[, source_type := fifelse(
  family == "HUMAnN pathway",
  "HUMAnN",
  "AGORA2"
)]
selected[, feature_label := clean_feature_label(feature, width=27)]
selected[, feature_key := sprintf("feature_%02d", seq_len(.N))]
selected[, positive_rho := pmax(rho, 0)]
selected[, negative_rho := pmax(-rho, 0)]
selected[, association_strength := -log10(pmax(p_value, 1e-300))]
selected[, nominal_support := as.numeric(p_value < 0.05)]
selected[, fdr_support := as.numeric(q_value < 0.05)]
selected[, shift_magnitude := abs(higher_minus_lower)]

scale01 <- function(value) {
  value <- as.numeric(value)
  keep <- is.finite(value)
  output <- rep(NA_real_, length(value))
  if (!any(keep)) return(output)
  limits <- range(value[keep], na.rm=TRUE)
  if (diff(limits) == 0) { output[keep] <- 0.5; return(output) }
  output[keep] <- (value[keep]-limits[[1]])/diff(limits)
  output
}

cluster_input <- data.frame(
  positive_rho=selected$positive_rho,
  negative_rho=selected$negative_rho,
  association_strength=selected$association_strength,
  nominal_support=selected$nominal_support,
  fdr_support=selected$fdr_support,
  prevalence=selected$prevalence,
  mean_value=selected$mean_value,
  shift_magnitude=selected$shift_magnitude,
  check.names=FALSE
)
rownames(cluster_input) <- selected$feature_key
for (column in names(cluster_input)) {
  value <- as.numeric(cluster_input[[column]])
  value[!is.finite(value)] <- if (any(is.finite(value))) median(value[is.finite(value)], na.rm=TRUE) else 0
  cluster_input[[column]] <- value
}
binary_columns <- c("nominal_support", "fdr_support")
variable_binary <- binary_columns[vapply(cluster_input[,binary_columns,drop=FALSE], function(value) all(c(0,1) %in% unique(value)), logical(1))]
asymmetric <- match(variable_binary, names(cluster_input))
if (length(asymmetric)) evidence_distance <- cluster::daisy(cluster_input, metric="gower", type=list(asymm=asymmetric)) else evidence_distance <- cluster::daisy(cluster_input, metric="gower")
evidence_tree <- hclust(evidence_distance, method="average")
evidence_order <- rownames(cluster_input)[evidence_tree$order]
cluster_id <- cutree(evidence_tree, k=min(4,nrow(cluster_input)))
cluster_export <- data.table(feature_id=rownames(cluster_input), cluster_id=as.integer(cluster_id), order_index=match(rownames(cluster_input), evidence_order), distance_method="Gower", linkage_method="average")
setorder(cluster_export, order_index)
fwrite(cluster_export, file.path(CLUSTERDIR, "Figure_3F_functional_evidence_cluster_order.tsv"), sep="\t")

wide_evidence <- data.table(
  feature_id=rownames(cluster_input),
  `Positive rho`=scale01(selected$positive_rho),
  `Negative rho`=scale01(selected$negative_rho),
  `Association strength`=scale01(selected$association_strength),
  `Nominal support`=selected$nominal_support,
  `FDR support`=selected$fdr_support,
  `Prevalence`=scale01(selected$prevalence),
  `Mean value`=scale01(selected$mean_value),
  `Higher–Lower shift`=scale01(selected$shift_magnitude)
)
metric_order <- c(
  "Positive rho",
  "Negative rho",
  "Association strength",
  "Nominal support",
  "FDR support",
  "Prevalence",
  "Mean value",
  "Higher–Lower shift"
)
feature_metadata <- selected[, .(
  feature_id=feature_key,
  display_label=feature_label,
  source_type
)]
evidence_long <- melt(
  wide_evidence,
  id.vars="feature_id",
  variable.name="metric",
  value.name="scaled_evidence"
)
evidence_long <- merge(
  evidence_long,
  feature_metadata,
  by="feature_id",
  all.x=TRUE,
  sort=FALSE
)
evidence_long[, metric := factor(metric, levels=metric_order)]
evidence_long[, feature_id := factor(feature_id, levels=evidence_order)]
evidence_long[, source_type := factor(
  source_type,
  levels=c("HUMAnN","AGORA2")
)]
evidence_long[, binary_symbol := fifelse(
  as.character(metric) %in% c("Nominal support","FDR support") &
    scaled_evidence >= 0.5,
  "●",
  ""
)]
fwrite(
  evidence_long,
  file.path(SRCDIR, "Figure_3F_clustered_functional_evidence.tsv"),
  sep="	"
)

segments_f <- hclust_segments(evidence_tree)

source_annotation <- unique(
  evidence_long[, .(feature_id, display_label, source_type)]
)
source_annotation[, feature_id := factor(
  as.character(feature_id),
  levels=evidence_order
)]
source_annotation[, source_type := factor(
  as.character(source_type),
  levels=c("HUMAnN","AGORA2")
)]
feature_label_map <- setNames(
  source_annotation$display_label,
  as.character(source_annotation$feature_id)
)


pE <- figure3_general_editorial(pE)

pE <- figure3_heatmap_editorial(pE)
pF_dend <- ggplot(segments_f) +
  geom_segment(
    aes(x=x, xend=xend, y=y, yend=yend),
    linewidth=0.42,
    colour=purple["grey_dark"]
  ) +
  scale_x_reverse(expand=c(0.02,0)) +
  scale_y_continuous(
    limits=c(0.5, length(evidence_order)+0.5),
    expand=c(0,0)
  ) +
  labs(tag="F") +
  theme_void() +
  theme(
    plot.tag=element_text(face="bold", size=13.0, colour="black"),
    plot.tag.position=c(0.02,0.985),
    plot.margin=margin(5,0,5,5)
  )

pF_source <- ggplot(
  source_annotation,
  aes(x=1, y=feature_id, shape=source_type)
) +
  geom_point(
    size=2.8,
    stroke=0.45,
    fill="white",
    colour=purple["deep"]
  ) +
  scale_shape_manual(
    values=source_shapes,
    name="Feature source"
  ) +
  scale_x_continuous(
    limits=c(0.75,1.25),
    expand=c(0,0)
  ) +
  scale_y_discrete(drop=FALSE) +
  theme_void() +
  theme(
    legend.position="bottom",
    legend.title=element_text(face="bold", size=8.5),
    legend.text=element_text(size=8.0),
    legend.key.width=unit(4.5,"mm"),
    plot.margin=margin(5,1,5,1)
  )

pF_heat <- ggplot(
  evidence_long,
  aes(x=metric, y=feature_id, fill=scaled_evidence)
) +
  geom_tile(colour="white", linewidth=0.38) +
  geom_text(
    aes(label=binary_symbol),
    colour="white",
    size=2.2
  ) +
  scale_y_discrete(
    labels=feature_label_map,
    position="right",
    drop=FALSE
  ) +
  scale_fill_gradientn(
    colours=c(
      "#E9E1F1",
      "#CCBBDD",
      purple["mid"],
      purple["deep"]
    ),
    limits=c(0,1),
    na.value=purple["missing"],
    name="Scaled evidence"
  ) +
  labs(x="Evidence metric", y="Functional feature") +
  base_theme(9.3) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(
      angle=25,
      hjust=1,
      size=7.8
    ),
    axis.text.y=element_text(
      size=7.2,
      lineheight=0.92,
      hjust=0
    ),
    axis.title.y=element_text(
      size=9.2,
      margin=margin(l=7)
    ),
    legend.position="bottom",
    plot.margin=margin(5,16,5,1)
  ) +
  guides(fill=guide_colourbar(
    direction="horizontal",
    title.position="top",
    barwidth=unit(40,"mm"),
    barheight=unit(3.2,"mm")
  ))

pF <- wrap_plots(
  pF_dend,
  pF_source,
  pF_heat,
  nrow=1,
  widths=c(0.16,0.07,1),
  guides="collect"
) & theme(legend.position="bottom")



# --- Figure 3 closeout: clean ordination legends for Panels C and D ---
class_fills <- c(
  "Lower" = unname(purple["mid"]),
  "Middle" = unname(purple["light"]),
  "Higher" = unname(purple["deep"]),
  "Ref.\n(no k)" = unname(purple["pale"]),
  "Ref. (no k)" = unname(purple["pale"]),
  "Reference (no k)" = unname(purple["pale"])
)

group_shapes <- c(
  "PFAS Exposed" = 21,
  "Reference" = 24
)

pC <- pC +
  scale_fill_manual(
    values = class_fills,
    name = "Elimination class"
  ) +
  scale_shape_manual(
    values = group_shapes,
    name = "Group"
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 3.2,
        colour = unname(purple["grey_dark"]),
        stroke = 0.35,
        alpha = 1
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "white",
        size = 3.2,
        colour = unname(purple["grey_dark"]),
        stroke = 0.35,
        alpha = 1
      )
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8.4),
    legend.text = element_text(size = 7.6)
  )

pD <- pD +
  scale_fill_manual(values = class_fills, name = "Elimination class") +
  scale_shape_manual(values = group_shapes, name = "Group") +
  theme(legend.position = "none")


row1 <- wrap_plots(
  pA,
  pB,
  nrow=1,
  widths=c(1,1)
)

row2 <- wrap_plots(
  pC,
  pD,
  pE,
  nrow=1,
  widths=c(1.08,1.08,0.92)
)

figure3 <- wrap_plots(
  row1,
  row2,
  pF,
  ncol=1,
  heights=c(1.38,1.00,1.72)
)

output_stem <- file.path(FIGDIR, "Figure_3_excellence_v1")
ggsave(
  paste0(output_stem,".pdf"),
  figure3,
  width=14.6,
  height=14.4,
  units="in",
  device=cairo_pdf,
  limitsize=FALSE
)
ggsave(
  paste0(output_stem,".png"),
  figure3,
  width=14.6,
  height=14.4,
  units="in",
  dpi=400,
  limitsize=FALSE
)
ggsave(
  paste0(output_stem,".tiff"),
  figure3,
  width=14.6,
  height=14.4,
  units="in",
  dpi=600,
  compression="lzw",
  limitsize=FALSE
)

# Promote only after all formats exist.
for (extension in c("pdf","png","tiff")) {
  source <- paste0(output_stem,".",extension)
  if (!file.exists(source) || file.info(source)$size <= 0) stop("Missing generated output: ", source, call.=FALSE)
  file.copy(source, file.path(FINALDIR, paste0("Figure_3.",extension)), overwrite=TRUE)
}

fwrite(all_stats, file.path(SRCDIR, "Figure_3_all_functional_association_statistics.tsv"), sep="\t")
input_manifest <- data.table(role=c("metadata","HUMAnN","AGORA2 subsystem"), path=c(meta_path,humann_path,agora_path), rows=c(nrow(meta),nrow(humann),nrow(agora)), columns=c(ncol(meta),ncol(humann),ncol(agora)), md5=as.character(tools::md5sum(c(meta_path,humann_path,agora_path))))
fwrite(input_manifest, file.path(SRCDIR, "Figure_3_input_manifest.tsv"), sep="\t")

audit <- data.table(
  check=c("METADATA_ROWS","HUMANN_SAMPLES","HUMANN_FEATURES","AGORA_SAMPLES","AGORA_FEATURES","HUMANN_HEAT_FEATURES","AGORA_HEAT_FEATURES","FUNCTIONAL_ASSOCIATIONS","PANEL_F_FEATURES","PANEL_F_METRICS","OUTPUT_PDF","OUTPUT_PNG","OUTPUT_TIFF"),
  value=c(nrow(meta),nrow(humann),ncol(humann),nrow(agora),ncol(agora),length(features_h),length(features_a),nrow(all_stats),nrow(cluster_input),length(metric_order),1,1,1)
)
fwrite(audit, file.path(REPORTDIR, "phase30b_figure3_excellence_audit.tsv"), sep="\t")
cat("=== PFAS PHASE 30B FIGURE 3 EXCELLENCE ===\n")
cat("HUMAnN samples/features:", nrow(humann), ncol(humann), "\n")
cat("AGORA2 samples/features:", nrow(agora), ncol(agora), "\n")
cat("Functional association rows:", nrow(all_stats), "\n")
cat("Panel F clustered features:", nrow(cluster_input), "\n")
cat("Output:", paste0(output_stem,".pdf"), "\n")
cat("Status: READY_FOR_FIGURE3_VISUAL_QC\n")


pF_heat <- figure3_general_editorial(pF_heat)

pF_heat <- figure3_heatmap_editorial(pF_heat)
