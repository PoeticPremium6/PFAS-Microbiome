# FIGURE4_PANELA_LEGEND_HOTFIX_V1_4
# FIGURE4_PANELB_FEATURE_LEVEL_DISPLAY_V1_3
# FIGURE4_PANELA_PANELB_HOTFIX_V1_2
# FIGURE4_READABILITY_EVIDENCE_PATCH_V1_1
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Repository root is required.", call.=FALSE)

ROOT <- normalizePath(args[[1]], mustWork=TRUE)
BUILD <- file.path(ROOT, "submission", "figure_builds", "Figure_4")
SRC <- file.path(BUILD, "source_data")
FIGDIR <- file.path(BUILD, "figures")
PREVIEWDIR <- file.path(BUILD, "previews")
REPORTDIR <- file.path(BUILD, "reports")

for (directory in c(FIGDIR, PREVIEWDIR, REPORTDIR)) {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
}

input_files <- c(
  A=file.path(SRC, "Figure_4A_low_p_enrichment_by_method_family_v6.tsv"),
  B=file.path(SRC, "Figure_4B_prevalence_association_strength_v6.tsv"),
  C=file.path(SRC, "Figure_4C_Phase11_permutation_null_diagnostics_v6.tsv"),
  D=file.path(SRC, "Figure_4D_ranked_exploratory_features_v6.tsv"),
  E=file.path(SRC, "Figure_4E_candidate_support_heatmap_v6.tsv")
)

missing <- input_files[!file.exists(input_files) | file.info(input_files)$size <= 0]
if (length(missing)) {
  stop(
    paste("Missing Figure 4 source data:", paste(missing, collapse=", ")),
    call.=FALSE
  )
}

A <- fread(input_files["A"], showProgress=FALSE)
B <- fread(input_files["B"], showProgress=FALSE)
C <- fread(input_files["C"], showProgress=FALSE)
D <- fread(input_files["D"], showProgress=FALSE)
E <- fread(input_files["E"], showProgress=FALSE)

required <- list(
  A=c("method","feature_family","fraction_raw_p05","n_q05"),
  B=c("feature_family","prevalence","neglog10p"),
  C=c("feature_family","observed","null_mean","null_q95"),
  D=c("feature_clean","feature_family","method","p_plot","neglog10p"),
  E=c("feature_clean","method","score","tile_label")
)

for (panel in names(required)) {
  object <- get(panel)
  absent <- setdiff(required[[panel]], names(object))
  if (length(absent)) {
    stop(
      sprintf("Panel %s source data is missing columns: %s",
              panel, paste(absent, collapse=", ")),
      call.=FALSE
    )
  }
}

purple <- c(
  deepest="#3B0764",
  deep="#4A148C",
  main="#5E2A9E",
  mid="#7B4EB2",
  soft="#A58AC8",
  light="#C8B7DD",
  pale="#E4DAEF",
  faint="#F6F2F9",
  ink="#2D2433",
  grey="#746C78",
  grid="#E9E3ED"
)

method_order <- c("ANCOM-BC2", "ALDEx2", "CLR sensitivity")
method_cols <- c(
  "ANCOM-BC2"=unname(purple["deepest"]),
  "ALDEx2"=unname(purple["mid"]),
  "CLR sensitivity"=unname(purple["soft"])
)
method_shapes <- c(
  "ANCOM-BC2"=21,
  "ALDEx2"=24,
  "CLR sensitivity"=22
)

family_cols <- c(
  "MetaPhlAn species"="#3B0764",
  "MetaPhlAn genus"="#54278F",
  "MAGs"="#6A3D9A",
  "HUMAnN pathways"="#7B4EB2",
  "AGORA2 reactions"="#8E63B8",
  "AGORA2 subsystems"="#AA8CC8",
  "AGORA2 exchanges"="#C7B5DE",
  "Other"="#8A818E"
)

theme_excellence <- function(base_size=11.2) {
  theme_bw(base_size=base_size, base_family="sans") +
    theme(
      panel.grid.major=element_line(colour=purple["grid"], linewidth=0.25),
      panel.grid.minor=element_blank(),
      panel.border=element_rect(colour="#6E6673", linewidth=0.38),
      axis.title=element_text(
        face="bold", colour=purple["ink"], size=base_size+0.8
      ),
      axis.text=element_text(
        colour=purple["ink"], size=base_size
      ),
      strip.background=element_rect(
        fill=purple["pale"], colour="#807587", linewidth=0.35
      ),
      strip.text=element_text(
        face="bold", colour=purple["ink"], size=base_size
      ),
      legend.title=element_text(
        face="bold", colour=purple["ink"], size=base_size
      ),
      legend.text=element_text(
        colour=purple["ink"], size=base_size-0.8
      ),
      plot.tag=element_text(
        face="bold", colour=purple["ink"], size=15
      ),
      plot.tag.position=c(0.01, 0.99),
      plot.margin=margin(4, 5, 4, 5)
    )
}

wrap_text <- function(x, width=34) {
  vapply(
    as.character(x),
    function(value) paste(strwrap(value, width=width), collapse="\n"),
    character(1)
  )
}

humanise_feature <- function(x) {
  x <- as.character(x)

  x <- gsub(
    "^Bacteroides\\.s\\s+",
    "",
    x
  )

  x <- gsub(
    "^MAG bin ([A-Za-z0-9]+) bin\\.?",
    "MAG \\1 bin ",
    x
  )

  x <- gsub(
    "Exchange: 2ddglcn [EX_2ddglcn(e)]",
    "2-dehydro-D-gluconate exchange",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "Exchange: butyrate [EX_but(e)]",
    "Butyrate exchange",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "Model reaction [EXarbt(e)]",
    "Arabitol exchange",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "Model reaction [EXsbtD(e)]",
    "D-sorbitol exchange",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "Exchange: thr L [EX_thr_L(e)]",
    "L-threonine exchange",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "Model reaction [EXmnl(e)]",
    "Mannitol exchange",
    x,
    fixed=TRUE
  )

  x
}

# Shared feature-family order: observed permutation count, then nominal enrichment.
family_stats <- A[, .(
  max_fraction=max(fraction_raw_p05, na.rm=TRUE),
  q_total=sum(n_q05, na.rm=TRUE)
), by=feature_family]

family_stats <- merge(
  family_stats,
  C[, .(feature_family, observed)],
  by="feature_family",
  all=TRUE
)
family_stats[!is.finite(observed), observed := 0]
family_stats[!is.finite(max_fraction), max_fraction := 0]
family_stats[!is.finite(q_total), q_total := 0]
setorder(family_stats, -observed, -max_fraction, -q_total)

family_order <- family_stats$feature_family
family_levels_y <- rev(family_order)

# -----------------------------
# A. Nominal enrichment by family and method
# -----------------------------
A[, method := factor(
  as.character(method),
  levels=method_order
)]

A[, feature_family := factor(
  as.character(feature_family),
  levels=family_levels_y
)]

pA <- ggplot(A, aes(y=feature_family)) +
  geom_vline(
    xintercept=0.05,
    linetype="dashed",
    colour=purple["grey"],
    linewidth=0.48
  ) +
  geom_segment(
    aes(
      x=0,
      xend=fraction_raw_p05,
      yend=feature_family,
      colour=method,
      group=method
    ),
    position=position_dodge(width=0.66),
    linewidth=1.05,
    alpha=0.82,
    lineend="round"
  ) +
  geom_point(
    aes(
      x=fraction_raw_p05,
      fill=method,
      shape=method,
      group=method
    ),
    position=position_dodge(width=0.66),
    size=4.2,
    colour=purple["ink"],
    stroke=0.40
  ) +
  scale_colour_manual(values=method_cols, drop=TRUE) +
  scale_fill_manual(values=method_cols, drop=TRUE) +
  scale_shape_manual(values=method_shapes, drop=TRUE) +
  scale_x_continuous(
    labels=percent_format(accuracy=1),
    expand=expansion(mult=c(0, 0.16))
  ) +
  labs(
    x="Tests with nominal p < 0.05",
    y="Feature family",
    tag="A"
  ) +
  theme_excellence(11.1) +
  theme(
    legend.position="none",
    axis.text.y=element_text(size=10.4)
  )

# -----------------------------
# B. Prevalence–association landscape
# -----------------------------

# FIGURE4_READABILITY_EVIDENCE_PATCH_V1_1: Panel A legend-only cleanup.
figure4_trim_panelA_shape_legend <- function(plot_object) {
  if (!inherits(plot_object, "ggplot")) return(plot_object)
  for (idx in seq_along(plot_object$scales$scales)) {
    sc <- plot_object$scales$scales[[idx]]
    aesthetics <- tryCatch(sc$aesthetics, error=function(e) character())
    if (!("shape" %in% aesthetics)) next
    breaks_now <- tryCatch(sc$get_breaks(), error=function(e) character())
    keep <- intersect(c("ALDEx2", "CLR sensitivity"), as.character(breaks_now))
    if (length(keep)) {
      plot_object$scales$scales[[idx]]$breaks <- keep
      plot_object$scales$scales[[idx]]$name <- "Method"
    }
  }
  plot_object
}
pA <- figure4_trim_panelA_shape_legend(pA)

# FIGURE4_PANELA_PANELB_HOTFIX_V1_2: Panel A legend cleanup.
# Keep only the two actually plotted method symbols.
pA <- pA +
  scale_shape_manual(
    name="Method",
    values=c("ALDEx2"=16, "CLR sensitivity"=17),
    breaks=c("ALDEx2", "CLR sensitivity"),
    drop=TRUE
  ) +
  guides(shape=guide_legend(order=1))

# FIGURE4_PANELA_LEGEND_HOTFIX_V1_4: only used Panel A method symbols.
pA <- pA +
  scale_shape_manual(
    name="Method",
    values=c("ALDEx2"=16, "CLR sensitivity"=17),
    breaks=c("ALDEx2", "CLR sensitivity"),
    labels=c("ALDEx2", "CLR sensitivity"),
    drop=TRUE
  ) +
  guides(
    shape=guide_legend(order=1),
    colour="none",
    fill="none",
    linetype="none",
    alpha="none"
  )

B[, feature_family := factor(as.character(feature_family), levels=family_order)]
if (!"q_plot" %in% names(B)) B[, q_plot := NA_real_]
if (!"p_plot" %in% names(B)) B[, p_plot := 10^(-neglog10p)]

B[, evidence := fifelse(
  is.finite(q_plot) & q_plot < 0.05,
  "FDR < 0.05",
  fifelse(p_plot < 0.05, "Nominal p < 0.05", "Other")
)]
B[, evidence := factor(
  evidence,
  levels=c("Other", "Nominal p < 0.05", "FDR < 0.05")
)]

B_summary <- B[, .(
  tested_n=uniqueN(feature_clean),
  fdr_n=uniqueN(feature_clean[is.finite(q_plot) & q_plot < 0.05]),
  y=max(neglog10p, na.rm=TRUE)
), by=feature_family]
B_summary[, x := 0.02]
B_summary[, label := paste0("n=", tested_n, "\nFDR=", fdr_n)]

evidence_cols <- c(
  "Other"=unname(purple["pale"]),
  "Nominal p < 0.05"=unname(purple["mid"]),
  "FDR < 0.05"="#26004D"
)
evidence_sizes <- c(
  "Other"=0.70,
  "Nominal p < 0.05"=1.05,
  "FDR < 0.05"=1.55
)
evidence_alpha <- c(
  "Other"=0.22,
  "Nominal p < 0.05"=0.64,
  "FDR < 0.05"=0.94
)


# FIGURE4_READABILITY_EVIDENCE_PATCH_V1_1: stronger evidence-state separation.
# Statistical classifications are unchanged.
evidence_cols["Nominal p < 0.05"] <- "#8A63B8"
evidence_cols["FDR < 0.05"] <- "#26004D"

pB <- ggplot(B, aes(x=prevalence, y=neglog10p)) +
  geom_hline(
    yintercept=-log10(0.05),
    linetype="dashed",
    colour=purple["grey"],
    linewidth=0.42
  ) +
  geom_point(
    aes(colour=evidence, size=evidence, alpha=evidence),
    stroke=0
  ) +
  geom_text(size = 3.4, 
    data=B_summary,
    aes(x=x, y=y, label=label),
    inherit.aes=FALSE,
    hjust=0,
    vjust=1.0,
    size=2.8,
    fontface="bold",
    colour=purple["ink"]
  ) +
  facet_wrap(
    ~feature_family,
    ncol=4,
    scales="free_y",
    labeller=labeller(
      feature_family=function(x) wrap_text(x, width=14)
    )
  ) +
  scale_x_continuous(
    labels=percent_format(accuracy=1),
    limits=c(0, 1),
    breaks=c(0, 0.25, 0.50, 0.75, 1)
  ) +
  scale_colour_manual(values=evidence_cols, drop=FALSE) +
  scale_size_manual(values=evidence_sizes, drop=FALSE) +
  scale_alpha_manual(values=evidence_alpha, drop=FALSE) +
  labs(
    x="Feature prevalence",
    y=expression(-log[10]~"raw p"),
    colour="Association evidence",
    tag="B"
  ) +
  theme_excellence(10.0) +
  theme(
    strip.text=element_text(size=9.7, lineheight=0.92),
    axis.text.x=element_text(size=8.1),
    axis.text.y=element_text(size=8.3),
    panel.spacing=unit(0.58, "lines"),
    legend.position="bottom"
  ) +
  guides(
    size="none",
    alpha="none",
    colour=guide_legend(
      title.position="top",
      nrow=1,
      override.aes=list(size=c(2.0, 2.5, 3.0), alpha=1)
    )
  )

# -----------------------------
# FIGURE4_PANELA_PANELB_HOTFIX_V1_2: Panel B facet/title overflow repair.
pB <- pB +
  scale_y_continuous(expand=expansion(mult=c(0.03, 0.30))) +
  theme(
    strip.text=element_text(size=10.1, face="bold", lineheight=0.92,
                            margin=margin(1.4, 1.6, 1.4, 1.6, "mm")),
    axis.text.x=element_text(size=8.9),
    axis.text.y=element_text(size=9.0),
    plot.margin=margin(2.5, 7.0, 2.0, 2.0, "mm")
  )

# C. Observed nominal counts versus permutation null
# -----------------------------
C[, feature_family := factor(as.character(feature_family), levels=family_levels_y)]
C[, enrichment_ratio := fifelse(
  is.finite(null_q95) & null_q95 > 0,
  observed / null_q95,
  NA_real_
)]
C[, observed_label := comma(round(observed, 0))]

null_points <- C[, .(
  feature_family,
  value=null_mean,
  component="Null mean"
)]
observed_points <- C[, .(
  feature_family,
  value=observed,
  component="Observed"
)]
point_data <- rbindlist(list(null_points, observed_points))
point_data[, component := factor(
  component,
  levels=c("Null mean", "Observed")
)]





pC <- ggplot(C, aes(y=feature_family)) +
  geom_segment(
    aes(
      x=null_mean,
      xend=null_q95,
      yend=feature_family
    ),
    colour=purple["soft"],
    linewidth=3.6,
    alpha=0.82,
    lineend="round"
  ) +
  geom_point(
    data=point_data,
    aes(x=value, fill=component, size=component),
    shape=21,
    colour=purple["ink"],
    stroke=0.34
  ) +
  geom_text(
    aes(x=observed, label=observed_label),
    nudge_x=max(c(C$observed, C$null_q95), na.rm=TRUE) * 0.022,
    hjust=0,
    vjust=0.5,
    size=3.05,
    lineheight=0.90,
    fontface="bold",
    colour=purple["ink"]
  ) +
  scale_fill_manual(
    values=c(
      "Null mean"=unname(purple["pale"]),
      "Observed"=unname(purple["deepest"])
    )
  ) +
  scale_size_manual(
    values=c("Null mean"=3.5, "Observed"=5.1)
  ) +
  scale_x_continuous(
    labels=comma,
    expand=expansion(mult=c(0, 0.18))
  ) +
  labs(
    x="Features with nominal p < 0.05\n(observed count labelled)",
    y="Feature family",
    tag="C"
  ) +
  theme_excellence(11.2) +
  theme(
    legend.position="none",
    axis.text.y=element_text(size=10.4)
  )

# -----------------------------
# D. Ranked exploratory features
# -----------------------------
D <- D[order(p_plot, neglog10p)]
if (nrow(D) > 14) D <- D[seq_len(14)]

if (!"feature_label_full" %in% names(D)) {
  D[, feature_label_full := feature_clean]
}
if (!"q_plot" %in% names(D)) D[, q_plot := NA_real_]

D[, feature_label_full := humanise_feature(feature_label_full)]
D[, display_label := wrap_text(feature_label_full, width=34)]
D[, display_label := make.unique(display_label)]
D[, display_label := factor(display_label, levels=rev(display_label))]
D[, method := factor(as.character(method), levels=method_order)]

family_present <- intersect(names(family_cols), unique(as.character(D$feature_family)))
method_present <- intersect(method_order, unique(as.character(D$method)))

pD <- ggplot(D, aes(x=neglog10p, y=display_label)) +
  geom_segment(
    aes(
      x=0,
      xend=neglog10p,
      yend=display_label
    ),
    colour=purple["light"],
    linewidth=1.15,
    alpha=0.82,
    lineend="round"
  ) +
  geom_point(
    aes(
      fill=method,
      shape=method
    ),
    size=4.5,
    colour=purple["ink"],
    stroke=0.42,
    alpha=0.98
  ) +
  geom_point(
    data=D[is.finite(q_plot) & q_plot < 0.05],
    aes(x=neglog10p, y=display_label),
    shape=21,
    size=6.0,
    fill=NA,
    colour=purple["deepest"],
    stroke=0.70,
    inherit.aes=FALSE
  ) +
  geom_text(
    data=D[is.finite(q_plot) & q_plot < 0.05],
    aes(
      x=neglog10p,
      y=display_label,
      label=paste0("q=", formatC(q_plot, format="fg", digits=2))
    ),
    inherit.aes=FALSE,
    hjust=-0.08,
    size=3.0,
    fontface="bold",
    colour=purple["deepest"]
  ) +
  scale_fill_manual(
    values=method_cols[method_present],
    breaks=method_present,
    drop=TRUE
  ) +
  scale_shape_manual(
    values=method_shapes[method_present],
    breaks=method_present,
    drop=TRUE
  ) +
  scale_x_continuous(
    expand=expansion(mult=c(0, 0.10))
  ) +
  labs(
    x=expression(-log[10]~"raw p"),
    y="Ranked exploratory feature",
    fill="Method",
    shape="Method",
    tag="D"
  ) +
  theme_excellence(10.2) +
  theme(
    axis.text.y=element_text(
      size=8.8,
      lineheight=0.94
    ),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.box.just="left"
  ) +
  guides(
    fill=guide_legend(
      title.position="top",
      nrow=1,
      override.aes=list(size=4.2)
    ),
    shape="none"
  )

# -----------------------------
# E. Candidate-support heatmap
# -----------------------------
if (!"feature_label_full" %in% names(E)) {
  E[, feature_label_full := feature_clean]
}
E[, method_display := fifelse(
  method == "CLR sensitivity",
  "CLR",
  as.character(method)
)]
method_levels <- c("ANCOM-BC2", "ALDEx2", "CLR", "Prevalence")
E[, method_display := factor(method_display, levels=method_levels)]

feature_order <- unique(
  D[, .(feature_clean, feature_label_full)]
)

E <- E[
  feature_clean %in% feature_order$feature_clean
]

E <- merge(
  E,
  feature_order,
  by="feature_clean",
  all.x=TRUE,
  suffixes=c("", "_ranked")
)

E[
  is.na(feature_label_full_ranked),
  feature_label_full_ranked := feature_label_full
]

E[
  ,
  feature_label_full_ranked :=
    humanise_feature(feature_label_full_ranked)
]

E[
  ,
  display_label :=
    wrap_text(feature_label_full_ranked, width=29)
]

rank_levels <- wrap_text(
  humanise_feature(feature_order$feature_label_full),
  width=29
)

E[, display_label := factor(
  display_label,
  levels=rev(rank_levels)
)]
E[, text_colour := fifelse(score >= 4.3, "white", purple["ink"])]

pE <- ggplot(E, aes(x=method_display, y=display_label, fill=score)) +
  geom_tile(
    colour="white",
    linewidth=0.60
  ) +
  geom_text(
    aes(label=tile_label, colour=text_colour),
    size=3.05,
    fontface="bold",
    show.legend=FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours=c(
      purple["faint"],
      purple["pale"],
      purple["soft"],
      purple["mid"],
      purple["deepest"]
    ),
    values=rescale(c(0, 1.5, 3.0, 5.0, 8.0)),
    limits=c(0, 8),
    oob=squish,
    name="Exploratory\nsupport score"
  ) +
  labs(
    x="Evidence layer",
    y="Exploratory feature",
    tag="E"
  ) +
  theme_excellence(9.9) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(angle=22, hjust=1, size=9.2),
    axis.text.y=element_text(size=8.3, lineheight=0.92),
    legend.position="right"
  ) +
  guides(
    fill=guide_colourbar(
      title.position="top",
      barheight=unit(38, "mm"),
      barwidth=unit(4.4, "mm")
    )
  )




# >>> PHASE_F4D_FIGURE4_LEGEND_SEMANTICS_CLEANUP_BEGIN
# Legend semantics cleanup only. Does not alter data or statistics.

figure4_label_text <- function(p) {
  if (!inherits(p, 'ggplot')) return('')
  vals <- lapply(
    p$labels,
    function(x) {
      tryCatch(
        paste(deparse(x), collapse = ' '),
        error = function(e) as.character(x)
      )
    }
  )
  paste(unlist(vals), collapse = ' | ')
}

figure4_find_plot_name <- function(patterns) {
  objects <- ls(envir = .GlobalEnv)
  for (nm in objects) {
    obj <- tryCatch(get(nm, envir = .GlobalEnv), error = function(e) NULL)
    if (!inherits(obj, 'ggplot')) next
    txt <- figure4_label_text(obj)
    hit <- any(vapply(patterns, function(pattern) grepl(pattern, txt, fixed = TRUE), logical(1)))
    if (hit) return(nm)
  }
  NULL
}

figure4_remove_scale <- function(p, aesthetic) {
  if (!inherits(p, 'ggplot')) return(p)
  keep <- vapply(p$scales$scales, function(scale) !(aesthetic %in% scale$aesthetics), logical(1))
  p$scales$scales <- p$scales$scales[keep]
  p
}

figure4_drop_geom_guides <- function(p) {
  p + ggplot2::guides(
    shape = 'none',
    colour = 'none',
    color = 'none',
    fill = 'none',
    linetype = 'none',
    alpha = 'none',
    size = 'none'
  )
}

figure4_add_manual_shape_legend <- function(p, title, labels, values, filled = NULL, colours = NULL) {
  helper <- data.frame(
    .legend_x = rep(Inf, length(labels)),
    .legend_y = rep(Inf, length(labels)),
    .legend_group = factor(labels, levels = labels)
  )

  p <- figure4_remove_scale(p, 'shape')
  p <- figure4_drop_geom_guides(p)

  if (is.null(filled)) filled <- rep(NA, length(labels))
  if (is.null(colours)) colours <- rep('black', length(labels))

  p +
    ggplot2::geom_point(
      data = helper,
      mapping = ggplot2::aes(
        x = .legend_x,
        y = .legend_y,
        shape = .legend_group
      ),
      inherit.aes = FALSE,
      alpha = 0,
      size = 3
    ) +
    ggplot2::scale_shape_manual(
      name = title,
      values = stats::setNames(values, labels),
      drop = FALSE
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        title = title,
        override.aes = list(
          alpha = 1,
          size = 3.4,
          fill = filled,
          colour = colours,
          stroke = 0.8
        ),
        order = 1
      )
    )
}

figure4_readable_legend <- function(p, position = 'bottom') {
  p + ggplot2::theme(
    legend.position = position,
    legend.direction = if (position == 'right') 'vertical' else 'horizontal',
    legend.box = if (position == 'right') 'vertical' else 'horizontal',
    legend.title = ggplot2::element_text(size = 10.8, face = 'bold'),
    legend.text = ggplot2::element_text(size = 9.8),
    legend.key.height = grid::unit(11, 'pt'),
    legend.key.width = grid::unit(11, 'pt'),
    plot.margin = ggplot2::margin(7, 9, 7, 9)
  )
}

figure4_remove_direct_text_layers <- function(p) {
  if (!inherits(p, 'ggplot')) return(p)
  keep <- vapply(
    p$layers,
    function(layer) {
      cls <- class(layer$geom)
      !any(cls %in% c('GeomText', 'GeomLabel'))
    },
    logical(1)
  )
  p$layers <- p$layers[keep]
  p
}

figure4_A_name <- figure4_find_plot_name(c('Tests with nominal p < 0.05'))
figure4_C_name <- figure4_find_plot_name(c('Features with nominal p < 0.05'))
figure4_D_name <- figure4_find_plot_name(c('Ranked exploratory feature'))

if (is.null(figure4_A_name)) stop('F4D could not resolve Panel A')
if (is.null(figure4_C_name)) stop('F4D could not resolve Panel C')
if (is.null(figure4_D_name)) stop('F4D could not resolve Panel D')

# Panel A: correct feature-class legend only.
figure4_A <- get(figure4_A_name, envir = .GlobalEnv)
figure4_A <- figure4_add_manual_shape_legend(
  figure4_A,
  title = 'Method',
  labels = c(),
  values = c(15, 17),
  filled = c('white', 'white'),
  colours = c('#6F6F6F', '#6F6F6F')
)
figure4_A <- figure4_readable_legend(figure4_A, 'bottom')
assign(figure4_A_name, figure4_A, envir = .GlobalEnv)

# Panel C: remove misleading association-evidence legend and show only component.
figure4_C <- get(figure4_C_name, envir = .GlobalEnv)
figure4_C <- figure4_add_manual_shape_legend(
  figure4_C,
  title = 'Component',
  labels = c('Null mean', 'Observed'),
  values = c(1, 16),
  filled = c(NA, '#3b2257'),
  colours = c('#7d74a8', '#3b2257')
)
figure4_C <- figure4_readable_legend(figure4_C, 'bottom')
assign(figure4_C_name, figure4_C, envir = .GlobalEnv)

# Panel D: remove blurred direct value callouts and show correct method shapes.
figure4_D <- get(figure4_D_name, envir = .GlobalEnv)
figure4_D <- figure4_remove_direct_text_layers(figure4_D)
figure4_D <- figure4_add_manual_shape_legend(
  figure4_D,
  title = 'Method',
  labels = c('ALDEx2', 'CLR sensitivity'),
  values = c(17, 15),
  filled = c('#3b2257', '#7d74a8'),
  colours = c('#3b2257', '#3b2257')
)
figure4_D <- figure4_readable_legend(figure4_D, 'bottom')
assign(figure4_D_name, figure4_D, envir = .GlobalEnv)

cat('F4D_RESOLVED_PANEL_OBJECTS=', figure4_A_name, ',', figure4_C_name, ',', figure4_D_name, '\n', sep = '')
# <<< PHASE_F4D_FIGURE4_LEGEND_SEMANTICS_CLEANUP_END


# >>> PHASE_F4E_FIGURE4_FINAL_MARKER_AND_Q_ANNOTATION_BEGIN
# Final two-item visual-QC repair only.

figure4_f4e_label_text <- function(p) {
  if (!inherits(p, "ggplot")) return("")
  vals <- lapply(
    p$labels,
    function(x) {
      tryCatch(
        paste(deparse(x), collapse = " "),
        error = function(e) as.character(x)
      )
    }
  )
  paste(unlist(vals), collapse = " | ")
}

figure4_f4e_find_plot_name <- function(patterns) {
  objects <- ls(envir = .GlobalEnv)
  for (nm in objects) {
    obj <- tryCatch(get(nm, envir = .GlobalEnv), error = function(e) NULL)
    if (!inherits(obj, "ggplot")) next
    txt <- figure4_f4e_label_text(obj)
    hit <- any(vapply(
      patterns,
      function(pattern) grepl(pattern, txt, fixed = TRUE),
      logical(1)
    ))
    if (hit) return(nm)
  }
  NULL
}

figure4_f4e_remove_scale <- function(p, aesthetic) {
  if (!inherits(p, "ggplot")) return(p)
  keep <- vapply(
    p$scales$scales,
    function(scale) !(aesthetic %in% scale$aesthetics),
    logical(1)
  )
  p$scales$scales <- p$scales$scales[keep]
  p
}

figure4_f4e_drop_helper_legend_layers <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  keep <- vapply(
    p$layers,
    function(layer) {
      df <- layer$data
      if (!is.data.frame(df)) return(TRUE)
      !all(c(".legend_x", ".legend_y", ".legend_group") %in% names(df))
    },
    logical(1)
  )
  p$layers <- p$layers[keep]
  p
}

figure4_f4e_find_method_column <- function(p) {
  data_candidates <- list(p$data)
  for (layer in p$layers) {
    if (is.data.frame(layer$data)) {
      data_candidates[[length(data_candidates) + 1L]] <- layer$data
    }
  }

  for (df in data_candidates) {
    if (!is.data.frame(df) || nrow(df) == 0L) next
    for (col in names(df)) {
      values <- unique(as.character(df[[col]]))
      values <- values[!is.na(values)]
      joined <- tolower(paste(values, collapse = " | "))
      if (
        grepl("aldex", joined, fixed = TRUE) &&
        grepl("clr", joined, fixed = TRUE)
      ) {
        return(col)
      }
    }
  }
  NULL
}

figure4_f4e_restore_method_markers <- function(p) {
  if (!inherits(p, "ggplot")) return(p)

  p <- figure4_f4e_drop_helper_legend_layers(p)
  method_col <- figure4_f4e_find_method_column(p)

  if (is.null(method_col)) {
    stop("F4E could not find ALDEx2/CLR method column in Panel A")
  }

  shape_mapping <- ggplot2::aes(shape = .data[[method_col]])$shape

  if (is.data.frame(p$data) && method_col %in% names(p$data)) {
    p$mapping$shape <- shape_mapping
  }

  mapped_point_layers <- 0L

  for (i in seq_along(p$layers)) {
    layer <- p$layers[[i]]
    geom_classes <- class(layer$geom)
    is_point <- any(geom_classes %in% c("GeomPoint"))
    if (!is_point) next

    layer_has_method <- (
      is.data.frame(layer$data) &&
      method_col %in% names(layer$data)
    )

    inherits_method <- (
      !is.data.frame(layer$data) &&
      is.data.frame(p$data) &&
      method_col %in% names(p$data)
    )

    if (layer_has_method || inherits_method) {
      layer$mapping$shape <- shape_mapping
      layer$aes_params$shape <- NULL
      p$layers[[i]] <- layer
      mapped_point_layers <- mapped_point_layers + 1L
    }
  }

  if (mapped_point_layers == 0L) {
    stop("F4E found method column but no Panel A point layer could be remapped")
  }

  values <- character(0)

  if (is.data.frame(p$data) && method_col %in% names(p$data)) {
    values <- unique(as.character(p$data[[method_col]]))
  }

  if (length(values) == 0L) {
    for (layer in p$layers) {
      if (
        is.data.frame(layer$data) &&
        method_col %in% names(layer$data)
      ) {
        values <- unique(as.character(layer$data[[method_col]]))
        if (length(values) > 0L) break
      }
    }
  }

  values <- values[!is.na(values)]
  shape_values <- rep(16, length(values))
  names(shape_values) <- values
  display_labels <- values

  for (i in seq_along(values)) {
    low <- tolower(values[[i]])
    if (grepl("aldex", low, fixed = TRUE)) {
      shape_values[[i]] <- 17
      display_labels[[i]] <- "ALDEx2"
    } else if (grepl("clr", low, fixed = TRUE)) {
      shape_values[[i]] <- 15
      display_labels[[i]] <- "CLR sensitivity"
    }
  }

  p <- figure4_f4e_remove_scale(p, "shape")

  p +
    ggplot2::scale_shape_manual(
      name = "Method",
      values = shape_values,
      labels = display_labels,
      drop = FALSE
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        title = "Method",
        override.aes = list(
          alpha = 1,
          size = 3.4
        ),
        order = 1
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
}

figure4_A_name <- figure4_f4e_find_plot_name(
  c("Tests with nominal p < 0.05")
)

figure4_D_name <- figure4_f4e_find_plot_name(
  c("Ranked exploratory feature")
)

if (is.null(figure4_A_name)) stop("F4E could not resolve Panel A")
if (is.null(figure4_D_name)) stop("F4E could not resolve Panel D")

# PANEL A ---------------------------------------------------------------
figure4_A <- get(figure4_A_name, envir = .GlobalEnv)
figure4_A <- figure4_f4e_restore_method_markers(figure4_A)
assign(figure4_A_name, figure4_A, envir = .GlobalEnv)

# PANEL D ---------------------------------------------------------------
figure4_D <- get(figure4_D_name, envir = .GlobalEnv)

figure4_D <- figure4_D +
  ggplot2::annotate(
    "text",
    x = Inf,
    y = "Lachnospiraceae unclassified",
    label = "q<0.05",
    hjust = 1.08,
    vjust = -0.55,
    size = 3.3,
    fontface = "bold",
    colour = "#3b2257"
  )

assign(figure4_D_name, figure4_D, envir = .GlobalEnv)

cat(
  "F4E_RESOLVED_PANEL_OBJECTS=",
  figure4_A_name, ",",
  figure4_D_name, "\n",
  sep = ""
)
# <<< PHASE_F4E_FIGURE4_FINAL_MARKER_AND_Q_ANNOTATION_END

top_row <- pA + pB + plot_layout(widths=c(0.88, 1.52))
bottom_row <- pD + pE + plot_layout(widths=c(1.02, 1.08))
figure4 <- top_row / pC / bottom_row +
  plot_layout(heights=c(1.03, 0.76, 1.25))

output_stem <- file.path(FIGDIR, "Figure_4_excellence_candidate")
pdf_out <- paste0(output_stem, ".pdf")
png_out <- paste0(output_stem, ".png")
tiff_out <- paste0(output_stem, ".tiff")
preview_out <- file.path(PREVIEWDIR, "Figure_4_excellence_preview.png")


# FIGURE4_READABILITY_EVIDENCE_PATCH_V1_1: display-only readability pass.
figure4_readability_theme_v1_1 <- theme(
  text = element_text(size = 12.4),
  axis.title = element_text(size = 14.0, face = "bold"),
  axis.text = element_text(size = 11.5),
  strip.text = element_text(size = 11.5, face = "bold", lineheight = 0.95),
  legend.title = element_text(size = 11.8, face = "bold"),
  legend.text = element_text(size = 10.8),
  plot.tag = element_text(size = 17.5, face = "bold")
)

for (figure4_nm in c("pA","pB","pC","pD","pE")) {
  if (exists(figure4_nm, inherits=FALSE)) {
    figure4_obj <- get(figure4_nm, inherits=FALSE)
    if (inherits(figure4_obj, "ggplot")) {
      assign(figure4_nm, figure4_obj + figure4_readability_theme_v1_1)
    }
  }
}

if (exists("fig", inherits=FALSE) && inherits(fig, "patchwork")) {
  fig <- fig & figure4_readability_theme_v1_1
}
if (exists("figure4", inherits=FALSE) && inherits(figure4, "patchwork")) {
  figure4 <- figure4 & figure4_readability_theme_v1_1
}

# FIGURE4_PANELA_PANELB_HOTFIX_V1_2: display-only readability pass.
figure4_readability_theme_v1_2 <- theme(
  text = element_text(size = 12.4),
  axis.title = element_text(size = 14.0, face = "bold"),
  axis.text = element_text(size = 11.5),
  strip.text = element_text(size = 11.2, face = "bold", lineheight = 0.95),
  legend.title = element_text(size = 11.8, face = "bold"),
  legend.text = element_text(size = 10.8),
  plot.tag = element_text(size = 17.5, face = "bold")
)

for (figure4_nm in c("pA","pB","pC","pD","pE")) {
  if (exists(figure4_nm, inherits=FALSE)) {
    figure4_obj <- get(figure4_nm, inherits=FALSE)
    if (inherits(figure4_obj, "ggplot")) {
      assign(figure4_nm, figure4_obj + figure4_readability_theme_v1_2)
    }
  }
}

if (exists("fig", inherits=FALSE) && inherits(fig, "patchwork")) {
  fig <- fig & figure4_readability_theme_v1_2
}
if (exists("figure4", inherits=FALSE) && inherits(figure4, "patchwork")) {
  figure4 <- figure4 & figure4_readability_theme_v1_2
}

ggsave(
  pdf_out,
  figure4,
  width=16.4,
  height=13.2,
  units="in",
  device=cairo_pdf,
  limitsize=FALSE
)
ggsave(
  png_out,
  figure4,
  width=16.4,
  height=13.2,
  units="in",
  dpi=400,
  limitsize=FALSE
)
ggsave(
  tiff_out,
  figure4,
  width=16.4,
  height=13.2,
  units="in",
  dpi=600,
  compression="lzw",
  limitsize=FALSE
)
ggsave(
  preview_out,
  figure4,
  width=16.4,
  height=13.2,
  units="in",
  dpi=240,
  limitsize=FALSE
)

audit <- data.table(
  item=c(
    "version",
    "source_A_rows",
    "source_B_rows",
    "source_C_rows",
    "source_D_rows",
    "source_E_rows",
    "display_D_rows",
    "output_pdf",
    "output_png",
    "output_tiff",
    "status"
  ),
  value=c(
    "v1_excellence_final_tidy",
    nrow(A),
    nrow(B),
    nrow(C),
    nrow(fread(input_files["D"], showProgress=FALSE)),
    nrow(E),
    nrow(D),
    pdf_out,
    png_out,
    tiff_out,
    "READY_FOR_FIGURE4_VISUAL_QC"
  )
)

fwrite(
  audit,
  file.path(REPORTDIR, "Figure_4_excellence_audit.tsv"),
  sep="\t"
)

manifest <- data.table(
  panel=names(input_files),
  source_file=unname(input_files),
  source_size_bytes=file.info(input_files)$size
)
fwrite(
  manifest,
  file.path(REPORTDIR, "Figure_4_excellence_source_manifest.tsv"),
  sep="\t"
)

cat("=== PFAS FIGURE 4 EXCELLENCE ===\n")
cat("Source rows A-E:", nrow(A), nrow(B), nrow(C), nrow(D), nrow(E), "\n")
cat("Output:", pdf_out, "\n")
cat("Status: READY_FOR_FIGURE4_VISUAL_QC\n")


# >>> PHASE_F4A_FIGURE4_EDITORIAL_REFINEMENT_BEGIN
figure4_apply_editorial_theme <- function(p) {
  p + theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.margin = margin(6, 8, 6, 8),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )
}

figure4_method_shape_scale <- scale_shape_manual(
  name = "Method",
  values = c("ALDEx2" = 17, "CLR sensitivity" = 15),
  drop = FALSE
)

figure4_count_shape_scale <- scale_shape_manual(
  name = "Association count",
  values = c("Nominal p < 0.05" = 21, "FDR < 0.05" = 24),
  drop = FALSE
)

figure4_summary_shape_scale <- scale_shape_manual(
  name = "Statistic",
  values = c("Nominal p < 0.05" = 21, "FDR < 0.05" = 24),
  drop = FALSE
)
# <<< PHASE_F4A_FIGURE4_EDITORIAL_REFINEMENT_END
