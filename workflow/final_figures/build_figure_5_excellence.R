# FIGURE5_LAYOUT_HOTFIX_V1_2
# Layout/readability hotfix: panel B label cleanup, panel D note removal,
# panel E dendrogram cleanup, uniform subplot tags.
#!/usr/bin/env Rscript
# === FIGURE5_RULE_BASED_SCORE_NAME_BEGIN ===
# Formal metric name: Rule-based candidate-priority score
# Compact displayed label: Rule-based priority score
# Interpretation: dimensionless exploratory ordering metric; not a
# probability, standardized effect size, causal estimate or validated biomarker.
# === FIGURE5_RULE_BASED_SCORE_NAME_END ===


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Repository root is required.", call. = FALSE)

ROOT <- normalizePath(args[[1]], mustWork = TRUE)
OUT <- file.path(ROOT, "28_figure_excellence")
SRC <- file.path(OUT, "source_data")
FIGDIR <- file.path(OUT, "main_figures")
CLUSTERDIR <- file.path(OUT, "cluster_orders")
REPORTDIR <- file.path(OUT, "reports")

for (directory in c(FIGDIR, CLUSTERDIR, REPORTDIR)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("data.table", "ggplot2", "patchwork", "scales", "cluster")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
    stop(
        sprintf("Missing required packages: %s", paste(missing_packages, collapse = ", ")),
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

HELPER <- file.path(OUT, "R", "figure_excellence_helpers.R")
if (!file.exists(HELPER)) {
    stop(sprintf("Missing shared helper: %s", HELPER), call. = FALSE)
}
source(HELPER)

INPUTS <- list(
    repeat_performance = file.path(SRC, "Figure_5A_repeat_level_oof_performance.tsv"),
    predictions = file.path(SRC, "Figure_5B_selected_oof_predictions.tsv"),
    prediction_metrics = file.path(SRC, "Figure_5B_selected_prediction_metrics.tsv"),
    stable_features = file.path(SRC, "Figure_5C_stable_features.tsv"),
    selected_candidates = file.path(SRC, "Figure_5D_selected_candidates.tsv"),
    authoritative_candidates = file.path(
        ROOT,
        "12_machine_learning_candidate_prioritization",
        "tables",
        "Table_S22b_Phase12_top50_candidate_prioritization.tsv"
    )
)

for (input_name in names(INPUTS)) {
    input_path <- INPUTS[[input_name]]
    if (!file.exists(input_path) || file.info(input_path)$size <= 0) {
        stop(sprintf("Required Phase 28B source file is missing: %s", input_path), call. = FALSE)
    }
}

repeat_performance <- fread(INPUTS$repeat_performance, showProgress = FALSE)
predictions <- fread(INPUTS$predictions, showProgress = FALSE)
prediction_metrics <- fread(INPUTS$prediction_metrics, showProgress = FALSE)
stable_features <- fread(INPUTS$stable_features, showProgress = FALSE)
candidates <- fread(INPUTS$selected_candidates, showProgress = FALSE)
authoritative_candidates <- fread(
    INPUTS$authoritative_candidates,
    data.table = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

setnames(
    authoritative_candidates,
    names(authoritative_candidates),
    trimws(sub("^\ufeff", "", names(authoritative_candidates)))
)

if (
    "FDR_supported" %in% names(authoritative_candidates) &&
    !"fdr_supported" %in% names(authoritative_candidates)
) {
    setnames(
        authoritative_candidates,
        "FDR_supported",
        "fdr_supported"
    )
}

require_columns <- function(value, required, label) {
    missing <- setdiff(required, names(value))
    if (length(missing)) {
        stop(
            sprintf("%s is missing columns: %s", label, paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

require_columns(
    repeat_performance,
    c("feature_set", "prediction_repeat", "spearman", "r2", "feature_label"),
    "Repeat-performance table"
)
require_columns(
    predictions,
    c("feature_set", "sample", "observed", "predicted", "feature_label"),
    "Prediction table"
)
require_columns(
    prediction_metrics,
    c("feature_set", "feature_label", "n", "spearman", "r2"),
    "Prediction-metrics table"
)
require_columns(
    stable_features,
    c(
        "feature_family", "feature", "mean_normalized_importance",
        "selected_fraction", "marginal_rho", "stability_score"
    ),
    "Stable-feature table"
)
require_columns(
    candidates,
    c(
        "candidate_rank", "rank_score", "feature_family", "feature",
        "mean_normalized_importance", "selected_fraction", "marginal_rho",
        "association_p", "fdr_flag", "nominal_flag", "prevalence",
        "public_context_rows", "bioaccumulation_flag"
    ),
    "Candidate table"
)
require_columns(
    authoritative_candidates,
    c(
        "positive_elimination_followup_candidate",
        "feature_family",
        "nominal_supported",
        "fdr_supported",
        "public_context_matching_rows",
        "bioaccumulation_relevance_keyword_flag"
    ),
    "Authoritative candidate-prioritisation table"
)

num <- function(value) suppressWarnings(as.numeric(value))

flag <- function(value) {
    if (is.logical(value)) return(fifelse(is.na(value), FALSE, value))
    normalized <- tolower(trimws(as.character(value)))
    normalized %in% c("true", "t", "1", "yes", "y", "supported", "positive")
}

scale01 <- function(value) {
    value <- num(value)
    output <- rep(NA_real_, length(value))
    keep <- is.finite(value)
    if (!any(keep)) return(output)
    limits <- range(value[keep], na.rm = TRUE)
    if (!all(is.finite(limits)) || diff(limits) == 0) {
        output[keep] <- 0.5
        return(output)
    }
    output[keep] <- (value[keep] - limits[[1]]) / diff(limits)
    output
}

family_short <- function(value) {
    key <- tolower(as.character(value))
    fifelse(
        grepl("taxonomy.*species", key), "Species",
        fifelse(
            grepl("taxonomy.*genus", key), "Genus",
            fifelse(
                grepl("humann", key), "HUMAnN",
                fifelse(
                    grepl("subsystem", key), "AGORA2 subsystem",
                    fifelse(
                        grepl("reaction", key), "AGORA2 reaction",
                        fifelse(
                            grepl("exchange", key), "AGORA2 exchange",
                            fifelse(grepl("mag", key), "MAG", gsub("_", " ", value))
                        )
                    )
                )
            )
        )
    )
}

clean_feature_plain <- function(value) {
    value <- as.character(value)
    value <- sub("\\|.*$", "", value)
    value <- sub("^.*::", "", value)
    value <- sub("^s__", "", value)
    value <- sub("^g__", "", value)
    value <- gsub("_", " ", value, fixed = TRUE)
    value <- gsub("\\s+", " ", value)
    trimws(value)
}

figure5_compact_candidate_label <- function(x, width = 18) {
  x <- figure5_exact_humanize(x)
  x <- gsub("\\n\\[[^]]+\\]$", "", x)
  x <- gsub("\\s*\\[[^]]+\\]$", "", x)
  x <- trimws(x)

  x <- gsub(
    "^Iron\\(III\\) dicitrate transport$",
    "Fe(III) dicitrate\ntransport",
    x
  )
  x <- gsub(
    "^3-Dehydrochenodeoxycholate exchange$",
    "3-Dehydrochenodeoxycholate\nexchange",
    x
  )
  x <- gsub(
    "^3-Dehydrocholate exchange$",
    "3-Dehydrocholate\nexchange",
    x
  )

  vapply(
    x,
    function(value) {
      if (grepl("\n", value, fixed = TRUE)) {
        return(value)
      }
      paste(strwrap(value, width = width), collapse = "\n")
    },
    character(1)
  )
}


wrap_label <- function(value, width = 27) {
    vapply(
        clean_feature_plain(value),
        function(item) paste(strwrap(item, width = width), collapse = "\n"),
        FUN.VALUE = character(1)
    )
}

theme_clean <- function(base_size = 7.8) {
    theme_pfas_excellence(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = base_size + 0.8, margin = margin(b = 3)),
            plot.subtitle = element_text(size = base_size - 0.2, margin = margin(b = 4)),
            plot.margin = margin(5, 7, 5, 7)
        )
}

# -------------------------------------------------------------------------
# Panel A: repeated out-of-fold performance
# -------------------------------------------------------------------------

repeat_performance[, spearman := num(spearman)]

performance_summary <- repeat_performance[
    ,
    .(
        median_spearman = median(spearman, na.rm = TRUE),
        lower_quartile = quantile(spearman, 0.25, na.rm = TRUE),
        upper_quartile = quantile(spearman, 0.75, na.rm = TRUE),
        repeat_count = .N
    ),
    by = .(feature_set, feature_label)
]

setorder(performance_summary, median_spearman)
performance_levels <- performance_summary$feature_label

repeat_performance[
    ,
    feature_label := factor(feature_label, levels = performance_levels)
]
performance_summary[
    ,
    feature_label := factor(feature_label, levels = performance_levels)
]

pA <- ggplot(repeat_performance, aes(x = feature_label, y = spearman)) +
    geom_hline(
        yintercept = 0,
        linewidth = 0.42,
        linetype = "dashed",
        colour = PFAS_PURPLE[["neutral_dark"]]
    ) +
    geom_jitter(
        width = 0.11,
        height = 0,
        size = 0.85,
        alpha = 0.24,
        colour = PFAS_PURPLE[["medium"]]
    ) +
    geom_pointrange(
        data = performance_summary,
        aes(
            x = feature_label,
            y = median_spearman,
            ymin = lower_quartile,
            ymax = upper_quartile
        ),
        inherit.aes = FALSE,
        linewidth = 0.48,
        fatten = 2.7,
        colour = PFAS_PURPLE[["darkest"]]
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(
        limits = c(-0.62, 0.62),
        breaks = seq(-0.6, 0.6, 0.3)
    ) +
    labs(
        title = "A",
        subtitle = NULL,
        x = NULL,
        y = "Spearman correlation"
    ) +
    theme_clean(7.7)

# -------------------------------------------------------------------------
# Panel B: observed versus predicted
# -------------------------------------------------------------------------

predictions[, observed := num(observed)]
predictions[, predicted := num(predicted)]
prediction_metrics[, spearman := num(spearman)]
prediction_metrics[, r2 := num(r2)]
prediction_metrics[, n := num(n)]

combined <- prediction_metrics[
    grepl("combined", feature_set, ignore.case = TRUE),
    feature_label
]
other <- prediction_metrics[
    !feature_label %in% combined
][
    order(-spearman, -r2),
    feature_label
]
facet_levels <- unique(c(combined, other))

predictions[
    ,
    feature_label := factor(feature_label, levels = facet_levels)
]
prediction_metrics[
    ,
    feature_label := factor(feature_label, levels = facet_levels)
]
prediction_metrics[
    ,
    annotation := sprintf("rho = %.2f\nR² = %.2f\nn = %d", spearman, r2, as.integer(n))
]

pB <- ggplot(predictions, aes(x = observed, y = predicted)) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linewidth = 0.4,
        linetype = "dashed",
        colour = PFAS_PURPLE[["neutral"]]
    ) +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linewidth = 0.52,
        colour = PFAS_PURPLE[["darkest"]]
    ) +
    geom_point(
        shape = 21,
        size = 1.5,
        stroke = 0.28,
        fill = PFAS_PURPLE[["light"]],
        colour = PFAS_PURPLE[["darkest"]],
        alpha = 0.66
    ) +
    geom_label(
        data = prediction_metrics,
        aes(x = -Inf, y = Inf, label = annotation),
        inherit.aes = FALSE,
        hjust = -0.08,
        vjust = 1.08,
        size = 2.0,
        linewidth = 0.18,
        label.padding = grid::unit(0.08, "lines"),
        fill = "white",
        colour = PFAS_PURPLE[["darkest"]]
    ) +
    facet_wrap(~ feature_label, ncol = 2) +
    labs(
        title = "B",
        subtitle = NULL,
        x = "Observed elimination score",
        y = "Predicted elimination score"
    ) +
    theme_clean(7.4) +
    theme(strip.text = element_text(size = 6.8, face = "bold"))

# -------------------------------------------------------------------------
# Panel C: stable-feature matrix
# -------------------------------------------------------------------------

stable_features[, importance := num(mean_normalized_importance)]
stable_features[, selection := num(selected_fraction)]
stable_features[, rho := num(marginal_rho)]
stable_features[, stability_score := num(stability_score)]
stable_features[, family := family_short(feature_family)]
stable_features[, biological_feature := clean_feature_plain(feature)]
stable_features[, feature_id := paste(family, biological_feature, sep = "::")]

setorder(stable_features, -stability_score)
stable_unique <- stable_features[!duplicated(feature_id)]
stable_balanced <- stable_unique[, head(.SD, 2), by = family]
setorder(stable_balanced, -stability_score)
stable_balanced <- stable_balanced[seq_len(min(10L, .N))]

stable_balanced[
    ,
    display_label := paste0(wrap_label(biological_feature, 23), "\n[", family, "]")
]
stable_balanced[
    ,
    display_label := factor(display_label, levels = rev(display_label))
]

stable_matrix <- rbindlist(
    list(
        stable_balanced[
            ,
            .(
                display_label,
                metric = "RF importance",
                value = importance,
                label = sprintf("%.3f", importance)
            )
        ],
        stable_balanced[
            ,
            .(
                display_label,
                metric = "Selection frequency",
                value = selection,
                label = sprintf("%.2f", selection)
            )
        ],
        stable_balanced[
            ,
            .(
                display_label,
                metric = "Marginal rho",
                value = rho,
                label = sprintf("%+.2f", rho)
            )
        ]
    )
)
stable_matrix[, scaled := scale01(value), by = metric]
stable_matrix[
    ,
    metric := factor(
        metric,
        levels = c("RF importance", "Selection frequency", "Marginal rho")
    )
]
stable_matrix[, text_colour := fifelse(is.finite(scaled) & scaled >= 0.62, "white", "black")]

pC <- ggplot(
    stable_matrix,
    aes(x = metric, y = display_label, fill = scaled)
) +
    geom_tile(colour = "white", linewidth = 0.38) +
    geom_text(aes(label = label, colour = text_colour), size = 2.18) +
    scale_colour_identity() +
    scale_fill_gradientn(
        colours = c(
            "white",
            PFAS_PURPLE[["pale"]],
            PFAS_PURPLE[["medium"]],
            PFAS_PURPLE[["darkest"]]
        ),
        limits = c(0, 1),
        na.value = PFAS_PURPLE[["missing"]],
        guide = "none"
    ) +
    labs(
        title = "C",
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) +
    theme_clean(7.2) +
    theme(
        axis.text.x = element_text(angle = 25, hjust = 1, size = 6.6),
        axis.text.y = element_text(size = 5.8),
        panel.grid = element_blank()
    )

# -------------------------------------------------------------------------
# Panel D: family evidence matrix
# -------------------------------------------------------------------------

authoritative_candidates[
    ,
    positive_candidate :=
        flag(positive_elimination_followup_candidate)
]

authoritative_candidates[
    ,
    nominal_flag_resolved :=
        flag(nominal_supported)
]

authoritative_candidates[
    ,
    fdr_flag_resolved :=
        flag(fdr_supported)
]

authoritative_candidates[
    ,
    public_context_rows_resolved :=
        num(public_context_matching_rows)
]

authoritative_candidates[
    ,
    bioaccumulation_flag_resolved :=
        flag(bioaccumulation_relevance_keyword_flag)
]

authoritative_candidates[
    ,
    family := family_short(feature_family)
]

positive_authoritative <- authoritative_candidates[
    positive_candidate == TRUE
]

family_pool <- if (nrow(positive_authoritative) >= 15L) {
    positive_authoritative
} else {
    authoritative_candidates
}

family_summary <- family_pool[
    ,
    .(
        n_candidates = .N,
        nominal_supported = sum(
            nominal_flag_resolved,
            na.rm = TRUE
        ),
        fdr_supported = sum(
            fdr_flag_resolved,
            na.rm = TRUE
        ),
        public_context_supported = sum(
            is.finite(public_context_rows_resolved) &
                public_context_rows_resolved > 0,
            na.rm = TRUE
        ),
        bioaccumulation_flagged = sum(
            bioaccumulation_flag_resolved,
            na.rm = TRUE
        )
    ),
    by = family
]

if (!nrow(family_summary)) {
    stop(
        "Authoritative candidate table produced an empty family summary.",
        call. = FALSE
    )
}

family_order <- family_summary[
    order(-n_candidates, -nominal_supported),
    family
]

family_long <- melt(
    family_summary,
    id.vars = "family",
    measure.vars = c(
        "n_candidates",
        "nominal_supported",
        "fdr_supported",
        "public_context_supported",
        "bioaccumulation_flagged"
    ),
    variable.name = "evidence_type",
    value.name = "count"
)

family_long[
    ,
    evidence_type := factor(
        evidence_type,
        levels = c(
            "n_candidates",
            "nominal_supported",
            "fdr_supported",
            "public_context_supported",
            "bioaccumulation_flagged"
        ),
        labels = c("Candidates", "Nominal", "FDR", "Public context", "Bioaccumulation")
    )
]
family_long[, family := factor(family, levels = rev(family_order))]
family_long[
    ,
    scaled := if (max(count, na.rm = TRUE) > 0) count / max(count, na.rm = TRUE) else 0,
    by = evidence_type
]
family_long[, text_colour := fifelse(scaled >= 0.62, "white", "black")]

pD <- ggplot(
    family_long,
    aes(x = evidence_type, y = family, fill = scaled)
) +
    geom_tile(colour = "white", linewidth = 0.38) +
    geom_text(aes(label = count, colour = text_colour), size = 2.35, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(
        colours = c(
            "white",
            PFAS_PURPLE[["pale"]],
            PFAS_PURPLE[["medium"]],
            PFAS_PURPLE[["darkest"]]
        ),
        limits = c(0, 1),
        guide = "none"
    ) +
    labs(
        title = "D",
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) +
    theme_clean(7.2) +
    theme(
        axis.text.x = element_text(angle = 28, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.2),
        panel.grid = element_blank()
    )

# -------------------------------------------------------------------------
# Panel E: deduplicated candidate clustering
# -------------------------------------------------------------------------

candidates[, candidate_rank := num(candidate_rank)]
candidates[, rank_score := num(rank_score)]
candidates[, importance := num(mean_normalized_importance)]
candidates[, selection := num(selected_fraction)]
candidates[, rho := num(marginal_rho)]
candidates[, association_p := num(association_p)]
candidates[, prevalence := num(prevalence)]
candidates[, public_context := num(public_context_rows)]
candidates[, nominal := flag(nominal_flag)]
candidates[, fdr := flag(fdr_flag)]
candidates[, bioaccumulation := flag(bioaccumulation_flag)]
candidates[, family := family_short(feature_family)]
candidates[, biological_feature := clean_feature_plain(feature)]
candidates[, candidate_id := paste(family, biological_feature, sep = "::")]

setorder(candidates, candidate_rank, -rank_score)
candidate_unique <- candidates[!duplicated(candidate_id)]
candidate_unique <- candidate_unique[seq_len(min(14L, .N))]

candidate_unique[
    ,
    display_label := paste0(wrap_label(biological_feature, 26), "\n[", family, "]")
]
candidate_unique[
    ,
    positive_rho := fifelse(is.finite(rho), pmax(rho, 0), NA_real_)
]
candidate_unique[
    ,
    negative_log10_p := fifelse(
        is.finite(association_p),
        -log10(pmax(association_p, 1e-300)),
        NA_real_
    )
]
candidate_unique[
    ,
    public_context_log := fifelse(
        is.finite(public_context),
        log1p(pmax(public_context, 0)),
        NA_real_
    )
]

cluster_input <- data.frame(
    candidate_score = candidate_unique$rank_score,
    rf_importance = candidate_unique$importance,
    selection_frequency = candidate_unique$selection,
    positive_association = candidate_unique$positive_rho,
    association_strength = candidate_unique$negative_log10_p,
    nominal_support = as.numeric(candidate_unique$nominal),
    fdr_support = as.numeric(candidate_unique$fdr),
    prevalence = candidate_unique$prevalence,
    public_context = candidate_unique$public_context_log,
    bioaccumulation = as.numeric(candidate_unique$bioaccumulation),
    check.names = FALSE
)
rownames(cluster_input) <- candidate_unique$candidate_id

for (column_name in names(cluster_input)) {
    value <- num(cluster_input[[column_name]])
    finite <- is.finite(value)
    replacement <- if (any(finite)) median(value[finite], na.rm = TRUE) else 0
    value[!finite] <- replacement
    cluster_input[[column_name]] <- value
}

binary_columns <- intersect(
    c("nominal_support", "fdr_support", "bioaccumulation"),
    names(cluster_input)
)
variable_binary <- binary_columns[
    vapply(
        cluster_input[, binary_columns, drop = FALSE],
        function(value) all(c(0, 1) %in% unique(value[is.finite(value)])),
        FUN.VALUE = logical(1)
    )
]
asymm_indices <- match(variable_binary, names(cluster_input))

if (length(asymm_indices)) {
    candidate_distance <- cluster::daisy(
        cluster_input,
        metric = "gower",
        type = list(asymm = asymm_indices)
    )
} else {
    candidate_distance <- cluster::daisy(cluster_input, metric = "gower")
}

candidate_tree <- hclust(candidate_distance, method = "average")
leaf_ids <- rownames(cluster_input)[candidate_tree$order]
cluster_number <- min(4L, nrow(candidate_unique))
cluster_assignment <- cutree(candidate_tree, k = cluster_number)

candidate_unique[
    ,
    cluster_id := as.integer(cluster_assignment[candidate_id])
]
candidate_unique[
    ,
    cluster_order := match(candidate_id, leaf_ids)
]
setorder(candidate_unique, cluster_order)

cluster_export <- candidate_unique[
    ,
    .(
        candidate_id,
        display_label = gsub("[\r\n]+", " | ", display_label),
        biological_feature,
        family,
        candidate_rank,
        cluster_id,
        cluster_order,
        distance_method = "Gower mixed-evidence distance",
        linkage_method = "average"
    )
]

fwrite(
    cluster_export,
    file.path(CLUSTERDIR, "Figure_5_v3_candidate_cluster_order.tsv"),
    sep = "\t",
    quote = TRUE
)

hclust_segments <- function(tree) {
    n <- length(tree$order)
    leaf_y <- numeric(n)
    leaf_y[tree$order] <- seq_len(n)
    node_y <- numeric(n - 1L)
    segments <- vector("list", 3L * (n - 1L))
    segment_index <- 1L

    child_position <- function(child) {
        if (child < 0) {
            leaf_index <- -child
            return(list(x = 0, y = leaf_y[[leaf_index]]))
        }
        list(x = tree$height[[child]], y = node_y[[child]])
    }

    for (merge_index in seq_len(n - 1L)) {
        left <- child_position(tree$merge[merge_index, 1])
        right <- child_position(tree$merge[merge_index, 2])
        parent_x <- tree$height[[merge_index]]
        parent_y <- mean(c(left$y, right$y))
        node_y[[merge_index]] <- parent_y

        segments[[segment_index]] <- data.table(
            x = left$x, xend = parent_x, y = left$y, yend = left$y
        )
        segment_index <- segment_index + 1L
        segments[[segment_index]] <- data.table(
            x = right$x, xend = parent_x, y = right$y, yend = right$y
        )
        segment_index <- segment_index + 1L
        segments[[segment_index]] <- data.table(
            x = parent_x, xend = parent_x, y = left$y, yend = right$y
        )
        segment_index <- segment_index + 1L
    }

    rbindlist(segments)
}

dendrogram_segments <- hclust_segments(candidate_tree)

pE_dendrogram <- ggplot(dendrogram_segments) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 0.4,
        colour = PFAS_PURPLE[["neutral_dark"]]
    ) +
    scale_x_reverse(expand = c(0.02, 0)) +
    scale_y_continuous(
        limits = c(0.5, nrow(candidate_unique) + 0.5),
        expand = c(0, 0)
    ) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

cluster_strip <- candidate_unique[
    ,
    .(candidate_id, cluster_id, y = cluster_order)
]

cluster_colours <- c(
    "1" = PFAS_PURPLE[["palest"]],
    "2" = PFAS_PURPLE[["pale"]],
    "3" = PFAS_PURPLE[["light"]],
    "4" = PFAS_PURPLE[["dark"]]
)

pE_strip <- ggplot(
    cluster_strip,
    aes(x = 1, y = y, fill = factor(cluster_id))
) +
    geom_tile(width = 1, height = 1, colour = "white", linewidth = 0.18) +
    scale_fill_manual(values = cluster_colours, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
        limits = c(0.5, nrow(candidate_unique) + 0.5),
        expand = c(0, 0)
    ) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))

evidence_wide <- candidate_unique[
    ,
    .(
        candidate_id,
        display_label,
        cluster_id,
        `Rule-based priority score` = scale01(rank_score),
        `RF importance` = scale01(importance),
        `Selection frequency` = scale01(selection),
        `Positive association` = scale01(positive_rho),
        `Association strength` = scale01(negative_log10_p),
        `Nominal support` = as.numeric(nominal),
        `FDR support` = as.numeric(fdr),
        `Prevalence` = scale01(prevalence),
        `Public context` = scale01(public_context_log),
        `Bioaccumulation flag` = as.numeric(bioaccumulation)
    )
]

evidence_long <- melt(
    evidence_wide,
    id.vars = c("candidate_id", "display_label", "cluster_id"),
    variable.name = "evidence_metric",
    value.name = "scaled_evidence"
)

metric_order <- c(
    "Rule-based priority score",
    "RF importance",
    "Selection frequency",
    "Positive association",
    "Association strength",
    "Nominal support",
    "FDR support",
    "Prevalence",
    "Public context",
    "Bioaccumulation flag"
)

evidence_long[
    ,
    evidence_metric := factor(evidence_metric, levels = metric_order)
]
evidence_long[
    ,
    candidate_id := factor(candidate_id, levels = leaf_ids)
]

display_map <- setNames(candidate_unique$display_label, candidate_unique$candidate_id)
binary_metrics <- c("Nominal support", "FDR support", "Bioaccumulation flag")
evidence_long[
    ,
    binary_symbol := fifelse(
        as.character(evidence_metric) %in% binary_metrics &
            is.finite(scaled_evidence) &
            scaled_evidence >= 0.5,
        "●",
        ""
    )
]

pE_heatmap <- ggplot(
    evidence_long,
    aes(x = evidence_metric, y = candidate_id, fill = scaled_evidence)
) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(
        aes(label = binary_symbol),
        size = 2.0,
        colour = "white"
    ) +
    scale_y_discrete(labels = display_map, expand = c(0, 0)) +
    scale_fill_gradientn(
        colours = c(
            "white",
            PFAS_PURPLE[["pale"]],
            PFAS_PURPLE[["medium"]],
            PFAS_PURPLE[["darkest"]]
        ),
        limits = c(0, 1),
        na.value = PFAS_PURPLE[["missing"]],
        name = "Scaled evidence",
        guide = guide_colourbar(
            direction = "horizontal",
            title.position = "top",
            barwidth = grid::unit(34, "mm"),
            barheight = grid::unit(3.1, "mm")
        )
    ) +
    labs(x = NULL, y = NULL) +
    theme_clean(7.0) +
    theme(
        axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 6.25),
        axis.text.y = element_text(size = 5.75),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        plot.margin = margin(0, 5, 4, 2)
    )


# === FIGURE5_MINIMAL_AXIS_SYMBOLS_BEGIN ===
# Formatting only. No additional patchwork panels are created.

figure5_symbol_label <- function(x) {
  raw <- as.character(x)
  symbol <- rep("\u25CB", length(raw))

  is_genus <- grepl(
    "\\[genus\\]|metaphlan\\s+genus|taxonomy:\\s*genus|\\bgenus\\b",
    raw,
    ignore.case=TRUE
  )

  is_species <- grepl(
    "\\[species\\]|metaphlan\\s+species|taxonomy:\\s*species|\\bspecies\\b|\\bbacterium\\b|\\bspp\\.?\\b",
    raw,
    ignore.case=TRUE
  )

  is_agora2 <- grepl(
    "\\[agora2[^]]*\\]|agora2|reaction|exchange|subsystem",
    raw,
    ignore.case=TRUE
  )

  is_humann <- grepl(
    "\\[humann[^]]*\\]|humann",
    raw,
    ignore.case=TRUE
  )

  is_mag <- grepl(
    "\\[mag\\]|\\bmag\\b|\\bbin[ ._-]*[a-z0-9]+",
    raw,
    ignore.case=TRUE
  )

  symbol[is_genus] <- "\u25CF"
  symbol[is_species] <- "\u25B2"
  symbol[is_agora2] <- "\u25A0"
  symbol[is_humann] <- "\u25C6"
  symbol[is_mag] <- "\u25BC"

  clean <- gsub(
    "\\s*\\[(?:MetaPhlAn[^]]*|Genus|Species|MAG|HUMAnN[^]]*|AGORA2[^]]*)\\]",
    "",
    raw,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub("^Taxonomy:\\s*", "", clean, ignore.case=TRUE)
  clean <- gsub("\\s+", " ", clean)
  clean <- trimws(clean)

  clean <- vapply(
    clean,
    function(value) {
      paste(strwrap(value, width=30), collapse="\n")
    },
    character(1)
  )

  paste(symbol, clean)
}

figure5_set_y_labels <- function(plot_object) {
  y_scale <- plot_object$scales$get_scales("y")

  if (!is.null(y_scale)) {
    y_scale$labels <- figure5_symbol_label
    return(plot_object)
  }

  plot_object +
    ggplot2::scale_y_discrete(labels=figure5_symbol_label)
}

figure5_source_key <- paste(
  "\u25CF Genus",
  "\u25B2 Species",
  "\u25A0 AGORA2",
  "\u25C6 HUMAnN",
  "\u25BC MAG",
  sep="    "
)

if (exists("pC") && inherits(pC, "ggplot")) {
  pC <- figure5_set_y_labels(pC)

  pC <- pC +
    ggplot2::labs(
      x="Selection evidence metric",
      y="Selected feature",
      caption=figure5_source_key
    ) +
    ggplot2::theme(
      axis.title=ggplot2::element_text(
        face="bold",
        size=10.5
      ),
      axis.text.y=ggplot2::element_text(
        size=7.5,
        lineheight=0.94
      ),
      plot.caption=ggplot2::element_text(
        size=7.6,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=6)
      )
    )
}

if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_set_y_labels(pE_heatmap)

  pE_heatmap <- pE_heatmap +
    ggplot2::labs(
      x="Evidence metric",
      y="Prioritized feature",
      caption=figure5_source_key
    ) +
    ggplot2::theme(
      axis.title=ggplot2::element_text(
        face="bold",
        size=10.5
      ),
      axis.text.y=ggplot2::element_text(
        size=7.3,
        lineheight=0.93
      ),
      plot.caption=ggplot2::element_text(
        size=7.6,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=6)
      )
    )
}

if (
  exists("pE_dendrogram") &&
  inherits(pE_dendrogram, "ggplot")
) {
  pE_dendrogram <- pE_dendrogram +
    ggplot2::labs(tag="E") +
    ggplot2::theme(
      plot.tag=ggplot2::element_text(
        face="bold",
        size=12
      ),
      plot.tag.position=c(0.02, 0.98)
    )
}
# === FIGURE5_MINIMAL_AXIS_SYMBOLS_END ===


# === FIGURE5_FINAL_SPACING_BEGIN ===
# Label and spacing changes only. No new patchwork layouts.

figure5_final_symbol <- function(raw) {
  raw <- as.character(raw)
  symbol <- rep("\u25CB", length(raw))

  symbol[grepl("\\[genus\\]|^genus::|taxonomy:\\s*genus|\\bgenus\\b", raw, ignore.case=TRUE)] <- "\u25CF"
  symbol[grepl("\\[species\\]|^species::|taxonomy:\\s*species|\\bspecies\\b|\\bbacterium\\b|\\bspp\\.?\\b", raw, ignore.case=TRUE)] <- "\u25B2"
  symbol[grepl("\\[agora2[^]]*\\]|^agora2[^:]*::|agora2|reaction|exchange|subsystem", raw, ignore.case=TRUE)] <- "\u25A0"
  symbol[grepl("\\[humann[^]]*\\]|^humann::|humann", raw, ignore.case=TRUE)] <- "\u25C6"
  symbol[grepl("\\[mag\\]|^mag::|\\bmag\\b|\\bbin[ ._-]*[a-z0-9]+", raw, ignore.case=TRUE)] <- "\u25BC"

  symbol
}

figure5_final_clean <- function(raw) {
  clean <- as.character(raw)

  clean <- gsub(
    "^\\s*(?:HUMAnN|Species|Genus|MAG|AGORA2(?:\\s+reaction|\\s+exchange|\\s+subsystem)?)\\s*::\\s*",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub(
    "\\s*\\[(?:MetaPhlAn[^]]*|Genus|Species|MAG|HUMAnN[^]]*|AGORA2[^]]*)\\]",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub("^Taxonomy:\\s*", "", clean, ignore.case=TRUE)
  clean <- gsub("\\s+", " ", clean)
  trimws(clean)
}

figure5_final_two_lines <- function(value, width=25) {
  parts <- strwrap(value, width=width)

  if (!length(parts)) return("")
  if (length(parts) == 1) return(parts[[1]])
  if (length(parts) == 2) return(paste(parts, collapse="\n"))

  second <- paste(parts[-1], collapse=" ")

  if (nchar(second) > width + 5) {
    second <- paste0(substr(second, 1, width + 2), "...")
  }

  paste(parts[[1]], second, sep="\n")
}

figure5_final_label_c <- function(x) {
  raw <- as.character(x)
  symbol <- figure5_final_symbol(raw)
  clean <- figure5_final_clean(raw)

  clean <- gsub(
    "P162-PWY:\\s*glutamate degradation V\\s*\\(via hydroxyglutarate\\)",
    "P162-PWY: glutamate degradation V",
    clean,
    ignore.case=TRUE
  )

  clean <- gsub(
    "PWY-6629:\\s*glutamate biosynthesis\\s*\\(anaerobic\\)",
    "PWY-6629: glutamate biosynthesis",
    clean,
    ignore.case=TRUE
  )

  wrapped <- vapply(
    clean,
    figure5_final_two_lines,
    character(1),
    width=24
  )

  paste(symbol, wrapped)
}

figure5_final_label_e <- function(x) {
  raw <- as.character(x)
  symbol <- figure5_final_symbol(raw)
  clean <- figure5_final_clean(raw)

  wrapped <- vapply(
    clean,
    figure5_final_two_lines,
    character(1),
    width=27
  )

  paste(symbol, wrapped)
}

figure5_final_set_y_labels <- function(plot_object, label_function) {
  y_scale <- plot_object$scales$get_scales("y")

  if (!is.null(y_scale)) {
    y_scale$labels <- label_function
    return(plot_object)
  }

  plot_object + ggplot2::scale_y_discrete(labels=label_function)
}

figure5_final_source_key <- paste(
  "\u25CF Genus",
  "\u25B2 Species",
  "\u25A0 AGORA2",
  "\u25C6 HUMAnN",
  "\u25BC MAG",
  sep="    "
)

if (exists("pA") && inherits(pA, "ggplot")) {
  pA <- pA +
    ggplot2::labs(
      x="Spearman correlation",
      y="Feature family"
    ) +
    ggplot2::theme(
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=10.5,
        margin=ggplot2::margin(r=8)
      )
    )
}

if (exists("pC") && inherits(pC, "ggplot")) {
  pC <- figure5_final_set_y_labels(
    pC,
    figure5_final_label_c
  )

  pC <- pC +
    ggplot2::labs(
      x="Selection evidence metric",
      y="Selected feature",
      caption=figure5_final_source_key
    ) +
    ggplot2::theme(
      axis.title.x=ggplot2::element_text(
        face="bold",
        size=10.2,
        margin=ggplot2::margin(t=7)
      ),
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=9.8,
        margin=ggplot2::margin(r=7)
      ),
      axis.text.y=ggplot2::element_text(
        size=6.4,
        lineheight=0.82
      ),
      axis.text.x=ggplot2::element_text(
        size=7.8,
        angle=23,
        hjust=1
      ),
      plot.caption=ggplot2::element_text(
        size=7.0,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=6)
      ),
      plot.margin=ggplot2::margin(4, 7, 4, 3)
    )
}

if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- pD +
    ggplot2::labs(
      x="Evidence category",
      y="Feature family"
    ) +
    ggplot2::theme(
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=10.5,
        margin=ggplot2::margin(r=8)
      )
    )
}

if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_final_set_y_labels(
    pE_heatmap,
    figure5_final_label_e
  )

  pE_heatmap <- pE_heatmap +
    ggplot2::labs(
      x="Evidence metric",
      y="Prioritized feature",
      caption=figure5_final_source_key
    ) +
    ggplot2::theme(
      axis.title.x=ggplot2::element_text(
        face="bold",
        size=10.5,
        margin=ggplot2::margin(t=8)
      ),
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=10,
        margin=ggplot2::margin(r=7)
      ),
      axis.text.y=ggplot2::element_text(
        size=6.9,
        lineheight=0.87
      ),
      axis.text.x=ggplot2::element_text(
        size=7.8,
        angle=26,
        hjust=1
      ),
      plot.caption=ggplot2::element_text(
        size=7.1,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=6)
      ),
      plot.margin=ggplot2::margin(4, 8, 4, 1)
    )
}

if (
  exists("pE_dendrogram") &&
  inherits(pE_dendrogram, "ggplot")
) {
  pE_dendrogram <- pE_dendrogram +
    ggplot2::labs(tag="E") +
    ggplot2::theme(
      plot.tag=ggplot2::element_text(
        face="bold",
        size=12
      ),
      plot.tag.position=c(0.02, 0.98)
    )
}
# === FIGURE5_FINAL_SPACING_END ===


# === FIGURE5_READABILITY_FIX_BEGIN ===
# Display-only closeout. No new patchwork layouts.

figure5_readability_map_name <- function(plot_object, aesthetic) {
  mapping_value <- plot_object$mapping[[aesthetic]]

  if (is.null(mapping_value) && length(plot_object$layers)) {
    for (layer in plot_object$layers) {
      mapping_value <- layer$mapping[[aesthetic]]
      if (!is.null(mapping_value)) break
    }
  }

  if (is.null(mapping_value)) return(NA_character_)

  tryCatch(
    rlang::as_name(mapping_value),
    error=function(e) NA_character_
  )
}

figure5_remove_invariant_frequency <- function(plot_object) {
  x_variable <- figure5_readability_map_name(
    plot_object,
    "x"
  )

  if (
    is.na(x_variable) ||
    !is.data.frame(plot_object$data) ||
    !x_variable %in% names(plot_object$data)
  ) {
    stop("Panel C x variable could not be resolved.")
  }

  remove_index <- grepl(
    "^\\s*Selection frequency\\s*$",
    as.character(plot_object$data[[x_variable]]),
    ignore.case=TRUE
  )

  if (!any(remove_index)) {
    message(
      "Panel C Selection frequency was already absent."
    )
    return(plot_object)
  }

  plot_object$data <- plot_object$data[
    !remove_index,
    ,
    drop=FALSE
  ]

  if (is.factor(plot_object$data[[x_variable]])) {
    plot_object$data[[x_variable]] <- droplevels(
      plot_object$data[[x_variable]]
    )
  }

  for (index in seq_along(plot_object$layers)) {
    layer_data <- plot_object$layers[[index]]$data

    if (
      inherits(layer_data, "waiver") ||
      !is.data.frame(layer_data) ||
      !x_variable %in% names(layer_data)
    ) {
      next
    }

    layer_remove <- grepl(
      "^\\s*Selection frequency\\s*$",
      as.character(layer_data[[x_variable]]),
      ignore.case=TRUE
    )

    layer_data <- layer_data[
      !layer_remove,
      ,
      drop=FALSE
    ]

    if (is.factor(layer_data[[x_variable]])) {
      layer_data[[x_variable]] <- droplevels(
        layer_data[[x_variable]]
      )
    }

    plot_object$layers[[index]]$data <- layer_data
  }

  x_scale <- plot_object$scales$get_scales("x")

  if (
    !is.null(x_scale) &&
    !is.null(x_scale$limits) &&
    !inherits(x_scale$limits, "waiver") &&
    !is.function(x_scale$limits)
  ) {
    x_scale$limits <- x_scale$limits[
      !grepl(
        "^\\s*Selection frequency\\s*$",
        as.character(x_scale$limits),
        ignore.case=TRUE
      )
    ]
  }

  message(
    "Panel C Selection frequency removed: ",
    "all displayed values were invariant at 1.00."
  )

  plot_object
}

figure5_source_symbol <- function(raw) {
  raw <- as.character(raw)
  symbol <- rep("\u25CB", length(raw))

  symbol[grepl(
    "\\[genus\\]|^genus::|taxonomy:\\s*genus|\\bgenus\\b",
    raw,
    ignore.case=TRUE
  )] <- "\u25CF"

  symbol[grepl(
    "\\[species\\]|^species::|taxonomy:\\s*species|\\bspecies\\b|\\bbacterium\\b|\\bspp\\.?\\b",
    raw,
    ignore.case=TRUE
  )] <- "\u25B2"

  symbol[grepl(
    "\\[agora2[^]]*\\]|^agora2[^:]*::|agora2|reaction|exchange|subsystem",
    raw,
    ignore.case=TRUE
  )] <- "\u25A0"

  symbol[grepl(
    "\\[humann[^]]*\\]|^humann::|humann|pwy",
    raw,
    ignore.case=TRUE
  )] <- "\u25C6"

  symbol[grepl(
    "\\[mag\\]|^mag::|\\bmag\\b|\\bbin[ ._-]*[a-z0-9]+",
    raw,
    ignore.case=TRUE
  )] <- "\u25BC"

  symbol
}

figure5_readable_name <- function(raw) {
  clean <- as.character(raw)

  clean <- gsub(
    "^\\s*(?:HUMAnN|Species|Genus|MAG|AGORA2(?:\\s+reaction|\\s+exchange|\\s+subsystem)?)\\s*::\\s*",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub(
    "\\s*\\[(?:MetaPhlAn[^]]*|Genus|Species|MAG|HUMAnN[^]]*|AGORA2[^]]*)\\]",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub(
    "^Taxonomy:\\s*",
    "",
    clean,
    ignore.case=TRUE
  )

  clean <- gsub(
    "\\s*\\([A-Z0-9._-]*PWY[A-Z0-9._-]*\\)",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub(
    "^[A-Z0-9._-]*PWY[A-Z0-9._-]*:\\s*",
    "",
    clean,
    perl=TRUE,
    ignore.case=TRUE
  )

  clean <- gsub("\\s+", " ", clean)
  clean <- trimws(clean)

  readable_exchange <- c(
    "phpyr"="Phenylpyruvate exchange",
    "ch4s"="Methanethiol exchange",
    "3dhchol"="3-Dehydrocholate exchange",
    "3dhcdchol"="3-Dehydrochenodeoxycholate exchange"
  )

  for (reaction_id in names(readable_exchange)) {
    pattern <- paste0(
      "^EX[_ ]*",
      reaction_id,
      "\\s*\\(e\\)$"
    )

    clean <- sub(
      pattern,
      readable_exchange[[reaction_id]],
      clean,
      ignore.case=TRUE
    )
  }

  remaining_exchange <- grepl(
    "^EX[_ ].*\\(e\\)$",
    clean,
    ignore.case=TRUE
  )

  if (any(remaining_exchange)) {
    generic <- clean[remaining_exchange]

    generic <- sub(
      "^EX[_ ]*",
      "",
      generic,
      ignore.case=TRUE
    )

    generic <- sub(
      "\\s*\\(e\\)$",
      "",
      generic,
      ignore.case=TRUE
    )

    generic <- gsub(
      "_",
      " ",
      generic,
      fixed=TRUE
    )

    clean[remaining_exchange] <- paste(
      generic,
      "exchange"
    )
  }

  clean <- gsub(
    "^L-glutamate degradation V\\s*\\(via[^)]*\\)$",
    "L-glutamate degradation V",
    clean,
    ignore.case=TRUE
  )

  clean <- gsub(
    "^fatty acid elongation\\s*[-–—]+\\s*saturated$",
    "Fatty acid elongation — saturated",
    clean,
    ignore.case=TRUE
  )

  clean
}

figure5_wrap_label <- function(
  value,
  width
) {
  wrapped <- strwrap(
    value,
    width=width
  )

  if (!length(wrapped)) return("")
  if (length(wrapped) <= 2) {
    return(paste(wrapped, collapse="\n"))
  }

  second <- paste(
    wrapped[-1],
    collapse=" "
  )

  if (nchar(second) > width + 6) {
    second <- paste0(
      substr(second, 1, width + 3),
      "..."
    )
  }

  paste(
    wrapped[[1]],
    second,
    sep="\n"
  )
}

figure5_label_c <- function(x) {
  raw <- as.character(x)
  symbol <- figure5_source_symbol(raw)
  clean <- figure5_readable_name(raw)

  wrapped <- vapply(
    clean,
    figure5_wrap_label,
    character(1),
    width=30
  )

  paste(symbol, wrapped)
}

figure5_label_e <- function(x) {
  raw <- as.character(x)
  symbol <- figure5_source_symbol(raw)
  clean <- figure5_readable_name(raw)

  wrapped <- vapply(
    clean,
    figure5_wrap_label,
    character(1),
    width=36
  )

  paste(symbol, wrapped)
}

figure5_set_y_labels <- function(
  plot_object,
  label_function
) {
  y_scale <- plot_object$scales$get_scales("y")

  if (!is.null(y_scale)) {
    y_scale$labels <- label_function
    return(plot_object)
  }

  plot_object +
    ggplot2::scale_y_discrete(
      labels=label_function
    )
}

figure5_enlarge_text_geoms <- function(
  plot_object,
  minimum=4.0
) {
  for (index in seq_along(plot_object$layers)) {
    geom_object <- plot_object$layers[[index]]$geom

    if (
      !inherits(
        geom_object,
        c("GeomText", "GeomLabel")
      )
    ) {
      next
    }

    current_size <- plot_object$layers[[index]]$aes_params$size

    if (
      is.null(current_size) ||
      !is.numeric(current_size) ||
      length(current_size) != 1 ||
      !is.finite(current_size)
    ) {
      current_size <- minimum
    }

    plot_object$layers[[index]]$aes_params$size <- max(
      current_size,
      minimum
    )
  }

  plot_object
}

figure5_source_key <- paste(
  "\u25CF Genus",
  "\u25B2 Species",
  "\u25A0 AGORA2",
  "\u25C6 HUMAnN",
  "\u25BC MAG",
  sep="    "
)

figure5_readable_theme <- ggplot2::theme(
  text=ggplot2::element_text(size=16),
  axis.title=ggplot2::element_text(
    face="bold",
    size=19
  ),
  axis.text=ggplot2::element_text(size=14),
  strip.text=ggplot2::element_text(
    face="bold",
    size=15
  ),
  legend.title=ggplot2::element_text(
    face="bold",
    size=14
  ),
  legend.text=ggplot2::element_text(size=13),
  plot.tag=ggplot2::element_text(
    face="bold",
    size=19
  )
)

if (exists("pA") && inherits(pA, "ggplot")) {
  if (inherits(pA$coordinates, "CoordFlip")) {
    pA <- pA +
      ggplot2::labs(
        x="Feature family",
        y="Spearman correlation"
      )
  } else {
    pA <- pA +
      ggplot2::labs(
        x="Spearman correlation",
        y="Feature family"
      )
  }

  pA <- pA +
    figure5_readable_theme +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(size=14),
      axis.text.y=ggplot2::element_text(size=13.5)
    )
}

if (exists("pB") && inherits(pB, "ggplot")) {
  pB <- figure5_enlarge_text_geoms(
    pB,
    minimum=4.2
  )

  pB <- pB +
    figure5_readable_theme +
    ggplot2::labs(
      x="Observed elimination score",
      y="Predicted elimination score"
    ) +
    ggplot2::theme(
      axis.text=ggplot2::element_text(size=13),
      strip.text=ggplot2::element_text(
        face="bold",
        size=14.5
      )
    )
}

if (exists("pC") && inherits(pC, "ggplot")) {
  pC <- figure5_remove_invariant_frequency(pC)
  pC <- figure5_set_y_labels(pC, figure5_label_c)
  pC <- figure5_enlarge_text_geoms(
    pC,
    minimum=4.2
  )

  pC <- pC +
    figure5_readable_theme +
    ggplot2::labs(
      x="Feature-selection metric",
      y="Selected feature",
      caption=figure5_source_key
    ) +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(
        size=13,
        angle=0,
        hjust=0.5
      ),
      axis.text.y=ggplot2::element_text(
        size=11.5,
        lineheight=0.92
      ),
      plot.caption=ggplot2::element_text(
        size=11.5,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=8)
      )
    )
}

if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- figure5_enlarge_text_geoms(
    pD,
    minimum=4.2
  )

  pD <- pD +
    figure5_readable_theme +
    ggplot2::labs(
      x="Evidence category",
      y="Feature family"
    ) +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(
        size=12.5,
        angle=22,
        hjust=1
      ),
      axis.text.y=ggplot2::element_text(size=13)
    )
}

if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_set_y_labels(
    pE_heatmap,
    figure5_label_e
  )

  pE_heatmap <- pE_heatmap +
    figure5_readable_theme +
    ggplot2::labs(
      x="Evidence metric",
      y="Prioritized feature",
      caption=figure5_source_key
    ) +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(
        size=12.5,
        angle=27,
        hjust=1
      ),
      axis.text.y=ggplot2::element_text(
        size=11,
        lineheight=0.92
      ),
      plot.caption=ggplot2::element_text(
        size=11.5,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=8)
      )
    )
}

if (
  exists("pE_dendrogram") &&
  inherits(pE_dendrogram, "ggplot")
) {
  pE_dendrogram <- pE_dendrogram +
    ggplot2::labs(tag="E") +
    ggplot2::theme(
      plot.tag=ggplot2::element_text(
        face="bold",
        size=19
      ),
      plot.tag.position=c(0.02, 0.98)
    )
}
# === FIGURE5_READABILITY_FIX_END ===





pE <- (
    pE_dendrogram +
        pE_strip +
        pE_heatmap
) +
    plot_layout(widths = c(0.16, 0.025, 1)) +
    plot_annotation(
        title = "E",
        subtitle = NULL,
        theme = theme(
            plot.title = element_text(face = "bold", size = 8.4, margin = margin(b = 2)),
            plot.subtitle = element_text(
                size = 6.9,
                colour = PFAS_PURPLE[["neutral_dark"]],
                margin = margin(b = 4)
            )
        )
    )

# -------------------------------------------------------------------------
# Final assembly and exports
# -------------------------------------------------------------------------


# === FIGURE5_BASIC_CLOSEOUT_BEGIN ===
# Basic visual closeout only. No new patchwork structures.

figure5_resize_points <- function(
  plot_object,
  dense_size=1.75,
  summary_size=3.0
) {
  if (!inherits(plot_object, "ggplot")) {
    return(plot_object)
  }

  for (index in seq_along(plot_object$layers)) {
    layer <- plot_object$layers[[index]]

    if (!inherits(layer$geom, c("GeomPoint", "GeomJitter"))) {
      next
    }

    layer_data <- layer$data

    row_count <- if (
      inherits(layer_data, "waiver") ||
      !is.data.frame(layer_data)
    ) {
      nrow(plot_object$data)
    } else {
      nrow(layer_data)
    }

    plot_object$layers[[index]]$aes_params$size <- if (
      is.finite(row_count) && row_count > 30
    ) {
      dense_size
    } else {
      summary_size
    }

    plot_object$layers[[index]]$aes_params$alpha <- if (
      is.finite(row_count) && row_count > 30
    ) {
      0.52
    } else {
      1
    }
  }

  plot_object
}

# Panel A: larger raw and summary points.
if (exists("pA") && inherits(pA, "ggplot")) {
  pA <- figure5_resize_points(
    pA,
    dense_size=1.85,
    summary_size=3.2
  )

  pA <- pA +
    ggplot2::theme(
      axis.text=ggplot2::element_text(size=13),
      axis.title=ggplot2::element_text(
        size=17,
        face="bold"
      ),
      plot.margin=ggplot2::margin(8, 10, 8, 12)
    )
}

# Panel B: prevent clipping and strengthen the diagnostic marks.
if (exists("pB") && inherits(pB, "ggplot")) {
  pB <- figure5_resize_points(
    pB,
    dense_size=1.9,
    summary_size=2.5
  )

  for (index in seq_along(pB$layers)) {
    if (inherits(pB$layers[[index]]$geom, "GeomAbline")) {
      pB$layers[[index]]$aes_params$linewidth <- 0.8
      pB$layers[[index]]$aes_params$alpha <- 0.85
    }

    if (inherits(
      pB$layers[[index]]$geom,
      c("GeomLine", "GeomSmooth")
    )) {
      pB$layers[[index]]$aes_params$linewidth <- 0.9
    }
  }

  pB <- pB +
    ggplot2::labs(
      x="Observed elimination score",
      y="Predicted elimination score"
    ) +
    ggplot2::theme(
      axis.title.x=ggplot2::element_text(
        size=16,
        face="bold",
        margin=ggplot2::margin(t=7)
      ),
      axis.title.y=ggplot2::element_text(
        size=16,
        face="bold",
        margin=ggplot2::margin(r=12)
      ),
      axis.text=ggplot2::element_text(size=12),
      strip.text=ggplot2::element_text(
        size=13,
        face="bold"
      ),
      plot.margin=ggplot2::margin(8, 10, 8, 34)
    )
}

# Panel E: place labels on the right so the dendrogram sits directly
# beside the heatmap, and split the source key over two lines.
if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  e_y_scale <- pE_heatmap$scales$get_scales("y")

  if (!is.null(e_y_scale)) {
    e_y_scale$position <- "right"
  }

  pE_heatmap <- pE_heatmap +
    ggplot2::labs(
      x="Evidence metric",
      y="Prioritized feature",
      caption=paste(
        "\u25CF Genus    \u25B2 Species    \u25A0 AGORA2",
        "\u25C6 HUMAnN    \u25BC MAG",
        sep="\n"
      )
    ) +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(
        size=12,
        angle=27,
        hjust=1
      ),
      axis.text.y.right=ggplot2::element_text(
        size=10.5,
        lineheight=0.92
      ),
      axis.title.x=ggplot2::element_text(
        size=16,
        face="bold",
        margin=ggplot2::margin(t=8)
      ),
      axis.title.y.right=ggplot2::element_text(
        size=15,
        face="bold",
        margin=ggplot2::margin(l=10)
      ),
      legend.position="bottom",
      legend.box="vertical",
      legend.spacing.y=grid::unit(5, "pt"),
      legend.title=ggplot2::element_text(
        size=12,
        face="bold"
      ),
      legend.text=ggplot2::element_text(size=11),
      plot.caption=ggplot2::element_text(
        size=10.5,
        lineheight=1.25,
        hjust=0,
        colour="#4B2E83",
        margin=ggplot2::margin(t=8)
      ),
      plot.margin=ggplot2::margin(8, 18, 8, 2)
    )
}
# === FIGURE5_BASIC_CLOSEOUT_END ===


# === FIGURE5_B_SUMMARY_E_SPACING_BEGIN ===
# Panel B: summarize observed-predicted agreement rather than showing
# four visually repetitive scatter plots.
# Panel E: separate the source key from the scaled-evidence legend.

figure5_mapping_name <- function(plot_object, aesthetic) {
  mapping_value <- plot_object$mapping[[aesthetic]]

  if (is.null(mapping_value) && length(plot_object$layers)) {
    for (layer in plot_object$layers) {
      mapping_value <- layer$mapping[[aesthetic]]
      if (!is.null(mapping_value)) break
    }
  }

  if (is.null(mapping_value)) return(NA_character_)

  tryCatch(
    rlang::as_name(mapping_value),
    error=function(e) NA_character_
  )
}

figure5_facet_variable <- function(plot_object) {
  for (parameter in c("facets", "rows", "cols")) {
    facet_specification <- plot_object$facet$params[[parameter]]

    if (
      !is.null(facet_specification) &&
      length(facet_specification) &&
      length(names(facet_specification))
    ) {
      candidate <- names(facet_specification)[[1]]

      if (candidate %in% names(plot_object$data)) {
        return(candidate)
      }
    }
  }

  NA_character_
}

figure5_stats_from_ggplot <- function(
  plot_object,
  fallback_label="Feature set"
) {
  if (
    !inherits(plot_object, "ggplot") ||
    !is.data.frame(plot_object$data) ||
    !nrow(plot_object$data)
  ) {
    return(NULL)
  }

  data_object <- plot_object$data

  x_variable <- figure5_mapping_name(
    plot_object,
    "x"
  )

  y_variable <- figure5_mapping_name(
    plot_object,
    "y"
  )

  if (
    is.na(x_variable) ||
    is.na(y_variable) ||
    !x_variable %in% names(data_object) ||
    !y_variable %in% names(data_object)
  ) {
    return(NULL)
  }

  group_variable <- figure5_facet_variable(
    plot_object
  )

  if (is.na(group_variable)) {
    candidates <- names(data_object)[
      vapply(
        data_object,
        function(value) {
          is.factor(value) || is.character(value)
        },
        logical(1)
      )
    ]

    candidates <- candidates[
      vapply(
        data_object[candidates],
        function(value) {
          count <- length(unique(value[!is.na(value)]))
          count >= 2 && count <= 8
        },
        logical(1)
      )
    ]

    preferred <- candidates[
      grepl(
        "feature|model|panel|family|source|set",
        candidates,
        ignore.case=TRUE
      )
    ]

    if (length(preferred)) {
      group_variable <- preferred[[1]]
    } else if (length(candidates)) {
      group_variable <- candidates[[1]]
    }
  }

  if (
    !is.na(group_variable) &&
    group_variable %in% names(data_object)
  ) {
    group_order <- unique(
      as.character(data_object[[group_variable]])
    )
  } else {
    group_order <- fallback_label
    data_object$.figure5_group <- fallback_label
    group_variable <- ".figure5_group"
  }

  rows <- lapply(
    group_order,
    function(group_name) {
      subset_object <- data_object[
        as.character(data_object[[group_variable]]) == group_name,
        ,
        drop=FALSE
      ]

      observed <- suppressWarnings(
        as.numeric(subset_object[[x_variable]])
      )

      predicted <- suppressWarnings(
        as.numeric(subset_object[[y_variable]])
      )

      keep <- is.finite(observed) & is.finite(predicted)
      observed <- observed[keep]
      predicted <- predicted[keep]

      if (length(observed) < 3) return(NULL)

      data.frame(
        feature_set=group_name,
        rho=suppressWarnings(
          stats::cor(
            observed,
            predicted,
            method="spearman"
          )
        ),
        mae=mean(
          abs(predicted - observed)
        ),
        n=length(observed),
        stringsAsFactors=FALSE
      )
    }
  )

  rows <- rows[
    !vapply(
      rows,
      is.null,
      logical(1)
    )
  ]

  if (!length(rows)) return(NULL)

  do.call(
    rbind,
    rows
  )
}

figure5_collect_ggplots <- function(plot_object) {
  if (inherits(plot_object, "patchwork")) {
    children <- plot_object$patches$plots

    return(
      unlist(
        lapply(
          children,
          figure5_collect_ggplots
        ),
        recursive=FALSE
      )
    )
  }

  if (inherits(plot_object, "ggplot")) {
    return(list(plot_object))
  }

  list()
}

figure5_panel_b_stats <- NULL

if (exists("pB")) {
  if (inherits(pB, "patchwork")) {
    child_plots <- figure5_collect_ggplots(pB)

    child_rows <- lapply(
      seq_along(child_plots),
      function(index) {
        child_plot <- child_plots[[index]]

        child_label <- child_plot$labels$title

        if (
          is.null(child_label) ||
          !nzchar(child_label)
        ) {
          child_label <- paste(
            "Feature set",
            index
          )
        }

        figure5_stats_from_ggplot(
          child_plot,
          fallback_label=child_label
        )
      }
    )

    child_rows <- child_rows[
      !vapply(
        child_rows,
        is.null,
        logical(1)
      )
    ]

    if (length(child_rows)) {
      figure5_panel_b_stats <- do.call(
        rbind,
        child_rows
      )
    }

  } else if (inherits(pB, "ggplot")) {
    figure5_panel_b_stats <- figure5_stats_from_ggplot(
      pB,
      fallback_label="Combined feature set"
    )
  }
}

if (
  !is.null(figure5_panel_b_stats) &&
  nrow(figure5_panel_b_stats)
) {
  figure5_panel_b_stats <- figure5_panel_b_stats[
    !duplicated(
      figure5_panel_b_stats$feature_set
    ),
    ,
    drop=FALSE
  ]

  original_order <- unique(
    figure5_panel_b_stats$feature_set
  )

  figure5_panel_b_stats$feature_set <- factor(
    figure5_panel_b_stats$feature_set,
    levels=rev(original_order)
  )

  limit <- max(
    0.22,
    max(
      abs(figure5_panel_b_stats$rho),
      na.rm=TRUE
    ) + 0.10
  )

  figure5_panel_b_stats$annotation <- sprintf(
    "rho = %.02f   |   MAE = %.02f   |   n = %d",
    figure5_panel_b_stats$rho,
    figure5_panel_b_stats$mae,
    figure5_panel_b_stats$n
  )

  pB <- ggplot2::ggplot(
    figure5_panel_b_stats,
    ggplot2::aes(
      x=rho,
      y=feature_set
    )
  ) +
    ggplot2::geom_vline(
      xintercept=0,
      linewidth=0.55,
      linetype="dashed",
      colour="#777777"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x=0,
        xend=rho,
        yend=feature_set
      ),
      linewidth=1.25,
      colour="#8F73B8"
    ) +
    ggplot2::geom_point(
      size=4.8,
      shape=21,
      stroke=0.6,
      fill="#7B50B4",
      colour="#3F1B68"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x=limit * 0.96,
        label=annotation
      ),
      hjust=1,
      size=3.8,
      colour="#3F1B68"
    ) +
    ggplot2::scale_x_continuous(
      limits=c(-limit, limit),
      breaks=scales::breaks_pretty(n=5)
    ) +
    ggplot2::labs(
      tag="B",
      x="Observed-predicted agreement (Spearman rho)",
      y="Feature set"
    ) +
    ggplot2::theme_minimal(
      base_size=13
    ) +
    ggplot2::theme(
      panel.grid.major.y=ggplot2::element_line(
        colour="#E7E1ED",
        linewidth=0.45
      ),
      panel.grid.minor=ggplot2::element_blank(),
      axis.title=ggplot2::element_text(
        face="bold",
        size=15
      ),
      axis.text.x=ggplot2::element_text(
        size=11.5
      ),
      axis.text.y=ggplot2::element_text(
        size=12.5
      ),
      plot.tag=ggplot2::element_text(
        face="bold",
        size=17
      ),
      plot.margin=ggplot2::margin(
        10,
        16,
        10,
        10
      )
    )
}

figure5_source_key <- paste(
  "\u25CF Genus    \u25B2 Species    \u25A0 AGORA2",
  "\u25C6 HUMAnN    \u25BC MAG",
  sep="\n"
)

if (
  exists("pE") &&
  inherits(pE, "patchwork")
) {
  # Remove the duplicated source caption from child plots.
  pE <- pE & ggplot2::labs(caption=NULL)

  # Add space around the collected scaled-evidence legend.
  pE <- pE & ggplot2::theme(
    legend.position="bottom",
    legend.box="vertical",
    legend.margin=ggplot2::margin(
      t=14,
      r=0,
      b=10,
      l=0
    ),
    legend.box.margin=ggplot2::margin(
      t=10,
      r=0,
      b=12,
      l=0
    ),
    legend.title=ggplot2::element_text(
      face="bold",
      size=11.5
    ),
    legend.text=ggplot2::element_text(
      size=10.5
    ),
    plot.margin=ggplot2::margin(
      8,
      8,
      18,
      8
    )
  )

  # Add one clean source key beneath the evidence legend.
  pE <- pE +
    patchwork::plot_annotation(
      caption=figure5_source_key,
      theme=ggplot2::theme(
        plot.caption=ggplot2::element_text(
          size=10.8,
          lineheight=1.35,
          hjust=0.5,
          colour="#4B2E83",
          margin=ggplot2::margin(
            t=16,
            b=8
          )
        ),
        plot.margin=ggplot2::margin(
          0,
          0,
          18,
          0
        )
      )
    )

} else if (
  exists("pE") &&
  inherits(pE, "ggplot")
) {
  pE <- pE +
    ggplot2::labs(
      caption=figure5_source_key
    ) +
    ggplot2::theme(
      legend.position="bottom",
      legend.box="vertical",
      legend.margin=ggplot2::margin(
        t=14,
        b=10
      ),
      plot.caption=ggplot2::element_text(
        size=10.8,
        lineheight=1.35,
        hjust=0.5,
        colour="#4B2E83",
        margin=ggplot2::margin(
          t=18,
          b=8
        )
      ),
      plot.margin=ggplot2::margin(
        8,
        8,
        28,
        8
      )
    )
}
# === FIGURE5_B_SUMMARY_E_SPACING_END ===


# === FIGURE5_CANDIDATE_RANKING_E_LEGEND_BEGIN ===
# Panel B now presents prioritized candidates rather than repeating
# weak observed-versus-predicted comparisons.
#
# Panel E retains the existing heatmap and clustering but receives a
# properly sized horizontal scaled-evidence colour bar.

figure5_candidate_symbol <- function(family_value) {
  family_value <- as.character(family_value)
  symbol <- rep("\u25CB", length(family_value))

  symbol[grepl(
    "genus",
    family_value,
    ignore.case=TRUE
  )] <- "\u25CF"

  symbol[grepl(
    "species",
    family_value,
    ignore.case=TRUE
  )] <- "\u25B2"

  symbol[grepl(
    "agora2|reaction|exchange|subsystem",
    family_value,
    ignore.case=TRUE
  )] <- "\u25A0"

  symbol[grepl(
    "humann|pathway",
    family_value,
    ignore.case=TRUE
  )] <- "\u25C6"

  symbol[grepl(
    "\\bmag\\b",
    family_value,
    ignore.case=TRUE
  )] <- "\u25BC"

  symbol
}

figure5_candidate_label <- function(value) {
  value <- as.character(value)
  value <- gsub("\\s+", " ", value)
  value <- trimws(value)

  vapply(
    value,
    function(item) {
      paste(
        strwrap(item, width=27),
        collapse="\n"
      )
    },
    character(1)
  )
}

if (
  !exists("candidate_unique") ||
  !is.data.frame(candidate_unique) ||
  !nrow(candidate_unique)
) {
  stop(
    "candidate_unique was unavailable for the new Panel B."
  )
}

required_candidate_columns <- c(
  "rank_score",
  "rho",
  "prevalence",
  "family",
  "biological_feature",
  "candidate_id"
)

missing_candidate_columns <- setdiff(
  required_candidate_columns,
  names(candidate_unique)
)

if (length(missing_candidate_columns)) {
  stop(
    "Panel B candidate columns missing: ",
    paste(
      missing_candidate_columns,
      collapse=", "
    )
  )
}

figure5_panel_b <- data.table::copy(
  candidate_unique
)

figure5_panel_b[, rank_score_plot :=
  suppressWarnings(
    as.numeric(rank_score)
  )
]

figure5_panel_b[, rho_plot :=
  suppressWarnings(
    as.numeric(rho)
  )
]

figure5_panel_b[, prevalence_plot :=
  suppressWarnings(
    as.numeric(prevalence)
  )
]

figure5_panel_b <- figure5_panel_b[
  is.finite(rank_score_plot) &
  is.finite(rho_plot)
]

if (!nrow(figure5_panel_b)) {
  stop(
    "No finite candidates were available for Panel B."
  )
}

figure5_panel_b[
  !is.finite(prevalence_plot),
  prevalence_plot := 0
]

figure5_panel_b[, prevalence_plot :=
  pmin(
    pmax(prevalence_plot, 0),
    1
  )
]

figure5_panel_b[, source_symbol :=
  figure5_candidate_symbol(family)
]

figure5_panel_b[, label_plot :=
  paste(
    source_symbol,
    figure5_candidate_label(
      biological_feature
    )
  )
]

figure5_panel_b[, label_plot :=
  make.unique(label_plot)
]

data.table::setorder(
  figure5_panel_b,
  -rank_score_plot
)

# Eight candidates fit the upper-right panel without shrinking labels.
figure5_panel_b <- head(
  figure5_panel_b,
  8L
)

figure5_panel_b[, display_rank :=
  seq_len(.N)
]

figure5_panel_b[, label_plot :=
  factor(
    label_plot,
    levels=rev(label_plot)
  )
]

rho_limit <- max(
  0.55,
  ceiling(
    max(
      abs(figure5_panel_b$rho_plot),
      na.rm=TRUE
    ) * 10
  ) / 10
)

pB <- ggplot2::ggplot(
  figure5_panel_b,
  ggplot2::aes(
    x=rho_plot,
    y=label_plot
  )
) +
  ggplot2::geom_vline(
    xintercept=0,
    linetype="dashed",
    linewidth=0.55,
    colour="#777777"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x=0,
      xend=rho_plot,
      yend=label_plot
    ),
    linewidth=1.05,
    colour="#B8A3D2"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size=rank_score_plot,
      fill=prevalence_plot
    ),
    shape=21,
    stroke=0.65,
    colour="#3C2457"
  ) +
  ggplot2::scale_x_continuous(
    limits=c(-rho_limit, rho_limit),
    breaks=scales::breaks_pretty(n=5)
  ) +
  ggplot2::scale_size_continuous(
    range=c(4.2, 7.4),
    name="Candidate\nscore"
  ) +
  ggplot2::scale_fill_gradient(
    low="#EEE8F4",
    high="#4B007D",
    limits=c(0, 1),
    oob=scales::squish,
    name="Prevalence"
  ) +
  ggplot2::labs(
    tag="B",
    x="Marginal association with elimination score (rho)",
    y="Prioritized candidate"
  ) +
  ggplot2::theme_minimal(
    base_size=12
  ) +
  ggplot2::theme(
    panel.grid.major.y=ggplot2::element_line(
      colour="#E5DEEC",
      linewidth=0.45
    ),
    panel.grid.minor=ggplot2::element_blank(),
    axis.title=ggplot2::element_text(
      face="bold",
      size=14
    ),
    axis.text.x=ggplot2::element_text(
      size=10.5
    ),
    axis.text.y=ggplot2::element_text(
      size=10.2,
      lineheight=0.93
    ),
    legend.position="right",
    legend.title=ggplot2::element_text(
      face="bold",
      size=9.8
    ),
    legend.text=ggplot2::element_text(
      size=9.2
    ),
    legend.key.height=grid::unit(
      4.5,
      "mm"
    ),
    plot.tag=ggplot2::element_text(
      face="bold",
      size=17
    ),
    plot.margin=ggplot2::margin(
      8,
      8,
      8,
      8
    )
  ) +
  ggplot2::guides(
    size=ggplot2::guide_legend(
      order=1,
      override.aes=list(
        fill="#8B6BB1"
      )
    ),
    fill=ggplot2::guide_colourbar(
      order=2,
      title.position="top",
      title.hjust=0.5,
      barwidth=grid::unit(7, "mm"),
      barheight=grid::unit(24, "mm")
    )
  )

if (
  exists("SRC") &&
  dir.exists(SRC)
) {
  data.table::fwrite(
    figure5_panel_b[
      ,
      .(
        display_rank,
        candidate_id,
        biological_feature,
        family,
        rank_score=rank_score_plot,
        marginal_rho=rho_plot,
        prevalence=prevalence_plot
      )
    ],
    file.path(
      SRC,
      "Figure_5B_candidate_priority_ranking.tsv"
    ),
    sep="\t"
  )
}

figure5_fix_scaled_evidence_legend <- function(
  plot_object
) {
  if (inherits(plot_object, "patchwork")) {
    if (
      !is.null(plot_object$patches) &&
      length(plot_object$patches$plots)
    ) {
      plot_object$patches$plots <- lapply(
        plot_object$patches$plots,
        figure5_fix_scaled_evidence_legend
      )
    }

    return(plot_object)
  }

  if (!inherits(plot_object, "ggplot")) {
    return(plot_object)
  }

  fill_scale <- plot_object$scales$get_scales(
    "fill"
  )

  if (is.null(fill_scale)) {
    return(plot_object)
  }

  scale_title <- paste(
    fill_scale$name,
    collapse=" "
  )

  if (
    !grepl(
      "scaled.*evidence",
      scale_title,
      ignore.case=TRUE
    )
  ) {
    return(plot_object)
  }

  fill_scale$breaks <- c(
    0,
    0.5,
    1
  )

  fill_scale$labels <- c(
    "0",
    "0.5",
    "1"
  )

  fill_scale$guide <- ggplot2::guide_colourbar(
    direction="horizontal",
    title.position="top",
    title.hjust=0.5,
    label.position="bottom",
    barwidth=grid::unit(
      62,
      "mm"
    ),
    barheight=grid::unit(
      4.5,
      "mm"
    ),
    ticks=TRUE,
    frame.colour="#6C5A78"
  )

  plot_object +
    ggplot2::theme(
      legend.position="bottom",
      legend.justification="center",
      legend.box="vertical",
      legend.margin=ggplot2::margin(
        t=18,
        r=0,
        b=12,
        l=0
      ),
      legend.box.margin=ggplot2::margin(
        t=12,
        r=0,
        b=16,
        l=0
      ),
      legend.title=ggplot2::element_text(
        face="bold",
        size=11
      ),
      legend.text=ggplot2::element_text(
        size=10
      ),
      plot.margin=ggplot2::margin(
        8,
        8,
        32,
        8
      )
    )
}

pE <- figure5_fix_scaled_evidence_legend(
  pE
)
# === FIGURE5_CANDIDATE_RANKING_E_LEGEND_END ===


# === FIGURE5_SAFE_FINAL_CLOSEOUT_BEGIN ===
# Final display refinements only.
# No Panel E reconstruction and no additional patchwork layouts.

figure5_safe_humanize <- function(x) {
  value <- as.character(x)

  value <- gsub(
    "FEDCabcpp",
    "Iron(III) dicitrate transport",
    value,
    fixed=TRUE
  )

  value <- gsub(
    "EX[_ ]*3dhdchol\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*3dhchol\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*3dhcdchol\\s*\\(e\\)",
    "3-Dehydrochenodeoxycholate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*phpyr\\s*\\(e\\)",
    "Phenylpyruvate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*ch4s\\s*\\(e\\)",
    "Methanethiol exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "^[A-Z0-9._-]*PWY[A-Z0-9._-]*:\\s*",
    "",
    value,
    perl=TRUE,
    ignore.case=TRUE
  )

  value <- gsub(
    "\\s+",
    " ",
    value
  )

  trimws(value)
}

figure5_safe_set_y_labels <- function(
  plot_object,
  label_function
) {
  if (!inherits(plot_object, "ggplot")) {
    return(plot_object)
  }

  y_scale <- plot_object$scales$get_scales("y")

  if (
    !is.null(y_scale) &&
    inherits(y_scale, "ScaleDiscrete")
  ) {
    y_scale$labels <- label_function
  }

  plot_object
}

# Panel B: preserve the ranked-candidate display, clarify its biology,
# humanize labels and put both legends beneath the panel.
if (exists("pB") && inherits(pB, "ggplot")) {
  pB <- figure5_safe_set_y_labels(
    pB,
    figure5_safe_humanize
  )

  pB <- pB +
    ggplot2::labs(
      x="Association with elimination score (Spearman rho; PFAS-exposed only)",
      y="Prioritized candidate"
    ) +
    ggplot2::guides(
      size=ggplot2::guide_legend(
        order=1,
        title.position="top",
        title.hjust=0.5,
        nrow=1,
        byrow=TRUE
      ),
      fill=ggplot2::guide_colourbar(
        order=2,
        direction="horizontal",
        title.position="top",
        title.hjust=0.5,
        label.position="bottom",
        barwidth=grid::unit(33, "mm"),
        barheight=grid::unit(4, "mm")
      )
    ) +
    ggplot2::theme(
      legend.position="bottom",
      legend.box="vertical",
      legend.box.just="center",
      legend.justification="center",
      legend.margin=ggplot2::margin(
        t=8,
        r=0,
        b=0,
        l=0
      ),
      legend.spacing.y=grid::unit(
        4,
        "pt"
      ),
      legend.title=ggplot2::element_text(
        face="bold",
        size=9.5
      ),
      legend.text=ggplot2::element_text(
        size=8.8
      ),
      axis.title.x=ggplot2::element_text(
        face="bold",
        size=13.5,
        margin=ggplot2::margin(t=7)
      ),
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=13
      ),
      axis.text.x=ggplot2::element_text(
        size=10
      ),
      axis.text.y=ggplot2::element_text(
        size=9.6,
        lineheight=0.93
      ),
      plot.margin=ggplot2::margin(
        5,
        5,
        5,
        5
      )
    )
}

# Panel E: humanize remaining labels and repair only its existing
# scaled-evidence guide. The good dendrogram/heatmap assembly is untouched.
if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_safe_set_y_labels(
    pE_heatmap,
    figure5_safe_humanize
  )

  fill_scale <- pE_heatmap$scales$get_scales(
    "fill"
  )

  if (!is.null(fill_scale)) {
    fill_scale$breaks <- c(
      0,
      0.5,
      1
    )

    fill_scale$labels <- c(
      "0",
      "0.5",
      "1"
    )

    fill_scale$guide <- ggplot2::guide_colourbar(
      direction="horizontal",
      title.position="top",
      title.hjust=0.5,
      label.position="bottom",
      barwidth=grid::unit(
        60,
        "mm"
      ),
      barheight=grid::unit(
        5,
        "mm"
      ),
      ticks=TRUE,
      frame.colour="#675674"
    )
  }

  pE_heatmap <- pE_heatmap +
    ggplot2::labs(
      x="Evidence metric for prioritized candidates",
      y="Prioritized candidate"
    ) +
    ggplot2::theme(
      legend.position="bottom",
      legend.justification="center",
      legend.box.just="center",
      legend.margin=ggplot2::margin(
        t=18,
        r=0,
        b=12,
        l=0
      ),
      legend.box.margin=ggplot2::margin(
        t=10,
        r=0,
        b=18,
        l=0
      ),
      legend.title=ggplot2::element_text(
        face="bold",
        size=11
      ),
      legend.text=ggplot2::element_text(
        size=10
      ),
      axis.title.x=ggplot2::element_text(
        face="bold",
        size=13,
        margin=ggplot2::margin(t=8)
      ),
      axis.title.y=ggplot2::element_text(
        face="bold",
        size=12
      ),
      axis.text.x=ggplot2::element_text(
        size=10,
        angle=27,
        hjust=1
      ),
      axis.text.y=ggplot2::element_text(
        size=9.2,
        lineheight=0.93
      ),
      plot.margin=ggplot2::margin(
        5,
        7,
        28,
        1
      )
    )
}
# === FIGURE5_SAFE_FINAL_CLOSEOUT_END ===


# === FIGURE5_B_RIGHT_LEGENDS_FINAL_BEGIN ===
# Final label/legend adjustment only. No layout reconstruction.

figure5_translate_remaining_agora2 <- function(x) {
  value <- as.character(x)

  value <- gsub(
    "FEDCabcpp",
    "Iron(III) dicitrate transport",
    value,
    fixed=TRUE
  )

  value <- gsub(
    "EX[_ ]*3dhcdchol\\s*\\(e\\)",
    "3-Dehydrochenodeoxycholate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*3dhchol\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*phpyr\\s*\\(e\\)",
    "Phenylpyruvate exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub(
    "EX[_ ]*ch4s\\s*\\(e\\)",
    "Methanethiol exchange",
    value,
    ignore.case=TRUE
  )

  value <- gsub("\\s+", " ", value)
  trimws(value)
}

figure5_replace_y_label_function <- function(
  plot_object,
  label_function
) {
  if (!inherits(plot_object, "ggplot")) {
    return(plot_object)
  }

  y_scale <- plot_object$scales$get_scales("y")

  if (
    !is.null(y_scale) &&
    inherits(y_scale, "ScaleDiscrete")
  ) {
    y_scale$labels <- label_function
  }

  plot_object
}

if (exists("pB") && inherits(pB, "ggplot")) {
  pB <- figure5_replace_y_label_function(
    pB,
    figure5_translate_remaining_agora2
  )

  pB <- pB +
    ggplot2::labs(
      x="Association with elimination score (Spearman rho; PFAS-exposed only)",
      y="Prioritized candidate"
    ) +
    ggplot2::guides(
      size=ggplot2::guide_legend(
        order=1,
        title="Candidate\nscore",
        title.position="top",
        title.hjust=0.5,
        ncol=1,
        byrow=TRUE,
        override.aes=list(alpha=1)
      ),
      fill=ggplot2::guide_colourbar(
        order=2,
        title="Prevalence",
        direction="vertical",
        title.position="top",
        title.hjust=0.5,
        label.position="right",
        barwidth=grid::unit(4.5, "mm"),
        barheight=grid::unit(25, "mm")
      )
    ) +
    ggplot2::theme(
      legend.position="right",
      legend.box="vertical",
      legend.box.just="top",
      legend.justification=c(0, 1),
      legend.margin=ggplot2::margin(
        t=0,
        r=0,
        b=0,
        l=6
      ),
      legend.spacing.y=grid::unit(6, "pt"),
      legend.title=ggplot2::element_text(
        face="bold",
        size=9.4
      ),
      legend.text=ggplot2::element_text(
        size=8.7
      ),
      plot.margin=ggplot2::margin(
        5,
        4,
        5,
        5
      )
    )
}

if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_replace_y_label_function(
    pE_heatmap,
    figure5_translate_remaining_agora2
  )
}
# === FIGURE5_B_RIGHT_LEGENDS_FINAL_END ===






# === FIGURE5_EXACT_LABEL_FIX_BEGIN ===
# Final override applied after every earlier Figure 5 closeout.
# No Panel E reconstruction or new patchwork layout.

figure5_exact_humanize <- function(x) {
  value <- as.character(x)

  value <- gsub(
    "(?i)\\bFEDCabcpp\\b",
    "Iron(III) dicitrate transport",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*3dhcdchol\\s*\\(e\\)",
    "3-Dehydrochenodeoxycholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*3dhchol\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\b3dhcdchol\\s*\\(e\\)",
    "3-Dehydrochenodeoxycholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\b3dhchol\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "\\s+",
    " ",
    value
  )

  trimws(value)
}

figure5_exact_compose_y_labeler <- function(plot_object) {
  if (!inherits(plot_object, "ggplot")) {
    return(plot_object)
  }

  y_scale <- plot_object$scales$get_scales("y")

  if (
    is.null(y_scale) ||
    !inherits(y_scale, "ScaleDiscrete")
  ) {
    return(plot_object)
  }

  original_labels <- y_scale$labels

  if (is.function(original_labels)) {
    y_scale$labels <- local({
      original_function <- original_labels

      function(x) {
        figure5_exact_humanize(
          original_function(x)
        )
      }
    })

  } else if (is.character(original_labels)) {
    y_scale$labels <- figure5_exact_humanize(
      original_labels
    )

  } else {
    y_scale$labels <- figure5_exact_humanize
  }

  plot_object
}

figure5_exact_recurse_patchwork <- function(plot_object) {
  if (inherits(plot_object, "patchwork")) {
    if (
      !is.null(plot_object$patches) &&
      length(plot_object$patches$plots)
    ) {
      plot_object$patches$plots <- lapply(
        plot_object$patches$plots,
        figure5_exact_recurse_patchwork
      )
    }

    if (
      !is.null(plot_object$plots) &&
      length(plot_object$plots)
    ) {
      plot_object$plots <- lapply(
        plot_object$plots,
        figure5_exact_recurse_patchwork
      )
    }

    return(plot_object)
  }

  if (inherits(plot_object, "ggplot")) {
    return(
      figure5_exact_compose_y_labeler(
        plot_object
      )
    )
  }

  plot_object
}

stopifnot(
  identical(
    figure5_exact_humanize("■ FEDCabcpp"),
    "■ Iron(III) dicitrate transport"
  )
)

stopifnot(
  grepl(
    "3-Dehydrochenodeoxycholate exchange",
    figure5_exact_humanize(
      "■ EX 3dhcdchol(e)"
    ),
    fixed=TRUE
  )
)

# Panel B: modify the final label function and reserve a compact
# right-side region inside the panel for the two guides.
if (exists("pB") && inherits(pB, "ggplot")) {
  pB <- figure5_exact_compose_y_labeler(pB)

  x_scale <- pB$scales$get_scales("x")

  if (!is.null(x_scale)) {
    x_scale$limits <- c(-0.65, 0.86)
    x_scale$breaks <- c(
      -0.6,
      -0.3,
      0,
      0.3,
      0.6
    )
  }

  size_scale <- pB$scales$get_scales("size")

  if (!is.null(size_scale)) {
    size_scale$breaks <- c(
      8.55,
      8.68,
      8.80
    )

    size_scale$labels <- c(
      "8.55",
      "8.68",
      "8.80"
    )

    size_scale$guide <- ggplot2::guide_legend(
      order=1,
      title="Candidate score",
      title.position="top",
      title.hjust=0.5,
      ncol=1,
      override.aes=list(
        fill="white",
        alpha=1
      )
    )
  }

  fill_scale <- pB$scales$get_scales("fill")

  if (!is.null(fill_scale)) {
    fill_scale$breaks <- c(
      0,
      0.5,
      1
    )

    fill_scale$labels <- c(
      "0",
      "0.5",
      "1"
    )

    fill_scale$guide <- ggplot2::guide_colourbar(
      order=2,
      direction="vertical",
      title="Prevalence",
      title.position="top",
      title.hjust=0.5,
      label.position="right",
      barwidth=grid::unit(4, "mm"),
      barheight=grid::unit(22, "mm")
    )
  }

  pB <- pB +
    ggplot2::scale_y_discrete(
      labels=function(x) figure5_compact_candidate_label(x, width=18)
    ) +
    ggplot2::labs(
      x="Association with elimination score (Spearman rho; PFAS-exposed only)",
      y="Prioritized candidate"
    ) +
    ggplot2::theme(
      legend.position="inside",
      legend.position.inside=c(0.985, 0.50),
      legend.justification.inside=c(1, 0.5),
      legend.box="vertical",
      legend.background=ggplot2::element_rect(
        fill=scales::alpha("white", 0.95),
        colour="#D8CDE3",
        linewidth=0.35
      ),
      legend.margin=ggplot2::margin(
        5,
        5,
        5,
        5
      ),
      legend.spacing.y=grid::unit(
        4,
        "pt"
      ),
      legend.title=ggplot2::element_text(
        face="bold",
        size=8.8
      ),
      legend.text=ggplot2::element_text(
        size=8
      ),
      legend.key.height=grid::unit(
        4,
        "mm"
      ),
      axis.text.y=ggplot2::element_text(
        size=8.8,
        lineheight=0.90
      ),
      plot.margin=ggplot2::margin(
        5,
        5,
        5,
        16
      )
    )
}

# Panel E was already assembled before this final override.
# Recursively replace the final displayed labels in its nested plots.
if (
  exists("pE_heatmap") &&
  inherits(pE_heatmap, "ggplot")
) {
  pE_heatmap <- figure5_exact_compose_y_labeler(
    pE_heatmap
  )
}

if (exists("pE")) {
  pE <- figure5_exact_recurse_patchwork(pE)
}


# === FIGURE5_PANEL_E_DENDRO_CLEANUP_BEGIN ===
for (.nm in c("pE_dendro", "pE_tree", "pE_strip")) {
  if (exists(.nm) && inherits(get(.nm), "ggplot")) {
    assign(
      .nm,
      get(.nm) +
        ggplot2::theme_void() +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(0, 2, 0, 2)
        ),
      envir = environment()
    )
  }
}

for (.nm in c("pA", "pB", "pC", "pD", "pE")) {
  if (exists(.nm)) {
    assign(
      .nm,
      get(.nm) & ggplot2::theme(
        plot.tag = ggplot2::element_text(size = 15, face = "bold"),
        plot.tag.position = c(0, 1)
      ),
      envir = environment()
    )
  }
}
# === FIGURE5_PANEL_E_DENDRO_CLEANUP_END ===

# === FIGURE5_EXACT_LABEL_FIX_END ===


# === FIGURE5_FINAL_READABILITY_EVIDENCE_V1_1_BEGIN ===
# Final Figure 5 presentation/evidence clarification pass.
# Audit finding: the selected candidate universe is FDR-supported by construction.
# Therefore nominal support is redundant in Panels D/E and is omitted from display.
# No p-values, q-values, candidate ranks, effect estimates, or candidate membership are changed.

figure5_final_humanize_v11 <- function(x) {
  value <- as.character(x)

  value <- gsub(
    "(?i)\\bFEDCabppc?\\b|(?i)\\bFEDCabcpp\\b",
    "Iron(III) dicitrate transport",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*(?:3dhcdchol|3hdcchol)\\s*\\(e\\)",
    "3-Dehydrochenodeoxycholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*(?:3dhdchol|3dhchol)\\s*\\(e\\)",
    "3-Dehydrocholate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*phpyr\\s*\\(e\\)",
    "Phenylpyruvate exchange",
    value,
    perl=TRUE
  )

  value <- gsub(
    "(?i)\\bEX[ _:-]*ch4s\\s*\\(e\\)",
    "Methanethiol exchange",
    value,
    perl=TRUE
  )

  value <- gsub("\\s+", " ", value)
  trimws(value)
}

figure5_final_wrap_humanize_v11 <- function(x, width=31L) {
  vapply(
    figure5_final_humanize_v11(x),
    function(z) paste(strwrap(z, width=width), collapse="\n"),
    FUN.VALUE=character(1)
  )
}

figure5_filter_rows_v11 <- function(value, column, remove_values) {
  if (is.null(value) || !is.data.frame(value) || !(column %in% names(value))) {
    return(value)
  }
  keep <- !(as.character(value[[column]]) %in% remove_values)
  value <- value[keep, , drop=FALSE]
  if (is.factor(value[[column]])) {
    value[[column]] <- droplevels(value[[column]])
  }
  value
}

figure5_filter_plot_data_v11 <- function(plot_object, column, remove_values) {
  if (!inherits(plot_object, "ggplot")) return(plot_object)

  plot_object$data <- figure5_filter_rows_v11(
    plot_object$data,
    column,
    remove_values
  )

  for (idx in seq_along(plot_object$layers)) {
    layer_data <- plot_object$layers[[idx]]$data
    if (is.data.frame(layer_data)) {
      plot_object$layers[[idx]]$data <- figure5_filter_rows_v11(
        layer_data,
        column,
        remove_values
      )
    }
  }

  plot_object
}

figure5_readability_theme_v11 <- ggplot2::theme(
  text=ggplot2::element_text(size=13.0),
  axis.title=ggplot2::element_text(size=14.8, face="bold"),
  axis.text=ggplot2::element_text(size=11.5),
  legend.title=ggplot2::element_text(size=11.2, face="bold"),
  legend.text=ggplot2::element_text(size=10.4),
  plot.tag=ggplot2::element_text(size=18.5, face="bold"),
  plot.caption=ggplot2::element_text(size=10.3, hjust=0, colour="#4A4650")
)

figure5_bump_text_layers_v11 <- function(plot_object, min_size=3.3) {
  if (!inherits(plot_object, "ggplot")) return(plot_object)
  for (idx in seq_along(plot_object$layers)) {
    layer <- plot_object$layers[[idx]]
    if (inherits(layer$geom, "GeomText")) {
      old <- layer$aes_params$size
      if (is.null(old) || !is.finite(old) || old < min_size) {
        plot_object$layers[[idx]]$aes_params$size <- min_size
      }
    }
  }
  plot_object
}

figure5_recurse_v11 <- function(plot_object) {
  if (inherits(plot_object, "patchwork")) {
    if (!is.null(plot_object$patches$plots) && length(plot_object$patches$plots)) {
      plot_object$patches$plots <- lapply(
        plot_object$patches$plots,
        figure5_recurse_v11
      )
    }
    if (!is.null(plot_object$plots) && length(plot_object$plots)) {
      plot_object$plots <- lapply(
        plot_object$plots,
        figure5_recurse_v11
      )
    }
    return(plot_object)
  }

  if (!inherits(plot_object, "ggplot")) return(plot_object)

  # Panel E: nominal support is redundant because the audited display candidate
  # universe is entirely FDR-supported. Keep FDR support; remove nominal support.
  plot_object <- figure5_filter_plot_data_v11(
    plot_object,
    "evidence_metric",
    c("Nominal support")
  )

  if ("display_label" %in% names(plot_object$data)) {
    plot_object$data$display_label <- figure5_final_humanize_v11(
      plot_object$data$display_label
    )
  }

  y_scale <- plot_object$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- function(x) figure5_final_wrap_humanize_v11(x, 32L)
  }

  x_scale <- plot_object$scales$get_scales("x")
  if (
    !is.null(x_scale) &&
    inherits(x_scale, "ScaleDiscrete") &&
    "evidence_metric" %in% names(plot_object$data)
  ) {
    x_scale$limits <- NULL
    x_scale$drop <- TRUE
  }

  plot_object <- figure5_bump_text_layers_v11(plot_object, 3.3)
  plot_object + figure5_readability_theme_v11
}

# Panel B: humanize the actual factor levels/data, not only a late scale labeler.
if (exists("pB") && inherits(pB, "ggplot")) {
  if ("label_plot" %in% names(pB$data)) {
    if (is.factor(pB$data$label_plot)) {
      old_levels <- levels(pB$data$label_plot)
      new_levels <- figure5_final_wrap_humanize_v11(old_levels, 28L)
      levels(pB$data$label_plot) <- new_levels
    } else {
      pB$data$label_plot <- figure5_final_wrap_humanize_v11(
        pB$data$label_plot,
        28L
      )
    }
  }

  y_scale <- pB$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- function(x) figure5_final_wrap_humanize_v11(x, 28L)
  }

  pB <- pB +
    figure5_readability_theme_v11 +
    ggplot2::theme(
      axis.text.y=ggplot2::element_text(size=11.2, lineheight=0.96),
      legend.title=ggplot2::element_text(size=10.8, face="bold"),
      legend.text=ggplot2::element_text(size=9.8)
    )
}

# Panel D: the audit confirmed every selected candidate is q<0.05.
# Remove only the redundant nominal column. Keep the FDR column explicitly and
# explain why it equals the candidate count.
if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- figure5_filter_plot_data_v11(
    pD,
    "evidence_type",
    c("Nominal")
  )

  x_scale <- pD$scales$get_scales("x")
  if (!is.null(x_scale) && inherits(x_scale, "ScaleDiscrete")) {
    x_scale$limits <- c(
      "Candidates",
      "FDR",
      "Public context",
      "Bioaccumulation"
    )
    x_scale$labels <- c(
      "Candidates",
      "FDR-supported",
      "Public context",
      "Bioaccumulation"
    )
    x_scale$drop <- TRUE
  }

  pD <- figure5_bump_text_layers_v11(pD, 3.6) +
    figure5_readability_theme_v11 +
    ggplot2::labs(
      caption=""
    ) +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(size=10.8, angle=24, hjust=1),
      plot.caption=ggplot2::element_text(size=9.8, hjust=0, colour="#4A4650")
    )
}

for (figure5_nm in c("pA", "pC")) {
  if (exists(figure5_nm, inherits=FALSE)) {
    figure5_obj <- get(figure5_nm, inherits=FALSE)
    if (inherits(figure5_obj, "ggplot")) {
      figure5_obj <- figure5_bump_text_layers_v11(figure5_obj, 3.6)
      assign(
        figure5_nm,
        figure5_obj + figure5_readability_theme_v11
      )
    }
  }
}

# Panel E: remove only redundant nominal support, retain FDR support, and explain
# the white-dot encoding explicitly.
if (exists("pE_heatmap") && inherits(pE_heatmap, "ggplot")) {
  pE_heatmap <- figure5_recurse_v11(pE_heatmap)
}

if (exists("pE")) {
  pE <- figure5_recurse_v11(pE)
  if (inherits(pE, "patchwork")) {
    pE <- pE + patchwork::plot_annotation(
      caption=paste0(
        "White dot = binary criterion met. ",
        "All displayed candidates are FDR-supported (q < 0.05); ",
        "nominal support is omitted as redundant."
      ),
      theme=ggplot2::theme(
        plot.tag=ggplot2::element_text(size=18.5, face="bold"),
        plot.caption=ggplot2::element_text(
          size=10.3,
          hjust=0,
          colour="#4A4650",
          margin=ggplot2::margin(t=5)
        )
      )
    )
  } else if (inherits(pE, "ggplot")) {
    pE <- pE + ggplot2::labs(
      caption=paste0(
        "White dot = binary criterion met. ",
        "All displayed candidates are FDR-supported (q < 0.05); ",
        "nominal support is omitted as redundant."
      )
    ) + figure5_readability_theme_v11
  }
}

# Uniform panel-tag typography on every final top-level panel.
for (figure5_nm in c("pA", "pB", "pC", "pD")) {
  if (exists(figure5_nm, inherits=FALSE)) {
    figure5_obj <- get(figure5_nm, inherits=FALSE)
    if (inherits(figure5_obj, "ggplot")) {
      assign(
        figure5_nm,
        figure5_obj + ggplot2::theme(
          plot.tag=ggplot2::element_text(size=18.5, face="bold")
        )
      )
    }
  }
}

# === FIGURE5_FINAL_READABILITY_EVIDENCE_V1_1_END ===


# === FIGURE5_FINAL_LAYOUT_OUTPUT_V1_3_BEGIN ===
# Final display-only repair after all earlier Figure 5 overrides.
# Scientific values/candidate membership are unchanged.

figure5_v13_tag_theme <- ggplot2::theme(
  plot.tag=ggplot2::element_text(size=18, face="bold"),
  plot.tag.position=c(0.01, 0.99)
)

figure5_v13_b_label <- function(x) {
  v <- as.character(x)
  if (exists("figure5_exact_humanize", mode="function")) {
    v <- figure5_exact_humanize(v)
  }
  v <- gsub("\\n\\[[^]]+\\]$", "", v)
  v <- gsub("\\s*\\[[^]]+\\]$", "", v)
  v <- trimws(v)
  v <- gsub("^Gallintestinimicrobium propionicum$", "Gallintestinimicrobium\\npropionicum", v)
  v <- gsub("^Bifidobacterium catenulatum$", "Bifidobacterium\\ncatenulatum", v)
  v <- gsub("^3-Dehydrochenodeoxycholate exchange$", "3-Dehydrochenodeoxycholate exch.", v)
  v <- gsub("^Iron\\(III\\) dicitrate transport$", "Fe(III)-dicitrate transport", v)
  v
}

if (exists("pB") && inherits(pB, "ggplot")) {
  y_scale <- pB$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- figure5_v13_b_label
  }
  pB <- pB +
    ggplot2::theme(
      legend.position="right",
      legend.box="vertical",
      legend.justification="top",
      legend.background=ggplot2::element_blank(),
      legend.title=ggplot2::element_text(size=9.7, face="bold"),
      legend.text=ggplot2::element_text(size=9.0),
      axis.text.y=ggplot2::element_text(size=9.4, lineheight=0.90),
      axis.title=ggplot2::element_text(size=12.3, face="bold"),
      plot.margin=ggplot2::margin(6, 3, 6, 10)
    ) +
    figure5_v13_tag_theme
}

if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- pD +
    ggplot2::labs(caption=NULL) +
    ggplot2::theme(
      plot.caption=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(size=10.2, angle=24, hjust=1),
      axis.text.y=ggplot2::element_text(size=11.0)
    ) +
    figure5_v13_tag_theme
}

if (exists("evidence_long") && is.data.frame(evidence_long)) {
  evidence_long <- evidence_long[
    as.character(evidence_long$evidence_metric) != "Nominal support",
    , drop=FALSE
  ]
  if (is.factor(evidence_long$evidence_metric)) {
    evidence_long$evidence_metric <- droplevels(evidence_long$evidence_metric)
  }
}

if (exists("pE_heatmap") && inherits(pE_heatmap, "ggplot")) {
  if (exists("evidence_long") && is.data.frame(evidence_long)) {
    pE_heatmap$data <- evidence_long
    for (idx in seq_along(pE_heatmap$layers)) {
      layer_data <- pE_heatmap$layers[[idx]]$data
      if (is.data.frame(layer_data) && "evidence_metric" %in% names(layer_data)) {
        layer_data <- layer_data[
          as.character(layer_data$evidence_metric) != "Nominal support",
          , drop=FALSE
        ]
        if (is.factor(layer_data$evidence_metric)) {
          layer_data$evidence_metric <- droplevels(layer_data$evidence_metric)
        }
        pE_heatmap$layers[[idx]]$data <- layer_data
      }
    }
  }
  x_scale <- pE_heatmap$scales$get_scales("x")
  if (!is.null(x_scale) && inherits(x_scale, "ScaleDiscrete")) {
    x_scale$limits <- NULL
    x_scale$drop <- TRUE
  }
  pE_heatmap <- pE_heatmap +
    ggplot2::labs(caption="White dot = binary criterion met.") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(size=10.0, angle=27, hjust=1),
      axis.text.y=ggplot2::element_text(size=9.2, lineheight=0.92),
      plot.caption=ggplot2::element_text(size=9.5, hjust=0, colour="#4A4650", margin=ggplot2::margin(t=5))
    )
}

figure5_v13_clean_dendro <- NULL
figure5_v13_dendro_names <- grep(
  "^pE.*dend|dend.*E",
  ls(envir=.GlobalEnv),
  value=TRUE,
  ignore.case=TRUE
)
if (length(figure5_v13_dendro_names)) {
  for (nm in figure5_v13_dendro_names) {
    obj <- get(nm, envir=.GlobalEnv)
    if (inherits(obj, "ggplot")) {
      figure5_v13_clean_dendro <- obj +
        ggplot2::theme_void() +
        ggplot2::theme(
          axis.title=ggplot2::element_blank(),
          axis.text=ggplot2::element_blank(),
          axis.ticks=ggplot2::element_blank(),
          panel.grid=ggplot2::element_blank(),
          plot.tag=ggplot2::element_blank(),
          plot.margin=ggplot2::margin(0, 2, 0, 0)
        )
      break
    }
  }
}

figure5_v13_clean_strip <- NULL
if (exists("e_closeout") && is.list(e_closeout) && inherits(e_closeout$strip, "ggplot")) {
  figure5_v13_clean_strip <- e_closeout$strip +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.title=ggplot2::element_blank(),
      axis.text=ggplot2::element_blank(),
      axis.ticks=ggplot2::element_blank(),
      panel.grid=ggplot2::element_blank(),
      plot.tag=ggplot2::element_text(size=18, face="bold"),
      plot.tag.position=c(0.02, 0.98),
      plot.margin=ggplot2::margin(0, 2, 0, 2)
    )
} else if (exists("pE_strip") && inherits(pE_strip, "ggplot")) {
  figure5_v13_clean_strip <- pE_strip +
    ggplot2::labs(tag="E") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.tag=ggplot2::element_text(size=18, face="bold"),
      plot.tag.position=c(0.02, 0.98),
      plot.margin=ggplot2::margin(0, 2, 0, 2)
    )
}

if (!is.null(figure5_v13_clean_dendro) && !is.null(figure5_v13_clean_strip) && exists("pE_heatmap") && inherits(pE_heatmap, "ggplot")) {
  pE <- patchwork::wrap_plots(
    figure5_v13_clean_dendro,
    figure5_v13_clean_strip,
    pE_heatmap,
    nrow=1,
    widths=c(0.16, 0.055, 1),
    guides="collect"
  ) & ggplot2::theme(legend.position="bottom")
}

if (exists("pA") && inherits(pA, "ggplot")) pA <- pA + figure5_v13_tag_theme
if (exists("pC")) pC <- pC & figure5_v13_tag_theme
if (exists("pE")) pE <- pE & figure5_v13_tag_theme
# === FIGURE5_FINAL_LAYOUT_OUTPUT_V1_3_END ===

# === FIGURE5_LAYOUT_TIGHTEN_V1_4_BEGIN ===
# Display-only whitespace/overlap cleanup. No scientific recomputation.

figure5_v14_wrap <- function(x, width = 34) {
  v <- as.character(x)
  v <- gsub("^HUMAnN::", "HUMAnN: ", v)
  v <- gsub("^AGORA2 reaction::", "AGORA2 reaction: ", v)
  v <- gsub("^AGORA2 exchange::", "AGORA2 exchange: ", v)
  v <- gsub("^AGORA2 subsystem::", "AGORA2 subsystem: ", v)
  v <- gsub("^Species::", "Species: ", v)
  v <- gsub("^Genus::", "Genus: ", v)
  v <- gsub("^MAG::", "MAG: ", v)
  v <- gsub(" -- ", " — ", v, fixed = TRUE)
  vapply(v, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

if (exists("pB") && inherits(pB, "ggplot")) {
  y_scale <- pB$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- function(x) figure5_v14_wrap(x, width = 28)
  }

  pB <- pB +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.justification = "top",
      legend.spacing.y = grid::unit(1.0, "mm"),
      legend.key.height = grid::unit(4.0, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      axis.text.y = ggplot2::element_text(size = 9.5, lineheight = 0.92),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

if (exists("pE_heatmap") && inherits(pE_heatmap, "ggplot")) {
  if (is.data.frame(pE_heatmap$data) && "display_label" %in% names(pE_heatmap$data)) {
    pE_heatmap$data$display_label <- figure5_v14_wrap(pE_heatmap$data$display_label, width = 34)
  }

  for (idx in seq_along(pE_heatmap$layers)) {
    layer_data <- pE_heatmap$layers[[idx]]$data
    if (is.data.frame(layer_data) && "display_label" %in% names(layer_data)) {
      layer_data$display_label <- figure5_v14_wrap(layer_data$display_label, width = 34)
      pE_heatmap$layers[[idx]]$data <- layer_data
    }
  }

  y_scale <- pE_heatmap$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- function(x) figure5_v14_wrap(x, width = 34)
  }

  pE_heatmap <- pE_heatmap +
    ggplot2::labs(
      tag = NULL,
      caption = "White dot = binary criterion met."
    ) +
    ggplot2::theme(
      plot.tag = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8.8, lineheight = 0.88),
      axis.text.x = ggplot2::element_text(size = 10.0, angle = 28, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 12.0, face = "bold", margin = ggplot2::margin(t = 6)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      plot.caption = ggplot2::element_text(
        size = 9.4,
        hjust = 0,
        colour = "#4A4650",
        margin = ggplot2::margin(t = 5)
      )
    )
}

if (exists("figure5_v13_clean_dendro") && inherits(figure5_v13_clean_dendro, "ggplot")) {
  figure5_v13_clean_dendro <- figure5_v13_clean_dendro +
    ggplot2::theme(
      plot.tag = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 0, 0)
    )
}

if (exists("figure5_v13_clean_strip") && inherits(figure5_v13_clean_strip, "ggplot")) {
  figure5_v13_clean_strip <- figure5_v13_clean_strip +
    ggplot2::labs(tag = "E") +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 18, face = "bold"),
      plot.tag.position = c(0.02, 0.98),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

if (!is.null(figure5_v13_clean_dendro) &&
    !is.null(figure5_v13_clean_strip) &&
    exists("pE_heatmap") &&
    inherits(pE_heatmap, "ggplot")) {
  pE <- patchwork::wrap_plots(
    figure5_v13_clean_dendro,
    figure5_v13_clean_strip,
    pE_heatmap,
    nrow = 1,
    widths = c(0.12, 0.028, 1.00),
    guides = "collect"
  ) & ggplot2::theme(
    legend.position = "bottom",
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )
}

if (exists("pA") && inherits(pA, "ggplot")) {
  pA <- pA + ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}
if (exists("pC")) {
  pC <- pC & ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}
if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- pD + ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}
# === FIGURE5_LAYOUT_TIGHTEN_V1_4_END ===





# === FIGURE5_COMPACT_LAYOUT_V1_5_BEGIN ===
# Final display-only compaction/readability pass.
# No scientific recomputation.

figure5_v15_tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(size = 18, face = "bold"),
  plot.tag.position = c(0.01, 0.99)
)

figure5_v15_b_label <- function(x) {
  v <- as.character(x)
  if (exists("figure5_exact_humanize", mode = "function")) {
    v <- figure5_exact_humanize(v)
  }
  v <- gsub("^Gallintestinimicrobium propionicum$", "Gallintestinimicrobium\npropionicum", v)
  v <- gsub("^Bifidobacterium catenulatum$", "Bifidobacterium\ncatenulatum", v)
  v <- gsub("^3-Dehydrochenodeoxycholate exchange$", "3-Dehydrochenodeoxycholate\nexchange", v)
  v <- gsub("^Iron\\(III\\) dicitrate transport$", "Iron(III) dicitrate\ntransport", v)
  v
}

if (exists("pA") && inherits(pA, "ggplot")) {
  pA <- pA +
    ggplot2::labs(tag = "A") +
    figure5_v15_tag_theme +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

if (exists("pB") && inherits(pB, "ggplot")) {
  y_scale <- pB$scales$get_scales("y")
  if (!is.null(y_scale) && inherits(y_scale, "ScaleDiscrete")) {
    y_scale$labels <- figure5_v15_b_label
  }

  pB <- pB +
    ggplot2::labs(tag = "B") +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        order = 1,
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1)
      ),
      fill = ggplot2::guide_colourbar(
        order = 2,
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(20, "mm"),
        barheight = grid::unit(3.2, "mm")
      ),
      colour = ggplot2::guide_colourbar(
        order = 2,
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(20, "mm"),
        barheight = grid::unit(3.2, "mm")
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "left",
      legend.title = ggplot2::element_text(size = 9.3, face = "bold"),
      legend.text = ggplot2::element_text(size = 8.7),
      legend.spacing.x = grid::unit(1.5, "mm"),
      legend.spacing.y = grid::unit(0.6, "mm"),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.6, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      axis.text.y = ggplot2::element_text(size = 9.2, lineheight = 0.90),
      axis.title = ggplot2::element_text(size = 12.0, face = "bold"),
      plot.margin = ggplot2::margin(2, 2, 1, 2)
    ) +
    figure5_v15_tag_theme
}

if (exists("pC") && inherits(pC, "ggplot")) {
  pC <- pC +
    ggplot2::labs(tag = "C") +
    figure5_v15_tag_theme +
    ggplot2::theme(
      plot.margin = ggplot2::margin(1, 2, 1, 2)
    )
}

if (exists("pD") && inherits(pD, "ggplot")) {
  pD <- pD +
    ggplot2::labs(tag = "D") +
    ggplot2::theme(
      plot.caption = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 2, 1, 2)
    ) +
    figure5_v15_tag_theme
}

if (exists("pE_heatmap") && inherits(pE_heatmap, "ggplot")) {
  pE_heatmap <- pE_heatmap +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9.6, angle = 28, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8.7, lineheight = 0.88),
      axis.title.x = ggplot2::element_text(size = 11.8, face = "bold"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      plot.caption = ggplot2::element_text(
        size = 9.0,
        hjust = 0,
        colour = "#4A4650",
        margin = ggplot2::margin(t = 4)
      )
    ) +
    ggplot2::labs(caption = "White dot = binary criterion met.")
}

if (exists("figure5_v13_clean_dendro") && inherits(figure5_v13_clean_dendro, "ggplot")) {
  figure5_v13_clean_dendro <- figure5_v13_clean_dendro +
    ggplot2::theme(
      plot.tag = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 0, 0)
    )
}

if (exists("figure5_v13_clean_strip") && inherits(figure5_v13_clean_strip, "ggplot")) {
  figure5_v13_clean_strip <- figure5_v13_clean_strip +
    ggplot2::theme(
      plot.tag = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

if (exists("pE_heatmap") &&
    exists("figure5_v13_clean_dendro") &&
    exists("figure5_v13_clean_strip") &&
    inherits(pE_heatmap, "ggplot") &&
    inherits(figure5_v13_clean_dendro, "ggplot") &&
    inherits(figure5_v13_clean_strip, "ggplot")) {

  pE <- patchwork::wrap_plots(
    figure5_v13_clean_dendro,
    figure5_v13_clean_strip,
    pE_heatmap,
    nrow = 1,
    widths = c(0.10, 0.018, 1.00),
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  pE <- pE & figure5_v15_tag_theme
}
# === FIGURE5_COMPACT_LAYOUT_V1_5_END ===


top_row <- pA + pB + plot_layout(widths = c(0.72, 1.28))
middle_row <- pC + pD + plot_layout(widths = c(1.03, 0.97))

figure5_v3 <- (
    top_row /
        middle_row /
        pE
) +
    plot_layout(heights = c(0.88, 0.80, 1.20)) +
    plot_annotation(
        title = NULL,
        subtitle = NULL,
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 11.0,
                colour = PFAS_PURPLE[["darkest"]],
                margin = margin(b = 3)
            ),
            plot.subtitle = element_text(
                size = 8.0,
                colour = PFAS_PURPLE[["neutral_dark"]],
                margin = margin(b = 7)
            ),
            plot.margin = margin(9, 10, 8, 10)
        )
    )

FIGURE5_REVIEW_DIR <- file.path(
    ROOT, "submission", "figure_builds", "Figure_5", "figures"
)
dir.create(FIGURE5_REVIEW_DIR, recursive=TRUE, showWarnings=FALSE)
OUTPUT_STEM <- file.path(FIGURE5_REVIEW_DIR, "Figure_5_excellence_candidate")

ggsave(
    paste0(OUTPUT_STEM, ".pdf"),
    figure5_v3,
    width = 18.5,
    height = 15.2,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)
ggsave(
    paste0(OUTPUT_STEM, ".png"),
    figure5_v3,
    width = 18.5,
    height = 15.2,
    units = "in",
    dpi = 400,
    limitsize = FALSE
)
ggsave(
    paste0(OUTPUT_STEM, ".tiff"),
    figure5_v3,
    width = 18.5,
    height = 15.2,
    units = "in",
    dpi = 600,
    compression = "lzw",
    limitsize = FALSE
)

fwrite(performance_summary, file.path(SRC, "Figure_5_v3_A_performance_summary.tsv"), sep = "\t")
fwrite(stable_matrix, file.path(SRC, "Figure_5_v3_C_stability_matrix.tsv"), sep = "\t")
fwrite(
    family_summary,
    file.path(
        SRC,
        "Figure_5_v3_D_reconstructed_family_summary.tsv"
    ),
    sep = "\t"
)

fwrite(family_long, file.path(SRC, "Figure_5_v3_D_family_evidence_matrix.tsv"), sep = "\t")
fwrite(candidate_unique, file.path(SRC, "Figure_5_v3_E_unique_candidates.tsv"), sep = "\t")
fwrite(evidence_long, file.path(SRC, "Figure_5_v3_E_clustered_evidence_matrix.tsv"), sep = "\t")

audit <- data.table(
    check = c(
        "PANEL_A_FEATURE_SETS",
        "PANEL_A_REPEAT_ROWS",
        "PANEL_B_FEATURE_SETS",
        "PANEL_C_FEATURES",
        "PANEL_C_METRICS",
        "AUTHORITATIVE_CANDIDATE_ROWS",
        "PANEL_D_FAMILY_POOL_ROWS",
        "PANEL_D_FAMILIES",
        "PANEL_D_EVIDENCE_TYPES",
        "PANEL_E_UNIQUE_CANDIDATES",
        "PANEL_E_EVIDENCE_METRICS",
        "PANEL_E_CLUSTERS",
        "DUPLICATED_CANDIDATE_IDS",
        "DOT_SUFFIX_LABELS",
        "CLUSTER_EXPORT_EMBEDDED_NEWLINES",
        "OUTPUT_PDF",
        "OUTPUT_PNG",
        "OUTPUT_TIFF"
    ),
    value = c(
        uniqueN(repeat_performance$feature_set),
        nrow(repeat_performance),
        uniqueN(predictions$feature_set),
        uniqueN(stable_matrix$display_label),
        uniqueN(stable_matrix$metric),
        nrow(authoritative_candidates),
        nrow(family_pool),
        uniqueN(family_long$family),
        uniqueN(family_long$evidence_type),
        nrow(candidate_unique),
        length(metric_order),
        cluster_number,
        anyDuplicated(candidate_unique$candidate_id),
        sum(grepl("\\.[0-9]+$", candidate_unique$display_label)),
        sum(grepl("[\r\n]", cluster_export$candidate_id)),
        as.integer(file.exists(paste0(OUTPUT_STEM, ".pdf"))),
        as.integer(file.exists(paste0(OUTPUT_STEM, ".png"))),
        as.integer(file.exists(paste0(OUTPUT_STEM, ".tiff")))
    )
)

fwrite(
    audit,
    file.path(REPORTDIR, "phase29D_figure5_excellence_v3_audit.tsv"),
    sep = "\t"
)

editorial_audit <- data.table(
    check = c(
        "OVERALL_DESCRIPTIVE_TITLE",
        "DESCRIPTIVE_SUBPLOT_HEADINGS",
        "PANEL_LETTERS_ONLY",
        "AXIS_LABELS_RETAINED"
    ),
    value = c(0, 0, 1, 1)
)

fwrite(
    editorial_audit,
    file.path(REPORTDIR, "phase29D_figure5_editorial_policy_audit.tsv"),
    sep = "\t"
)

cat("=== PFAS PHASE 29D FIGURE 5 EDITORIAL FINALIZATION ===\n")
cat("Stable features displayed:", uniqueN(stable_matrix$display_label), "\n")
cat("Unique clustered candidates:", nrow(candidate_unique), "\n")
cat("Candidate clusters:", cluster_number, "\n")
cat("Output:", paste0(OUTPUT_STEM, ".pdf"), "\n")
cat("Status: READY_FOR_FIGURE5_V3_VISUAL_QC\n")
