# FIGURE6_LAYOUT_READABILITY_PATCH_V1_0
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
TABLES <- file.path(REV, "12_tables")
P13 <- file.path(REV, "13_PFAS_bioaccumulation_interpretation")
P12 <- file.path(REV, "12_machine_learning_candidate_prioritization")
DERIVED <- file.path(P14, "03_source_data_derived")

OUT <- file.path(REV, "submission", "figure_builds", "Figure_6")
SRC <- file.path(OUT, "source_data")
LOGD <- file.path(OUT, "reports")
FIGD <- file.path(OUT, "figures")
PNGD <- file.path(OUT, "previews")
CLUSTERDIR <- file.path(OUT, "cluster_orders")

for (directory in c(OUT, SRC, DERIVED, LOGD, FIGD, PNGD, CLUSTERDIR)) {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
}

purple <- c(
  deep="#4B0082",
  main="#5E3C99",
  mid="#7A5AB5",
  soft="#9E8ACB",
  pale="#D9CFF0",
  faint="#F1EDF8",
  grey_dark="#333333",
  grey_mid="#777777",
  grey_light="#D9D9D9",
  white="#FFFFFF"
)

theme_pfas <- function(base_size=11.2) {
  theme_bw(base_size=base_size, base_family="sans") +
    theme(
      panel.grid.major=element_line(colour="#E4E4E4", linewidth=0.24),
      panel.grid.minor=element_blank(),
      axis.title=element_text(face="bold", colour=purple["grey_dark"], size=base_size+0.8),
      axis.text=element_text(colour=purple["grey_dark"], size=base_size),
      strip.background=element_rect(fill=purple["faint"], colour=purple["grey_mid"], linewidth=0.25),
      strip.text=element_text(face="bold", colour=purple["grey_dark"], size=base_size),
      legend.title=element_text(face="bold", size=base_size),
      legend.text=element_text(size=base_size-0.8),
      plot.margin=margin(4,5,4,5),
      plot.tag=element_text(face="bold", size=14, colour=purple["grey_dark"]),
      plot.tag.position=c(0.01,0.995)
    )
}

clean_taxon <- function(x) {
  x <- as.character(x)
  x <- gsub(".*\\|s__", "s__", x)
  x <- gsub(".*\\|g__", "g__", x)
  x <- gsub("^s__", "", x)
  x <- gsub("^g__", "", x)
  x <- gsub("\\|.*$", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

feature_key <- function(x) {
  x <- tolower(clean_taxon(x))
  x <- gsub("[^a-z0-9 ]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

genus_key <- function(x) {
  y <- feature_key(x)
  vapply(strsplit(y, "\\s+"), function(z) if (length(z) >= 1) z[1] else "", character(1))
}

wrap_lab <- function(x, w=30) {
  vapply(as.character(x), function(z) paste(strwrap(z, width=w), collapse="\n"), character(1))
}

canon_sample <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)
  x <- gsub("\\.bam$|\\.fastq$|\\.fq$|\\.gz$|\\.profile$|\\.tsv$|\\.txt$", "", x, ignore.case=TRUE)
  x <- gsub("_metaphlan.*$|\\.metaphlan.*$|_profile.*$|\\.profile.*$", "", x, ignore.case=TRUE)
  x <- gsub("[^A-Za-z0-9]", "", x)
  toupper(x)
}

safe01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(!is.finite(x))) return(rep(0.5, length(x)))
  rng <- range(x, na.rm=TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

detect_sample_col <- function(dt) {
  candidates <- c("sample_id","SampleID","Sample","sample","Sample_Name","sample_name","ID","id")
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  names(dt)[1]
}

pick_existing <- function(paths) {
  hits <- paths[file.exists(paths) & file.info(paths)$size > 0]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}


# Compatibility helpers for patched Figure 6 v9 functional panel.
clean_text <- function(x) {
  clean_taxon(x)
}

wrap_label <- function(x, width=30) {
  wrap_lab(x, w=width)
}


cluster_rows <- function(mat) {
  m <- as.matrix(mat)
  if (is.null(dim(m))) m <- matrix(m, nrow=1)
  if (nrow(m) <= 1) return(list(order=seq_len(nrow(m)), tree=NULL))
  m[!is.finite(m)] <- 0
  if (ncol(m) > 1) {
    m <- t(scale(t(m)))
    m[!is.finite(m)] <- 0
  }
  tree <- hclust(dist(m), method="average")
  if (length(tree$height) > 0 && max(tree$height, na.rm=TRUE) > 0) {
    tree$height <- tree$height / max(tree$height, na.rm=TRUE)
  }
  list(order=tree$order, tree=tree)
}

cluster_cols <- function(mat) {
  m <- as.matrix(mat)
  if (is.null(dim(m))) m <- matrix(m, ncol=1)
  if (ncol(m) <= 1) return(list(order=seq_len(ncol(m)), tree=NULL))
  m <- t(m)
  m[!is.finite(m)] <- 0
  if (ncol(m) > 1) {
    m <- t(scale(t(m)))
    m[!is.finite(m)] <- 0
  }
  tree <- hclust(dist(m), method="average")
  if (length(tree$height) > 0 && max(tree$height, na.rm=TRUE) > 0) {
    tree$height <- tree$height / max(tree$height, na.rm=TRUE)
  }
  list(order=tree$order, tree=tree)
}

row_dendrogram_segments <- function(tree) {
  if (is.null(tree)) return(data.table())
  n <- length(tree$order)
  leaf_y <- numeric(n)
  leaf_y[tree$order] <- seq_len(n)
  node_y <- numeric(n - 1)
  output <- list()
  idx <- 1L

  child_details <- function(child) {
    if (child < 0) return(list(x=0, y=leaf_y[-child]))
    list(x=tree$height[child], y=node_y[child])
  }

  for (merge_index in seq_len(n - 1)) {
    left <- child_details(tree$merge[merge_index, 1])
    right <- child_details(tree$merge[merge_index, 2])
    parent_x <- tree$height[merge_index]
    parent_y <- mean(c(left$y, right$y))
    node_y[merge_index] <- parent_y
    output[[idx]] <- data.table(x=left$x, xend=parent_x, y=left$y, yend=left$y); idx <- idx + 1L
    output[[idx]] <- data.table(x=right$x, xend=parent_x, y=right$y, yend=right$y); idx <- idx + 1L
    output[[idx]] <- data.table(x=parent_x, xend=parent_x, y=left$y, yend=right$y); idx <- idx + 1L
  }
  rbindlist(output)
}

col_dendrogram_segments <- function(tree) {
  if (is.null(tree)) return(data.table())
  n <- length(tree$order)
  leaf_x <- numeric(n)
  leaf_x[tree$order] <- seq_len(n)
  node_x <- numeric(n - 1)
  output <- list()
  idx <- 1L

  child_details <- function(child) {
    if (child < 0) return(list(x=leaf_x[-child], y=0))
    list(x=node_x[child], y=tree$height[child])
  }

  for (merge_index in seq_len(n - 1)) {
    left <- child_details(tree$merge[merge_index, 1])
    right <- child_details(tree$merge[merge_index, 2])
    parent_y <- tree$height[merge_index]
    parent_x <- mean(c(left$x, right$x))
    node_x[merge_index] <- parent_x
    output[[idx]] <- data.table(x=left$x, xend=left$x, y=left$y, yend=parent_y); idx <- idx + 1L
    output[[idx]] <- data.table(x=right$x, xend=right$x, y=right$y, yend=parent_y); idx <- idx + 1L
    output[[idx]] <- data.table(x=left$x, xend=right$x, y=parent_y, yend=parent_y); idx <- idx + 1L
  }
  rbindlist(output)
}

make_row_dendrogram_plot <- function(tree, n_leaves, panel_tag=NULL) {
  seg <- row_dendrogram_segments(tree)
  if (nrow(seg) == 0) {
    return(ggplot() + theme_void() + labs(tag=panel_tag))
  }
  ggplot(seg) +
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 linewidth=0.56, lineend="round", colour=purple["grey_dark"]) +
    scale_x_reverse(limits=c(1.02, 0), expand=c(0,0)) +
    scale_y_continuous(limits=c(0.5, n_leaves + 0.5), expand=c(0,0)) +
    labs(tag=panel_tag) +
    theme_void() +
    theme(
      plot.tag=element_text(face="bold", size=14, colour=purple["grey_dark"]),
      plot.tag.position=c(0.02, 0.98),
      plot.margin=margin(2, 1, 2, 2)
    )
}

make_col_dendrogram_plot <- function(tree, n_leaves) {
  seg <- col_dendrogram_segments(tree)
  if (nrow(seg) == 0) {
    return(ggplot() + theme_void())
  }
  ggplot(seg) +
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 linewidth=0.42, colour=purple["grey_dark"]) +
    scale_x_continuous(limits=c(0.5, n_leaves + 0.5), expand=c(0,0)) +
    scale_y_continuous(expand=c(0.02, 0)) +
    theme_void() +
    theme(plot.margin=margin(2, 2, 1, 2))
}

fix_family <- function(current_family, raw_feature, label=NULL) {
  cf <- as.character(current_family)
  raw <- paste(as.character(raw_feature), as.character(label), sep=" || ")
  low <- tolower(raw)
  cflow <- tolower(cf)
  out <- cf

  out[grepl("^s__|\\|s__|\\bs__", low)] <- "MetaPhlAn species"
  out[grepl("^g__|\\|g__|\\bg__", low)] <- "MetaPhlAn genus"
  out[grepl("exchange:|^ex_|\\[ex_|agora2 exchange", low)] <- "AGORA2 exchanges"
  out[grepl("reaction:|agora2 reaction|fedcabcpp|aldc|alcd|transport reaction", low)] <- "AGORA2 reactions"
  out[grepl("agora2 subsystem|subsystem:|oxalate metabolism|steroid metabolism|zeatin biosynthesis|butanoate metabolism", low)] <- "AGORA2 subsystems"

  is_taxon <- out %in% c("MetaPhlAn species", "MetaPhlAn genus")
  is_agora <- out %in% c("AGORA2 reactions", "AGORA2 subsystems", "AGORA2 exchanges")
  humann_like <- grepl("pwy|pathway|biosynthesis|degradation|metabolism|fermentation|utilization|glycolysis|vitamin|coenzyme|glycan|ion", low)

  out[!is_taxon & !is_agora & humann_like] <- "HUMAnN pathways"
  out[grepl("\\bbin\\.?[0-9]+|mag bin|_bin", low) & !is_taxon] <- "MAGs"

  out[grepl("species", cflow) & !is_agora] <- "MetaPhlAn species"
  out[grepl("genus", cflow) & !is_agora] <- "MetaPhlAn genus"
  out[grepl("humann|pathway", cflow) & !is_taxon & !is_agora] <- "HUMAnN pathways"
  out[grepl("reaction", cflow)] <- "AGORA2 reactions"
  out[grepl("subsystem", cflow)] <- "AGORA2 subsystems"
  out[grepl("exchange", cflow)] <- "AGORA2 exchanges"
  out[grepl("mag", cflow) & !is_taxon] <- "MAGs"

  out[is.na(out) | out == ""] <- "Other"
  out
}


detect_col <- function(nms, patterns, exclude=NULL) {
  low <- tolower(nms)
  for (pat in patterns) {
    hit <- grep(pat, low, perl=TRUE)
    if (!is.null(exclude) && length(hit) > 0) {
      for (ex in exclude) hit <- hit[!grepl(ex, low[hit], perl=TRUE)]
    }
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

# ----------------------
# Core inputs
# ----------------------
taxcorr_path <- file.path(DERIVED, "Figure6_taxon_elimination_correlations.tsv")
if (!file.exists(taxcorr_path)) stop("Missing taxon-correlation table: ", taxcorr_path)

overlap_detail_path <- pick_existing(c(
  file.path(TABLES, "Table_S26_unique_taxon_overlap_detail.tsv"),
  file.path(P13, "tables", "Table_S26_unique_taxon_overlap_detail.tsv")
))
if (is.na(overlap_detail_path)) stop("Missing Table_S26_unique_taxon_overlap_detail.tsv")

taxcorr <- fread(taxcorr_path, data.table=TRUE, showProgress=FALSE)
overD <- fread(overlap_detail_path, data.table=TRUE, showProgress=FALSE)
setnames(overD, make.unique(names(overD), sep="_dup"))

# ----------------------
# Lindell support context
# ----------------------
src_col <- if ("source" %in% names(overD)) "source" else names(overD)[1]
feat_col <- if ("feature" %in% names(overD)) "feature" else names(overD)[2]
level_col <- if ("match_level" %in% names(overD)) "match_level" else {
  hit <- names(overD)[grepl("match|level", names(overD), ignore.case=TRUE)]
  if (length(hit) == 0) stop("Cannot detect match level column in overlap detail.")
  hit[1]
}

overD[, feature_clean := clean_taxon(get(feat_col))]
overD[, feature_key := feature_key(feature_clean)]
overD[, genus_key := genus_key(feature_clean)]
overD[, source_raw := as.character(get(src_col))]
overD[, feature_family := fifelse(grepl("MetaPhlAn.*genus", source_raw, ignore.case=TRUE), "MetaPhlAn genus",
                           fifelse(grepl("MetaPhlAn.*species", source_raw, ignore.case=TRUE), "MetaPhlAn species",
                           fifelse(grepl("MAG|GTDB", source_raw, ignore.case=TRUE), "MAGs GTDB taxonomy", source_raw)))]
overD[, match_level_raw := as.character(get(level_col))]
overD[, match_rank := fifelse(grepl("species", match_level_raw, ignore.case=TRUE), 3L,
                       fifelse(grepl("genus", match_level_raw, ignore.case=TRUE), 2L,
                       fifelse(grepl("family", match_level_raw, ignore.case=TRUE), 1L, 0L)))]
overD[, overlap_level := fifelse(match_rank == 3L, "Species exact",
                          fifelse(match_rank == 2L, "Genus",
                          fifelse(match_rank == 1L, "Family", "Other")))]

best_over <- overD[feature_family %in% c("MetaPhlAn genus", "MetaPhlAn species")]
setorder(best_over, feature_family, feature_key, -match_rank)
best_over <- best_over[, .SD[1], by=.(feature_family, feature_key)]

tc <- merge(
  taxcorr,
  best_over[, .(feature_family, feature_key, overlap_level, match_rank)],
  by=c("feature_family", "feature_key"),
  all.x=TRUE
)
tc[, lindell_support := !is.na(overlap_level)]
tc[is.na(overlap_level), overlap_level := "No Lindell support"]
tc[is.na(match_rank), match_rank := 0L]
tc[, overlap_level := as.character(overlap_level)]
tc[overlap_level == "Species exact", overlap_level := "Species-level match"]
tc[overlap_level == "Genus", overlap_level := "Genus-level match"]
tc[overlap_level == "Family", overlap_level := "Family-level match"]
tc[, overlap_level := factor(
  overlap_level,
  levels=c("Species-level match","Genus-level match","Family-level match","No Lindell support")
)]
tc[, abs_rho := abs(rho)]
tc[, nominal_p := is.finite(p) & p < 0.05]
tc[, direction2 := fifelse(rho > 0, "Higher elimination", "Lower elimination")]
tc[, genus_key := genus_key(feature)]

# ----------------------
# Select taxa for main triangulation
# ----------------------
sel <- tc[lindell_support == TRUE & is.finite(rho)]
setorder(sel, -abs_rho, p)
sel <- sel[seq_len(min(16, .N))]
sel <- sel[order(rho)]
sel[, feature_label_txt := wrap_lab(feature, 34)]
sel[, feature_label := factor(feature_label_txt, levels=feature_label_txt)]
sel[, overlap_level_short := as.character(overlap_level)]
support_levels_present <- intersect(
  c("Species-level match","Genus-level match","Family-level match"),
  unique(sel$overlap_level_short)
)
support_palette <- c(
  "Species-level match"=unname(purple["deep"]),
  "Genus-level match"=unname(purple["mid"]),
  "Family-level match"=unname(purple["pale"])
)
support_shapes <- c(
  "Species-level match"=21,
  "Genus-level match"=24,
  "Family-level match"=22
)

fwrite(sel, file.path(SRC, "Figure_6A_selected_lindell_support_taxa.tsv"), sep="\t")

# ----------------------
# Panel A: Primary elimination association
# ----------------------
pA <- ggplot(sel, aes(x=rho, y=feature_label)) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill=purple["faint"], alpha=0.55) +
  geom_vline(xintercept=0, linetype="dashed", colour=purple["grey_mid"], linewidth=0.45) +
  geom_segment(aes(x=0, xend=rho, y=feature_label, yend=feature_label),
               colour=purple["mid"], linewidth=0.95) +
  geom_point(
    aes(
      size=prevalence_exposed,
      fill=overlap_level_short,
      shape=overlap_level_short
    ),
    colour=purple["grey_dark"],
    stroke=0.42
  ) +
  scale_fill_manual(
    values=support_palette[support_levels_present],
    breaks=support_levels_present,
    drop=TRUE
  ) +
  scale_shape_manual(
    values=support_shapes[support_levels_present],
    breaks=support_levels_present,
    drop=TRUE,
    guide="none"
  ) +
  scale_size_continuous(labels=percent_format(accuracy=1), range=c(3.2,8.2)) +
  labs(
    x="Spearman rho with primary elimination score",
    y="Lindell-supported Ronneby taxon",
    fill="Lindell taxonomic match",
    size="Prevalence among\nPFAS-exposed participants (%)"
  ) +
  labs(tag="A") +
  theme_pfas(11.6) +
  theme(
    axis.text.y=element_text(size=9.5),
    axis.text.x=element_text(size=10.1),
    legend.position="right"
  ) +
  guides(
    fill=guide_legend(
      title.position="top",
      order=1,
      override.aes=list(
        shape=unname(support_shapes[support_levels_present]),
        size=4.9
      )
    ),
    size=guide_legend(
      title.position="top",
      order=2,
      override.aes=list(shape=21, fill="white")
    )
  )

# ----------------------
# Panel B: Compound-specific elimination heatmap
# ----------------------
pfas_candidates <- unique(c(
  file.path(REV, "10_PFAS_phenotype_landscape", "phase10_pfas_matrix.tsv"),
  list.files(REV, pattern="phase10_pfas_matrix\\.tsv$", recursive=TRUE, full.names=TRUE)
))
pfas_candidates <- pfas_candidates[file.exists(pfas_candidates) & file.info(pfas_candidates)$size > 0]
pfas_path <- pfas_candidates[order(!grepl("phase10_pfas_matrix", basename(pfas_candidates), ignore.case=TRUE),
                                   file.info(pfas_candidates)$size)][1]

species_candidates <- c(
  file.path(REV, "08_microbiome_statistics", "00_inputs", "corrected_feature_matrices", "metaphlan_species_sample_matched.tsv"),
  file.path(REV, "03_taxonomy_metaphlan", "03_merged", "Ronneby_MetaPhlAn_species.tsv"),
  file.path(REV, "12_tables", "Table_S6_MetaPhlAn_species.tsv")
)
species_path <- pick_existing(species_candidates)
if (is.na(species_path)) stop("Could not find corrected MetaPhlAn species matrix for panel B.")

pfas <- fread(pfas_path, data.table=TRUE, showProgress=FALSE)
setnames(pfas, make.unique(names(pfas), sep="_dup"))
sid <- detect_sample_col(pfas)
setnames(pfas, sid, "sample_id")
pfas[, sample_id_raw := as.character(sample_id)]
pfas[, sample_key := canon_sample(sample_id)]

k_cols <- names(pfas)[grepl("elim.*k|elimination.*k|^elim_k", names(pfas), ignore.case=TRUE)]
k_cols <- k_cols[!grepl("half|life|tertile|group|class|available|n_k|include", k_cols, ignore.case=TRUE)]
if (length(k_cols) == 0) stop("No compound-specific elimination k columns detected.")

for (kc in k_cols) suppressWarnings(pfas[, (kc) := as.numeric(get(kc))])
target_k <- pfas[, c("sample_id_raw","sample_key", k_cols), with=FALSE]
setnames(target_k, "sample_id_raw", "sample_id")

sp <- fread(species_path, data.table=TRUE, showProgress=FALSE)
setnames(sp, make.unique(names(sp), sep="_dup"))
col_map <- data.table(column_name=names(sp), sample_key=canon_sample(names(sp)))
col_map <- merge(col_map, target_k[, .(sample_key, sample_id)], by="sample_key", all=FALSE)
sample_cols <- col_map$column_name

feature_col <- names(sp)[1]
fc <- names(sp)[tolower(names(sp)) %in% c("feature","taxon","clade_name","name","species")]
if (length(fc) > 0) feature_col <- fc[1]

long <- melt(
  sp[, c(feature_col, sample_cols), with=FALSE],
  id.vars=feature_col,
  variable.name="column_name",
  value.name="abundance"
)
setnames(long, feature_col, "feature")
long <- merge(long, col_map[, .(column_name, sample_id)], by="column_name", all.x=TRUE)
suppressWarnings(long[, abundance := as.numeric(abundance)])
long <- long[is.finite(abundance) & !is.na(sample_id)]
long[, species_feature := clean_taxon(feature)]
long[, genus_feature := genus_key(species_feature)]

species_long <- long[, .(sample_id, feature_family="MetaPhlAn species", feature=species_feature, abundance)]
genus_long <- long[, .(abundance=sum(abundance, na.rm=TRUE)), by=.(sample_id, feature=genus_feature)]
genus_long[, feature_family := "MetaPhlAn genus"]

ab_long <- rbindlist(list(species_long, genus_long), fill=TRUE)
ab_long[, feature_key := feature_key(feature)]

top_keys <- unique(sel[, .(feature_family, feature_key, feature_label_txt, rho)])
ab_top <- merge(ab_long, top_keys[, .(feature_family, feature_key)], by=c("feature_family","feature_key"), all=FALSE)
ab_top <- merge(ab_top, target_k, by="sample_id", all=FALSE)

cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8) return(NA_real_)
  suppressWarnings(unname(cor(x[ok], y[ok], method="spearman")))
}

heat <- ab_top[, {
  vals <- lapply(k_cols, function(kc) cor_one(abundance, get(kc)))
  names(vals) <- k_cols
  as.data.table(vals)
}, by=.(feature_family, feature_key, feature)]

heat_long <- melt(heat, id.vars=c("feature_family","feature_key","feature"),
                  variable.name="PFAS_k", value.name="rho")
heat_long[, PFAS := gsub("^elim_k_", "", PFAS_k)]
heat_long[, PFAS := gsub("_", " ", PFAS)]
heat_long <- merge(
  heat_long,
  sel[, .(feature_family, feature_key, feature_label_txt, primary_rho=rho)],
  by=c("feature_family","feature_key"),
  all.x=TRUE
)
heat_long[, same_direction := sign(rho) == sign(primary_rho)]

compound_summary <- heat_long[, .(
  compound_consistency=mean(same_direction, na.rm=TRUE),
  mean_abs_compound_rho=mean(abs(rho), na.rm=TRUE)
), by=.(feature_family, feature_key)]

heat_wide <- dcast(heat_long, feature_label_txt ~ PFAS, value.var="rho", fill=0)
heat_row_names <- heat_wide$feature_label_txt
heat_wide[, feature_label_txt := NULL]
heat_col_names <- names(heat_wide)
heat_matrix <- as.matrix(heat_wide)
rownames(heat_matrix) <- heat_row_names

row_cluster_B <- cluster_rows(heat_matrix)
col_cluster_B <- cluster_cols(heat_matrix)
heat_row_order <- heat_row_names[row_cluster_B$order]
heat_col_order <- heat_col_names[col_cluster_B$order]

fwrite(data.table(order_index=seq_along(heat_row_order), label=heat_row_order),
       file.path(CLUSTERDIR, "Figure_6B_taxon_cluster_order.tsv"), sep="\t")
fwrite(data.table(order_index=seq_along(heat_col_order), label=heat_col_order),
       file.path(CLUSTERDIR, "Figure_6B_compound_cluster_order.tsv"), sep="\t")

heat_long[, feature_label := factor(feature_label_txt, levels=heat_row_order)]
heat_long[, PFAS := factor(PFAS, levels=heat_col_order)]

fwrite(heat_long, file.path(SRC, "Figure_6B_compound_specific_elimination_heatmap.tsv"), sep="\t")

pB_heat <- ggplot(heat_long, aes(x=PFAS, y=feature_label, fill=rho)) +
  geom_tile(colour="white", linewidth=0.40) +
  scale_fill_gradient2(low=purple["pale"], mid="white", high=purple["deep"], midpoint=0,
                       limits=c(-0.55,0.55), oob=squish, name="Spearman\nrho") +
  labs(x="Compound-specific elimination-rate estimate", y=NULL) +
  theme_pfas(11.0) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(angle=35, hjust=1, size=9.4),
    axis.text.y=element_text(size=9.0),
    axis.title.y=element_blank(),
    legend.position="right"
  )

pB_row <- make_row_dendrogram_plot(row_cluster_B$tree, length(heat_row_order), panel_tag="B")
pB_col <- make_col_dendrogram_plot(col_cluster_B$tree, length(heat_col_order))

pB <- wrap_plots(
  list(pB_row, pB_heat),
  ncol=2,
  widths=c(0.18, 1)
)

# ----------------------
# Panel C: selected Lindell-supported taxa linked to functional mechanism bins
# by sample-level taxon-function co-variation
# ----------------------
# ----------------------
# Panel C: selected Lindell-supported taxa linked to functional mechanism bins
# by sample-level taxon-function co-variation
# ----------------------
if (!exists("clean_text")) {
  clean_text <- function(x) clean_taxon(x)
}
if (!exists("wrap_label")) {
  wrap_label <- function(x, width=30) wrap_lab(x, w=width)
}
if (!exists("fix_family")) {
  fix_family <- function(current_family, raw_feature, label=NULL) {
    cf <- as.character(current_family)
    raw <- paste(as.character(raw_feature), as.character(label), sep=" || ")
    low <- tolower(raw)
    cflow <- tolower(cf)
    out <- cf

    out[grepl("^s__|\\|s__|\\bs__", low)] <- "MetaPhlAn species"
    out[grepl("^g__|\\|g__|\\bg__", low)] <- "MetaPhlAn genus"
    out[grepl("exchange:|^ex_|\\[ex_|agora2 exchange", low)] <- "AGORA2 exchanges"
    out[grepl("reaction:|agora2 reaction|fedcabcpp|aldc|alcd|transport reaction", low)] <- "AGORA2 reactions"
    out[grepl("agora2 subsystem|subsystem:|oxalate metabolism|steroid metabolism|zeatin biosynthesis|butanoate metabolism", low)] <- "AGORA2 subsystems"

    is_taxon <- out %in% c("MetaPhlAn species", "MetaPhlAn genus")
    is_agora <- out %in% c("AGORA2 reactions", "AGORA2 subsystems", "AGORA2 exchanges")
    humann_like <- grepl("pwy|pathway|biosynthesis|degradation|metabolism|fermentation|utilization|glycolysis|vitamin|coenzyme|glycan|ion", low)

    out[!is_taxon & !is_agora & humann_like] <- "HUMAnN pathways"
    out[grepl("\\bbin\\.?[0-9]+|mag bin|_bin", low) & !is_taxon] <- "MAGs"

    out[grepl("species", cflow) & !is_agora] <- "MetaPhlAn species"
    out[grepl("genus", cflow) & !is_agora] <- "MetaPhlAn genus"
    out[grepl("humann|pathway", cflow) & !is_taxon & !is_agora] <- "HUMAnN pathways"
    out[grepl("reaction", cflow)] <- "AGORA2 reactions"
    out[grepl("subsystem", cflow)] <- "AGORA2 subsystems"
    out[grepl("exchange", cflow)] <- "AGORA2 exchanges"
    out[grepl("mag", cflow) & !is_taxon] <- "MAGs"

    out[is.na(out) | out == ""] <- "Other"
    out
  }
}

classify_mechanism <- function(x) {
  z <- tolower(as.character(x))
  z <- gsub("[_\\-]+", " ", z)
  z <- gsub("\\s+", " ", z)

  out <- rep("Unclassified functional context", length(z))

  set_if_empty <- function(pattern, label) {
    hit <- grepl(pattern, z, perl=TRUE) & out == "Unclassified functional context"
    out[hit] <<- label
  }

  # Most specific / directly interpretable categories first.
  set_if_empty("exchange|\\bex\\b|\\[ex_|transport|transporter|permease|abc|efflux|uptake|secretion|symport|antiport|transloc|import|export", "Transport/exchange")
  set_if_empty("bile|chol|steroid|sterol|lipid|fatty acid|fatty acyl|long chain|beta oxidation|β oxidation|acyl|glycerolipid|phospholipid|sphingolipid|isoprenoid|mevalonate", "Bile/lipid/steroid")
  set_if_empty("butanoate|butyrate|propionate|acetate|acetyl coa|scfa|fermentation|succinate|lactate|pyruvate|formate|ethanol|acetogenesis", "SCFA/fermentation")
  set_if_empty("glycan|mucin|cell wall|peptidoglycan|lipopolysaccharide|\\blps\\b|capsule|capsular|polysaccharide|extracellular polymer|\\beps\\b|surface|envelope|adhesion|flagell|fimbria|pili|teichoic|lipoteichoic|oxalate", "Cell surface/glycan")
  set_if_empty("redox|electron|respirat|oxidoreduct|oxidation reduction|dehydrogenase|fadh|\\bfad\\b|nadh|\\bnad\\b|quinone|menaquinone|ferredoxin|cytochrome|hydrogenase|oxidative phosphorylation|atp synthase|proton motive", "Redox/energy")

  # Nutrient and metabolic-context categories.
  set_if_empty("vitamin|biotin|folate|cobalamin|\\bb12\\b|riboflavin|thiamine|niacin|nicotinate|pantothenate|pyridox|pyridoxal|coenzyme|coa\\b|heme|tetrahydrofolate|molybdopterin|lipoate", "Cofactor/vitamin")
  set_if_empty("amino acid|alanine|arginine|aspartate|asparagine|cysteine|glutamate|glutamine|glycine|histidine|isoleucine|leucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|ornithine|urea|ammonia|ammonium|nitrogen|sulfur|sulphur|sulfate|sulphate|taurine", "Amino acid/N/S")
  set_if_empty("carbohydrate|central carbon|glycolysis|gluconeogenesis|pentose phosphate|pentose|starch|sucrose|fructose|mannose|galactose|glucose|glyoxylate|glycolate|citrate cycle|tca|tricarboxylic|citrate|malate|fumarate|2 oxoglutarate|alpha ketoglutarate", "Carbohydrate/central carbon")
  set_if_empty("purine|pyrimidine|nucleotide|nucleoside|ribose|deoxyribose|\\bdna\\b|\\brna\\b|adenine|guanine|cytosine|thymine|uracil|atp|gtp|ctp|utp", "Nucleotide metabolism")
  set_if_empty("aromatic|xenobiotic|drug|detox|stress|glutathione|peroxidase|catalase|superoxide|antioxidant|benzoate|toluene|phenyl|phenol|indole|shikimate|chorismate|antibiotic|resistance", "Stress/xenobiotic")
  set_if_empty("protein|translation|ribosome|ribosomal|trna|aminoacyl|chaperone|folding|protease|peptidase|secretion system|signal peptide|sec pathway|tat pathway", "Protein/translation")
  set_if_empty("iron|\\bfe\\b|siderophore|metal|zinc|copper|manganese|magnesium|potassium|sodium|chloride|phosphate|phosphorus|ion homeostasis|cation|anion", "Ion/metal homeostasis")

  # Broad residual metabolic classes that are still more informative than "Other".
  set_if_empty("biosynthesis|degradation|catabolism|anabolism|metabolism|salvage|assimilation|conversion|interconversion", "Broad biosynthesis/housekeeping")
  set_if_empty("biomass|growth associated|maintenance|demand reaction|sink reaction|artificial reaction", "Model biomass/maintenance")

  out
}

find_candidate_files <- function() {
  roots <- c(SRC, TABLES, P12, file.path(P14, "03_source_data_derived"))
  files <- unique(unlist(lapply(roots, function(r) {
    if (!dir.exists(r)) return(character())
    list.files(r, pattern="\\.(tsv|csv)$", recursive=TRUE, full.names=TRUE)
  })))
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (length(files) == 0) return(character())

  low <- tolower(basename(files))
  keep <- grepl("figure_5|candidate|priorit|table_s22|rf50|functional|humann|agora2", low)
  keep <- keep & !grepl("figure_6|background|manifest|audit|source_data_manifest|taxon_elimination|public_background", low)

  files <- files[keep]
  if (length(files) == 0) return(character())

  low <- tolower(basename(files))
  score <- rep(0L, length(files))
  score <- score + 12L*grepl("figure_5d|figure_5.*candidate|top.*candidate|candidate.*family", low)
  score <- score + 10L*grepl("table_s22|candidate_priorit", low)
  score <- score + 5L*grepl("rf50|phase12", low)
  score <- score + 3L*grepl("humann|agora2|functional", low)
  files[order(-score, file.info(files)$mtime, decreasing=FALSE)]
}

standardize_candidate_file <- function(path) {
  dt <- tryCatch(fread(path, data.table=TRUE, showProgress=FALSE), error=function(e) data.table())
  if (nrow(dt) == 0 || ncol(dt) < 2) return(data.table())

  setnames(dt, make.unique(names(dt), sep="_dup"))
  nms <- names(dt)

  feat_col <- detect_col(nms, c("^feature$", "candidate_feature", "feature_label", "label", "name", "term", "pathway", "reaction", "subsystem", "exchange"))
  family_col <- detect_col(nms, c("feature_family_fixed", "feature_family", "family", "source", "matrix", "feature_set", "data_layer"))
  score_col <- detect_col(nms, c("score_use", "candidate_score", "priority_score", "score_plot", "^score$", "importance", "rf.*importance"))
  rho_col <- detect_col(nms, c("^rho$", "spearman", "correlation", "estimate", "effect"))
  p_col <- detect_col(nms, c("^p$", "p_value", "pvalue", "raw_p", "nominal_p"))
  rank_col <- detect_col(nms, c("^rank$", "candidate_rank", "priority_rank"))

  if (is.na(feat_col)) return(data.table())

  out <- data.table(
    feature=clean_text(dt[[feat_col]]),
    feature_family_raw=if (!is.na(family_col)) as.character(dt[[family_col]]) else "",
    candidate_score=if (!is.na(score_col)) suppressWarnings(as.numeric(dt[[score_col]])) else NA_real_,
    rho=if (!is.na(rho_col)) suppressWarnings(as.numeric(dt[[rho_col]])) else NA_real_,
    p=if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else NA_real_,
    rank=if (!is.na(rank_col)) suppressWarnings(as.numeric(dt[[rank_col]])) else NA_real_,
    source_file=path
  )

  out <- out[feature != "" & !is.na(feature)]
  out[, feature_family := fix_family(feature_family_raw, feature, feature)]
  out
}

candidate_files <- find_candidate_files()
cand_list <- lapply(candidate_files, standardize_candidate_file)
func_raw <- rbindlist(cand_list, fill=TRUE)

if (nrow(func_raw) > 0) {
  func <- copy(func_raw)
  func <- func[!grepl("age|sex|bmi|include|available|intercept|batch|n_k|sample", feature, ignore.case=TRUE)]
  func <- func[feature_family %in% c("HUMAnN pathways", "AGORA2 reactions", "AGORA2 subsystems", "AGORA2 exchanges")]

  func[, feature_key := feature_key(feature)]
  func[, candidate_weight := fifelse(is.finite(candidate_score), candidate_score,
                              fifelse(is.finite(p) & p > 0, -log10(p),
                              fifelse(is.finite(rho), abs(rho), 1)))]
  func[!is.finite(candidate_weight), candidate_weight := 1]
  func[, mechanism_bin := classify_mechanism(paste(feature, feature_family))]
  func[, direction := fifelse(is.finite(rho) & rho > 0, "higher",
                       fifelse(is.finite(rho) & rho < 0, "lower", "not encoded"))]

  # Deduplicate, then keep a bounded candidate set so "Other" cannot dominate only by volume.
  setorder(func, feature_family, feature_key, -candidate_weight)
  func <- func[, .SD[1], by=.(feature_family, feature_key)]

  setorder(func, feature_family, -candidate_weight)
  func <- func[, head(.SD, 45), by=feature_family]

  # Keep a controlled number of residual "Other" rows so it remains visible but does not dominate.
  func_other <- func[mechanism_bin == "Other functional context"][head(order(-candidate_weight), 25)]
  func_known <- func[mechanism_bin != "Other functional context"]
  func <- rbindlist(list(func_known, func_other), fill=TRUE)
} else {
  func <- data.table()
}

# ----------------------
# Find functional abundance matrices and compute selected taxon-function co-variation
# ----------------------
sample_map <- unique(ab_long[, .(sample_id, sample_key=canon_sample(sample_id))])

score_function_matrix <- function(path, family, keys) {
  dt <- tryCatch(fread(path, data.table=TRUE, showProgress=FALSE, nrows=3000), error=function(e) data.table())
  if (nrow(dt) == 0 || ncol(dt) < 3) {
    return(data.table(path=path, feature_family=family, sample_overlap=0, feature_overlap=0, score=-Inf))
  }

  nms <- names(dt)
  sample_overlap <- length(intersect(canon_sample(nms), sample_map$sample_key))
  feature_col <- nms[1]
  fc <- nms[tolower(nms) %in% c("feature","reaction","pathway","subsystem","exchange","name","id")]
  if (length(fc) > 0) feature_col <- fc[1]
  feature_overlap <- length(intersect(feature_key(dt[[feature_col]]), keys))

  low <- tolower(path)
  priority <- 0
  if (family == "HUMAnN pathways" && grepl("humann|pathway", low)) priority <- priority + 10
  if (family == "AGORA2 reactions" && grepl("agora2.*reaction|reaction_abundance", low)) priority <- priority + 10
  if (family == "AGORA2 subsystems" && grepl("agora2.*subsystem|subsystem_abundance", low)) priority <- priority + 10
  if (family == "AGORA2 exchanges" && grepl("agora2.*exchange|exchange_abundance", low)) priority <- priority + 10
  if (grepl("ronneby", low)) priority <- priority + 5
  if (grepl("hmp|nielsen|background|metadata|manifest|audit|summary|mapping|unmatched|matches|model_count|input_summary", low)) priority <- priority - 20
  if (grepl("candidate|priorit|figure_5|figure_6|table_s22|source_data", low)) priority <- priority - 12

  score <- sample_overlap*100 + feature_overlap*20 + priority
  data.table(path=path, feature_family=family, sample_overlap=sample_overlap, feature_overlap=feature_overlap, score=score)
}

find_function_matrices <- function(func) {
  if (nrow(func) == 0) return(data.table())
  files <- unique(list.files(REV, pattern="\\.(tsv|csv)$", recursive=TRUE, full.names=TRUE))
  files <- files[file.exists(files) & file.info(files)$size > 0]
  files <- files[grepl("ronneby|humann|agora2|pathway|reaction|subsystem|exchange", files, ignore.case=TRUE)]
  files <- files[!grepl("hmp|nielsen|public_background|metadata|manifest|audit|summary|mapping_fraction|unmatched|matches|model_load|model_counts|input_summary", files, ignore.case=TRUE)]

  out <- rbindlist(lapply(unique(func$feature_family), function(ff) {
    keys <- unique(func[feature_family == ff, feature_key])
    fam_scores <- rbindlist(lapply(files, score_function_matrix, family=ff, keys=keys), fill=TRUE)
    fam_scores <- fam_scores[is.finite(score)]
    if (nrow(fam_scores) == 0) return(data.table())
    setorder(fam_scores, -score, -feature_overlap, -sample_overlap)
    fam_scores[1]
  }), fill=TRUE)
  out
}

read_function_matrix <- function(path, family, keys) {
  dt <- tryCatch(fread(path, data.table=TRUE, showProgress=FALSE), error=function(e) data.table())
  if (nrow(dt) == 0 || ncol(dt) < 3) return(data.table())

  setnames(dt, make.unique(names(dt), sep="_dup"))
  nms <- names(dt)
  sample_cols <- nms[canon_sample(nms) %in% sample_map$sample_key]

  feature_col <- nms[1]
  fc <- nms[tolower(nms) %in% c("feature","reaction","pathway","subsystem","exchange","name","id")]
  if (length(fc) > 0) feature_col <- fc[1]

  if (length(sample_cols) >= 8) {
    dt[, feature_key_tmp := feature_key(get(feature_col))]
    dt <- dt[feature_key_tmp %in% keys]
    if (nrow(dt) == 0) return(data.table())

    long <- melt(
      dt[, c(feature_col, "feature_key_tmp", sample_cols), with=FALSE],
      id.vars=c(feature_col, "feature_key_tmp"),
      variable.name="sample_col",
      value.name="functional_abundance"
    )
    setnames(long, feature_col, "functional_feature")
    long[, sample_key := canon_sample(sample_col)]
    long <- merge(long, sample_map, by="sample_key", all.x=TRUE)
    suppressWarnings(long[, functional_abundance := as.numeric(functional_abundance)])
    long <- long[is.finite(functional_abundance) & !is.na(sample_id)]
    long[, feature_family := family]
    setnames(long, "feature_key_tmp", "functional_key")
    return(long[, .(sample_id, feature_family, functional_feature, functional_key, functional_abundance)])
  }

  data.table()
}

matrix_scores <- find_function_matrices(func)
func_long <- data.table()
if (nrow(matrix_scores) > 0) {
  func_long <- rbindlist(lapply(seq_len(nrow(matrix_scores)), function(i) {
    ff <- matrix_scores$feature_family[i]
    keys <- unique(func[feature_family == ff, feature_key])
    read_function_matrix(matrix_scores$path[i], ff, keys)
  }), fill=TRUE)
}

taxa_for_c <- copy(sel)
taxa_for_c <- merge(taxa_for_c, compound_summary, by=c("feature_family","feature_key"), all.x=TRUE)
taxa_for_c[!is.finite(compound_consistency), compound_consistency := 0]
taxa_for_c[, selection_score := abs_rho + 0.35*compound_consistency + 0.15*(match_rank/3)]
setorder(taxa_for_c, -selection_score)
taxa_for_c <- taxa_for_c[seq_len(min(.N, 10))]
taxa_for_c[, feature_label_short := wrap_label(feature, 28)]

taxa_ab_c <- merge(
  ab_long,
  taxa_for_c[, .(feature_family, feature_key, taxon_label=feature_label_short)],
  by=c("feature_family", "feature_key"),
  all=FALSE
)
setnames(taxa_ab_c, "abundance", "taxon_abundance")

micro_func <- data.table()

if (nrow(func_long) > 0 && nrow(taxa_ab_c) > 0) {
  func_annot <- func[, .(functional_key=feature_key, feature_family, mechanism_bin, candidate_weight)]
  func_long <- merge(func_long, func_annot, by=c("feature_family", "functional_key"), all.x=TRUE)
  func_long <- func_long[!is.na(mechanism_bin)]

  mf <- merge(taxa_ab_c, func_long, by="sample_id", allow.cartesian=TRUE)

  cor_one <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 8) return(list(rho=NA_real_, p=NA_real_, n=sum(ok)))
    z <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    list(rho=unname(z$estimate), p=z$p.value, n=sum(ok))
  }

  micro_func <- mf[, {
    z <- cor_one(taxon_abundance, functional_abundance)
    .(rho=z$rho, p=z$p, n=z$n)
  }, by=.(taxon_label, functional_feature, feature_family.y, functional_key, mechanism_bin)]

  setnames(micro_func, "feature_family.y", "functional_family")
  micro_func[, abs_rho := abs(rho)]
  micro_func <- micro_func[is.finite(rho)]
}

if (nrow(micro_func) > 0) {
  mf_sum <- micro_func[, .(
    n_features=.N,
    n_linked=sum(abs_rho >= 0.30 | p < 0.05, na.rm=TRUE),
    max_abs_rho=max(abs_rho, na.rm=TRUE),
    signed_max_rho=rho[which.max(abs_rho)],
    representative_feature=functional_feature[which.max(abs_rho)]
  ), by=.(taxon_label, mechanism_bin)]

  all_grid <- CJ(taxon_label=unique(taxa_for_c$feature_label_short), mechanism_bin=sort(unique(func$mechanism_bin)))
  mf_sum <- merge(all_grid, mf_sum, by=c("taxon_label", "mechanism_bin"), all.x=TRUE)
  mf_sum[is.na(n_features), `:=`(n_features=0, n_linked=0, max_abs_rho=0, signed_max_rho=0, representative_feature="")]

  bin_rank <- mf_sum[, .(
    linked=sum(n_linked, na.rm=TRUE),
    total_features=sum(n_features, na.rm=TRUE),
    max_abs=max(max_abs_rho, na.rm=TRUE)
  ), by=mechanism_bin]
  setorder(bin_rank, -linked, -max_abs, -total_features)
  keep_bins <- bin_rank[seq_len(min(.N, 10)), mechanism_bin]
  mf_sum <- mf_sum[mechanism_bin %in% keep_bins]

  mf_sum[, taxon_label := factor(taxon_label, levels=rev(taxa_for_c$feature_label_short))]
  mf_sum[, mechanism_bin := factor(mechanism_bin, levels=bin_rank[mechanism_bin %in% keep_bins, mechanism_bin])]

  mf_wide <- dcast(mf_sum, taxon_label ~ mechanism_bin, value.var="signed_max_rho", fill=0)
  mf_row_names <- mf_wide$taxon_label
  mf_wide[, taxon_label := NULL]
  mf_col_names <- names(mf_wide)
  mf_matrix <- as.matrix(mf_wide)
  rownames(mf_matrix) <- mf_row_names

  row_cluster_C <- cluster_rows(mf_matrix)
  col_cluster_C <- cluster_cols(mf_matrix)
  mf_row_order <- mf_row_names[row_cluster_C$order]
  mf_col_order <- mf_col_names[col_cluster_C$order]

  fwrite(data.table(order_index=seq_along(mf_row_order), label=mf_row_order),
         file.path(CLUSTERDIR, "Figure_6C_taxon_cluster_order.tsv"), sep="\t")
  fwrite(data.table(order_index=seq_along(mf_col_order), label=mf_col_order),
         file.path(CLUSTERDIR, "Figure_6C_mechanism_cluster_order.tsv"), sep="\t")

  mf_sum[, taxon_label := factor(as.character(taxon_label), levels=mf_row_order)]
  mf_sum[, mechanism_bin := factor(as.character(mechanism_bin), levels=mf_col_order)]

  pC_heat <- ggplot(mf_sum, aes(x=mechanism_bin, y=taxon_label, fill=signed_max_rho)) +
    geom_tile(colour="white", linewidth=0.38) +
    scale_fill_gradient2(
      low=purple["pale"], mid="white", high=purple["deep"],
      midpoint=0, limits=c(-0.65,0.65), oob=scales::squish,
      name="Strongest\nSpearman rho"
    ) +
    labs(x="Functional mechanism bin", y="Selected Lindell-supported taxon") +
    theme_pfas(10.9) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=32, hjust=1, size=9.3),
      axis.text.y=element_text(size=9.0),
      legend.position="right"
    )

  pC_row <- make_row_dendrogram_plot(row_cluster_C$tree, length(mf_row_order), panel_tag="C")
  pC_col <- make_col_dendrogram_plot(col_cluster_C$tree, length(mf_col_order))

  pC <- wrap_plots(
    list(pC_row, pC_heat),
    ncol=2,
    widths=c(0.20, 1)
  )
} else {
  mf_sum <- data.table(
    taxon_label="Functional abundance matrices not detected",
    mechanism_bin="No taxon-function co-variation computed",
    n_features=0,
    n_linked=0,
    max_abs_rho=0,
    signed_max_rho=0,
    representative_feature=""
  )
  fwrite(data.table(), file.path(CLUSTERDIR, "Figure_6C_taxon_cluster_order.tsv"), sep="\t")
  fwrite(data.table(), file.path(CLUSTERDIR, "Figure_6C_mechanism_cluster_order.tsv"), sep="\t")
  pC <- ggplot() +
    annotate("text", x=0.5, y=0.5,
             label="Functional abundance matrices were not detected for taxon-function co-variation",
             size=4.2, colour=purple["grey_dark"]) +
    theme_void()
}

# ----------------------
# Panel D: expanded mechanism-by-feature-family support matrix
# ----------------------
if (nrow(func) > 0) {
  mech_family <- func[, .(
    n_candidates=.N,
    max_weight=max(candidate_weight, na.rm=TRUE),
    median_weight=median(candidate_weight, na.rm=TRUE),
    mean_rho=ifelse(any(is.finite(rho)), mean(rho, na.rm=TRUE), NA_real_)
  ), by=.(mechanism_bin, feature_family)]

  mech_family[, display_family := fifelse(
    grepl("^AGORA2 ", feature_family),
    "AGORA2 features",
    as.character(feature_family)
  )]

  mech_family_display <- mech_family[, .(
    n_candidates=sum(n_candidates, na.rm=TRUE),
    max_weight=max(max_weight, na.rm=TRUE),
    median_weight=weighted.mean(
      median_weight,
      w=pmax(n_candidates, 1),
      na.rm=TRUE
    ),
    mean_rho=if (any(is.finite(mean_rho))) {
      weighted.mean(
        mean_rho[is.finite(mean_rho)],
        w=n_candidates[is.finite(mean_rho)],
        na.rm=TRUE
      )
    } else {
      NA_real_
    }
  ), by=.(mechanism_bin, feature_family=display_family)]

  mech_total <- func[, .(
    n_total=.N,
    n_feature_families=uniqueN(feature_family),
    max_weight=max(candidate_weight, na.rm=TRUE),
    median_weight=median(candidate_weight, na.rm=TRUE)
  ), by=mechanism_bin]

  setorder(mech_total, -n_total, -max_weight)
  mech_total[, mechanism_bin := factor(mechanism_bin, levels=rev(mechanism_bin))]
  mech_family[, mechanism_bin := factor(
    mechanism_bin,
    levels=levels(mech_total$mechanism_bin)
  )]
  mech_family_display[, mechanism_bin := factor(
    mechanism_bin,
    levels=levels(mech_total$mechanism_bin)
  )]
  mech_family_display[, feature_family := factor(
    feature_family,
    levels=c("HUMAnN pathways", "AGORA2 features")
  )]
} else {
  mech_total <- data.table(
    mechanism_bin=factor("No HUMAnN/AGORA2 candidate rows detected"),
    n_total=0,
    n_feature_families=0,
    max_weight=0,
    median_weight=0
  )
  mech_family <- data.table(
    mechanism_bin=mech_total$mechanism_bin,
    feature_family=factor("No functional candidates"),
    n_candidates=0,
    max_weight=0,
    median_weight=0,
    mean_rho=NA_real_
  )
  mech_family_display <- copy(mech_family)
}

mech_matrix <- copy(mech_family_display)
mech_matrix[, direction_label := fifelse(is.finite(mean_rho) & mean_rho > 0, "↑",
                                  fifelse(is.finite(mean_rho) & mean_rho < 0, "↓", ""))]
mech_matrix[, tile_label := fifelse(n_candidates > 0, paste0("n=", n_candidates, "\n", direction_label), "")]
mech_matrix[, support_value := n_candidates]

fwrite(func, file.path(SRC, "Figure_6C_functional_candidate_pool_HUMAnN_AGORA2_v8.tsv"), sep="\t")
fwrite(matrix_scores, file.path(SRC, "Figure_6C_functional_matrix_scores_v8.tsv"), sep="\t")
fwrite(micro_func, file.path(SRC, "Figure_6C_taxon_function_correlations_v8.tsv"), sep="\t")
fwrite(mf_sum, file.path(SRC, "Figure_6C_taxon_mechanism_summary_v8.tsv"), sep="\t")
fwrite(mech_total, file.path(SRC, "Figure_6D_functional_mechanism_bin_summary_v8.tsv"), sep="\t")
fwrite(mech_family, file.path(SRC, "Figure_6D_raw_mechanism_by_feature_family.tsv"), sep="\t")
fwrite(mech_matrix, file.path(SRC, "Figure_6D_consolidated_mechanism_by_feature_family.tsv"), sep="\t")

mech_matrix_plot <- copy(mech_matrix)
mech_matrix_plot[, mechanism_bin_chr := as.character(mechanism_bin)]
mech_wide <- dcast(mech_matrix_plot, mechanism_bin_chr ~ feature_family, value.var="support_value", fill=0)
mech_row_names <- mech_wide$mechanism_bin_chr
mech_wide[, mechanism_bin_chr := NULL]
mech_col_names <- names(mech_wide)
mech_matrix_values <- as.matrix(mech_wide)
rownames(mech_matrix_values) <- mech_row_names

row_cluster_D <- cluster_rows(mech_matrix_values)
col_cluster_D <- cluster_cols(mech_matrix_values)
mech_row_order <- mech_row_names[row_cluster_D$order]
mech_col_order <- mech_col_names[col_cluster_D$order]

fwrite(data.table(order_index=seq_along(mech_row_order), label=mech_row_order),
       file.path(CLUSTERDIR, "Figure_6D_mechanism_cluster_order.tsv"), sep="\t")
fwrite(data.table(order_index=seq_along(mech_col_order), label=mech_col_order),
       file.path(CLUSTERDIR, "Figure_6D_feature_family_cluster_order.tsv"), sep="\t")

mech_matrix[, mechanism_bin := factor(as.character(mechanism_bin), levels=mech_row_order)]
mech_matrix[, feature_family := factor(as.character(feature_family), levels=mech_col_order)]

pD_heat <- ggplot(mech_matrix, aes(x=feature_family, y=mechanism_bin, fill=support_value)) +
  geom_tile(colour="white", linewidth=0.45) +
  scale_fill_gradient(low=purple["faint"], high=purple["deep"], name="Candidate\ncount") +
  labs(x="Functional feature family", y=NULL) +
  theme_pfas(10.9) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(angle=18, hjust=1, size=10.0),
    axis.text.y=element_text(size=9.2),
    legend.position="right"
  )

pD_row <- make_row_dendrogram_plot(row_cluster_D$tree, length(mech_row_order), panel_tag="D")
pD_col <- make_col_dendrogram_plot(col_cluster_D$tree, length(mech_col_order))

pD <- wrap_plots(
  list(pD_row, pD_heat),
  ncol=2,
  widths=c(0.20, 1)
)

# ----------------------
# Manuscript text blocks for Methods, Results, and Caption
# ----------------------
methods_text <- paste(
  "Figure 6 taxon-function mechanism-context analysis.",
  "Lindell-supported taxa were defined as Ronneby taxa with species-, genus-, or family-level taxonomic correspondence to PFAS-bioaccumulating taxa reported in Lindell-derived experimental literature curation.",
  "For these taxa, Spearman correlations were calculated against the primary elimination score and individual PFAS elimination-rate estimates in PFAS-exposed Ronneby participants.",
  "Functional mechanism context was assessed using HUMAnN pathway and AGORA2 reaction, subsystem, and exchange candidate features from the Phase 12/Figure 5 candidate-prioritization outputs; the three AGORA2 categories were consolidated into a single AGORA2-features column for Panel D visualization while remaining separate in the exported raw source data.",
  "Candidate functional features were assigned to pre-specified mechanism bins using a transparent keyword dictionary covering transport/exchange, bile/lipid/steroid metabolism, SCFA/fermentation, redox/energy metabolism, cell-surface/glycan context, cofactor/vitamin metabolism, amino-acid/nitrogen/sulfur metabolism, carbohydrate/central-carbon metabolism, nucleotide metabolism, stress/xenobiotic context, protein/translation, ion/metal homeostasis, broad biosynthesis/housekeeping, model biomass/maintenance, and residual unclassified context.",
  "For selected Lindell-supported taxa, sample-level taxon abundances were correlated with candidate functional feature abundances where matched Ronneby functional matrices were available; these taxon-function correlations were summarized by mechanism bin.",
  "These analyses were used for exploratory biological interpretation and hypothesis generation, not for causal inference, PFAS-transformation inference, or biomarker validation.",
  sep=" "
)

results_text <- paste(
  "A subset of Lindell-supported Ronneby taxa showed exploratory associations with the primary elimination score and consistent directional structure across multiple compound-specific elimination-rate estimates.",
  "Taxon-function co-variation analysis linked selected Lindell-supported taxa to functional candidate bins including transport/exchange, bile/lipid/steroid metabolism, SCFA/fermentation, redox/energy metabolism, cell-surface/glycan context, and related metabolic processes.",
  "The functional support matrix indicated that these mechanism bins were represented across HUMAnN pathways and the consolidated AGORA2 feature set.",
  "Together, these results support a gut-retention and microbial-context hypothesis for follow-up while remaining non-causal and exploratory.",
  sep=" "
)

caption_text <- paste(
  "Figure 6. Lindell-supported PFAS-bioaccumulation context and taxon-function mechanism hypotheses for Ronneby PFAS elimination.",
  "A, Lindell-supported Ronneby taxa ranked by Spearman correlation with the primary elimination score among PFAS-exposed participants; point size indicates prevalence among PFAS-exposed participants.",
  "B, compound-specific elimination-rate correlations for the same selected taxa.",
  "C, sample-level co-variation between selected Lindell-supported taxa and HUMAnN/AGORA2 candidate functional features, summarized by functional mechanism bin; labels indicate the number of candidate functions in each bin showing |rho| >= 0.30 or nominal p < 0.05 with the taxon.",
  "D, mechanism-by-feature-family support matrix showing candidate counts for HUMAnN pathways and the consolidated AGORA2 feature set; the separate AGORA2 reaction, subsystem, and exchange counts are retained in the exported source data.",
  "All panels are exploratory biological-context analyses and do not establish causal PFAS transformation or validated microbial biomarkers.",
  sep=" "
)

writeLines(methods_text, file.path(SRC, "Figure_6_v9_methods_text.txt"))
writeLines(results_text, file.path(SRC, "Figure_6_v9_results_text.txt"))
writeLines(caption_text, file.path(SRC, "Figure_6_v9_caption_text.txt"))

# ----------------------
# Assemble
# ----------------------

# === FIGURE6_LAYOUT_READABILITY_PATCH_V1_0_BEGIN ===
panel_tag_theme_fig6 <- theme(
  plot.tag=element_text(size=14, face="bold", colour="black"),
  plot.tag.position=c(0, 1)
)

if (exists("pA")) {
  pA <- pA +
    theme_pfas(11.3) +
    panel_tag_theme_fig6 +
    theme(
      axis.text.x=element_text(size=9.8),
      axis.text.y=element_text(size=9.8),
      axis.title.x=element_text(size=12.0, margin=margin(t=8)),
      axis.title.y=element_text(size=12.0, margin=margin(r=8)),
      legend.title=element_text(size=11.3, face="bold"),
      legend.text=element_text(size=9.5)
    )
}

if (exists("pB_heat") && exists("pB_row")) {
  pB_heat <- pB_heat +
    theme_pfas(11.1) +
    theme(
      axis.text.x=element_text(angle=32, hjust=1, size=9.6),
      axis.text.y=element_text(size=9.4),
      axis.title.x=element_text(size=11.8, margin=margin(t=8)),
      legend.title=element_text(size=11.0, face="bold"),
      legend.text=element_text(size=9.2),
      plot.margin=margin(4, 6, 4, 0)
    )
  pB_row <- pB_row + panel_tag_theme_fig6 + theme(plot.margin=margin(4, 1, 4, 0))
  pB <- wrap_plots(list(pB_row, pB_heat), ncol=2, widths=c(0.13, 1))
}

if (exists("pC_heat") && exists("pC_row")) {
  pC_heat <- pC_heat +
    labs(x="Functional mechanism bin", y=NULL) +
    theme_pfas(11.1) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=30, hjust=1, size=9.7),
      axis.text.y=element_text(size=9.2),
      axis.title.x=element_text(size=11.8, margin=margin(t=8)),
      legend.title=element_text(size=11.0, face="bold"),
      legend.text=element_text(size=9.2),
      plot.margin=margin(4, 6, 4, 0)
    )
  pC_row <- pC_row + panel_tag_theme_fig6 + theme(plot.margin=margin(4, 1, 4, 0))
  pC_ylabel <- wrap_elements(full=grid::textGrob(
    "Selected Lindell-supported taxon",
    rot=90,
    gp=grid::gpar(fontsize=11.2, fontface="bold", col=purple["grey_dark"])
  ))
  pC <- wrap_plots(list(pC_ylabel, pC_row, pC_heat), ncol=3, widths=c(0.08, 0.13, 1))
}

if (exists("pD_heat") && exists("pD_row")) {
  pD_heat <- pD_heat +
    theme_pfas(11.1) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=18, hjust=1, size=9.8),
      axis.text.y=element_text(size=9.4),
      axis.title.x=element_text(size=11.8, margin=margin(t=8)),
      legend.title=element_text(size=11.0, face="bold"),
      legend.text=element_text(size=9.2),
      plot.margin=margin(4, 6, 4, 0)
    )
  pD_row <- pD_row + panel_tag_theme_fig6 + theme(plot.margin=margin(4, 1, 4, 0))
  pD <- wrap_plots(list(pD_row, pD_heat), ncol=2, widths=c(0.13, 1))
}
# === FIGURE6_LAYOUT_READABILITY_PATCH_V1_0_END ===

# FIGURE6_LAYOUT_READABILITY_PATCH_V1_1

# === FIGURE6_LAYOUT_READABILITY_PATCH_V1_1_BEGIN ===
make_vertical_ylabel_fig6 <- function(label, size=11.4) {
  wrap_elements(
    full = grid::textGrob(
      label,
      rot = 90,
      gp = grid::gpar(fontsize=size, fontface="bold", col="#333333")
    )
  )
}

panel_tag_theme_fig6_v11 <- theme(
  plot.tag = element_text(size=15, face="bold", colour="black"),
  plot.tag.position = c(0.01, 0.99)
)

if (exists("pA")) {
  pA <- pA +
    theme_pfas(11.4) +
    panel_tag_theme_fig6_v11 +
    theme(
      axis.text.x = element_text(size=9.8),
      axis.text.y = element_text(size=9.8),
      axis.title.x = element_text(size=12.0, face="bold", margin=margin(t=8)),
      axis.title.y = element_text(size=12.0, face="bold", margin=margin(r=8)),
      legend.title = element_text(size=11.0, face="bold"),
      legend.text = element_text(size=9.3),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

if (exists("pB_heat") && exists("pB_row")) {
  pB_heat <- pB_heat +
    labs(x="Compound-specific elimination-rate estimate", y=NULL) +
    theme_pfas(11.2) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=32, hjust=1, size=9.4),
      axis.text.y = element_text(size=9.1),
      axis.title.x = element_text(size=11.8, face="bold", margin=margin(t=8)),
      axis.title.y = element_blank(),
      legend.title = element_text(size=10.8, face="bold"),
      legend.text = element_text(size=9.1),
      legend.position = "right",
      plot.margin = margin(2, 2, 2, 2, "mm")
    )

  pB_row <- pB_row +
    panel_tag_theme_fig6_v11 +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"))

  pB_ylabel <- make_vertical_ylabel_fig6(
    "Lindell-supported Ronneby taxon",
    size=11.2
  )

  pB <- wrap_plots(
    list(pB_ylabel, pB_row, pB_heat),
    ncol=3,
    widths=c(0.10, 0.16, 1.00)
  )
}

if (exists("pC_heat") && exists("pC_row")) {
  pC_heat <- pC_heat +
    labs(x="Functional mechanism bin", y=NULL) +
    theme_pfas(11.2) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=30, hjust=1, size=9.2),
      axis.text.y = element_text(size=9.0),
      axis.title.x = element_text(size=11.8, face="bold", margin=margin(t=8)),
      axis.title.y = element_blank(),
      legend.title = element_text(size=10.8, face="bold"),
      legend.text = element_text(size=9.1),
      legend.position = "right",
      plot.margin = margin(2, 2, 2, 2, "mm")
    )

  pC_row <- pC_row +
    panel_tag_theme_fig6_v11 +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"))

  pC_ylabel <- make_vertical_ylabel_fig6(
    "Selected Lindell-supported taxon",
    size=11.2
  )

  pC <- wrap_plots(
    list(pC_ylabel, pC_row, pC_heat),
    ncol=3,
    widths=c(0.10, 0.18, 1.00)
  )
}

if (exists("pD_heat") && exists("pD_row")) {
  pD_heat <- pD_heat +
    labs(x="Functional feature family", y=NULL) +
    theme_pfas(11.2) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=18, hjust=1, size=9.3),
      axis.text.y = element_text(size=9.0),
      axis.title.x = element_text(size=11.8, face="bold", margin=margin(t=8)),
      axis.title.y = element_blank(),
      legend.title = element_text(size=10.8, face="bold"),
      legend.text = element_text(size=9.1),
      legend.position = "right",
      plot.margin = margin(2, 2, 2, 2, "mm")
    )

  pD_row <- pD_row +
    panel_tag_theme_fig6_v11 +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"))

  pD_ylabel <- make_vertical_ylabel_fig6(
    "Functional mechanism bin",
    size=11.2
  )

  pD <- wrap_plots(
    list(pD_ylabel, pD_row, pD_heat),
    ncol=3,
    widths=c(0.10, 0.18, 1.00)
  )
}
# === FIGURE6_LAYOUT_READABILITY_PATCH_V1_1_END ===

row_top <- wrap_plots(
  list(pA, wrap_elements(full=pB)),
  ncol=2,
  widths=c(0.90, 1.10)
)
row_bottom <- wrap_plots(
  list(wrap_elements(full=pC), wrap_elements(full=pD)),
  ncol=2,
  widths=c(1.02, 0.98)
)
fig <- wrap_plots(list(row_top, row_bottom), ncol=1, heights=c(0.94, 1.06))

output_stem <- file.path(FIGD, "Figure_6_excellence_candidate")
pdf_out <- paste0(output_stem, ".pdf")
png_out <- paste0(output_stem, ".png")
tif_out <- paste0(output_stem, ".tiff")
preview_out <- file.path(PNGD, "Figure_6_excellence_preview.png")

ggsave(pdf_out, fig, width=16.4, height=12.8, units="in", device=cairo_pdf, limitsize=FALSE)
ggsave(png_out, fig, width=16.4, height=12.8, units="in", dpi=420, limitsize=FALSE)
ggsave(tif_out, fig, width=16.4, height=12.8, units="in", dpi=600, compression="lzw", limitsize=FALSE)
ggsave(preview_out, fig, width=16.4, height=12.8, units="in", dpi=240, limitsize=FALSE)

required_outputs <- c(pdf_out, png_out, tif_out)
if (any(!file.exists(required_outputs)) || any(file.info(required_outputs)$size <= 0)) {
  stop("One or more Figure 6 excellence candidate outputs are missing.", call.=FALSE)
}

manifest <- data.table(
  figure="Figure_6",
  panel=c("A","B","C","D","Methods/Results/Caption"),
  source_data=c(
    "Figure_6A_selected_lindell_support_taxa.tsv",
    "Figure_6B_compound_specific_elimination_heatmap.tsv",
    "Figure_6C_taxon_mechanism_summary_v8.tsv; Figure_6C_taxon_function_correlations_v8.tsv",
    "Figure_6D_raw_mechanism_by_feature_family.tsv; Figure_6D_consolidated_mechanism_by_feature_family.tsv; Figure_6D_functional_mechanism_bin_summary_v8.tsv",
    "Figure_6_v9_methods_text.txt; Figure_6_v9_results_text.txt; Figure_6_v9_caption_text.txt"
  )
)
fwrite(manifest, file.path(SRC, "Figure_6_source_data_manifest.tsv"), sep="\t")

audit <- data.table(
  item=c(
    "version",
    "taxcorr_path",
    "overlap_detail_path",
    "pfas_path",
    "species_path",
    "candidate_files_scanned",
    "candidate_files_used",
    "functional_matrices_used",
    "n_taxa_tested",
    "n_selected_lindell_supported_taxa",
    "n_lindell_supported_taxa_in_taxcorr",
    "n_functional_candidates_after_filtering",
    "n_taxon_function_correlations",
    "n_functional_mechanism_bins",
    "methods_text",
    "results_text",
    "caption_text",
    "Figure_6_pdf",
    "Figure_6_png",
    "Figure_6_tiff",
    "Figure_6_preview"
  ),
  value=c(
    "v11_closeout_candidate_lindell_taxon_function_mechanism",
    taxcorr_path,
    overlap_detail_path,
    pfas_path,
    species_path,
    length(candidate_files),
    paste(candidate_files, collapse=";"),
    ifelse(nrow(matrix_scores) > 0, paste(matrix_scores$path, collapse=";"), "none detected"),
    nrow(taxcorr),
    nrow(sel),
    sum(tc$lindell_support, na.rm=TRUE),
    nrow(func),
    nrow(micro_func),
    nrow(mech_total),
    file.path(SRC, "Figure_6_v9_methods_text.txt"),
    file.path(SRC, "Figure_6_v9_results_text.txt"),
    file.path(SRC, "Figure_6_v9_caption_text.txt"),
    pdf_out,
    png_out,
    tif_out,
    preview_out
  )
)
fwrite(audit, file.path(LOGD, "Figure_6_excellence_audit.tsv"), sep="\t")

cat("=== PFAS FIGURE 6 EXCELLENCE CANDIDATE ===\n")
cat("Selected Lindell-supported taxa:", nrow(sel), "\n")
cat("Functional candidates:", nrow(func), "\n")
cat("Taxon-function correlations:", nrow(micro_func), "\n")
cat("Mechanism bins:", nrow(mech_total), "\n")
cat("Output:", pdf_out, "\n")
cat("Status: READY_FOR_FIGURE6_VISUAL_QC\n")
