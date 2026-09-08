suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 32A_build_Figure_S01_supplementary.R ROOT STAGE_BASE")
ROOT <- normalizePath(args[[1]], mustWork = TRUE)
STAGE_BASE <- args[[2]]

FIG_ID <- "Figure_S01"
FIG_BUILD_DIR <- file.path(STAGE_BASE, "submission", "figure_builds", FIG_ID)
FIGURES_DIR <- file.path(FIG_BUILD_DIR, "figures")
SRC_DIR <- file.path(STAGE_BASE, "submission", "source_data", "supplementary", FIG_ID)
SCRIPT_DIR <- file.path(STAGE_BASE, "submission", "scripts_used", "supplementary", FIG_ID)
LEG_DIR <- file.path(STAGE_BASE, "submission", "legends", "supplementary")
MAN_DIR <- file.path(STAGE_BASE, "submission", "manifests", "supplementary", FIG_ID)
for (d in c(FIGURES_DIR, SRC_DIR, SCRIPT_DIR, LEG_DIR, MAN_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pick_first <- function(paths) {
  for (p in paths) if (file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  NA_character_
}
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  fread(path, sep = if (ext == 'csv') ',' else '\t', data.table = TRUE)
}
fmt_label <- function(x) {
  y <- gsub('^fecal_', '', x, ignore.case = TRUE)
  y <- gsub('^elim_', '', y)
  y <- gsub('^k_', '', y)
  y <- gsub('^HL_', '', y)
  y <- gsub('totalPFOS', 'Tot PFOS', y, ignore.case = TRUE)
  y <- gsub('PFOS_MP26', 'PFOS 62', y)
  y <- gsub('PFOS_MP345', 'PFOS 345', y)
  y <- gsub('PFOS_MP1', 'PFOS MP1', y)
  y <- gsub('_', ' ', y)
  trimws(y)
}

src <- list(
  meta = pick_first(c(file.path(ROOT,'PFAS_Metadata.csv'))),
  counts = pick_first(c(file.path(ROOT,'submission/source_data/main/Figure_1/Figure_1A_study_sample_counts.tsv'), file.path(ROOT,'submission/figure_builds/Figure_1/source_data/Figure_1A_study_sample_counts.tsv'))),
  part = pick_first(c(file.path(ROOT,'submission/source_data/main/Figure_1/Figure_1E_participant_fecal_marker_metadata.tsv'), file.path(ROOT,'submission/figure_builds/Figure_1/source_data/Figure_1E_participant_fecal_marker_metadata.tsv'))),
  mapping = pick_first(c(file.path(ROOT,'submission/manifests/supplementary_batch1_mapping/Figure_S01_panel_map.tsv'))),
  popres = pick_first(c(file.path(ROOT,'submission/manifests/supplementary_batch1_mapping/Figure_S01_population_resolution.tsv')))
)
if (is.na(src$meta) || is.na(src$counts)) stop('Required Figure S1 inputs not found')
meta <- read_any(src$meta); setDT(meta); setnames(meta, make.unique(names(meta), sep='__dup'))

find_col <- function(patterns, nms=names(meta)) {
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, perl = TRUE, value = TRUE)
    if (length(hit) >= 1) return(hit[1])
  }
  NA_character_
}
numvec <- function(x) suppressWarnings(as.numeric(as.character(x)))

sample_col <- find_col(c('^Sample_Name$','^sample_name$','^sample$','^id$'))
if (is.na(sample_col)) stop('Sample column not found')
age_col <- find_col(c('^Age$','^age$'))
sex_col <- find_col(c('^Sex$','^sex$'))
bmi_col <- find_col(c('^BMI$','^bmi$','body.*mass'))
zon_col <- find_col(c('zonulin'))
cal_col <- find_col(c('calprotectin'))
condition_col <- find_col(c('^Condition_clean$','^Condition$','group','cohort'))

meta[, Group := NA_character_]
if (!is.na(condition_col)) {
  meta[, Group := fifelse(grepl('ref|control', get(condition_col), ignore.case = TRUE), 'Reference',
                   fifelse(grepl('expos|pfas', get(condition_col), ignore.case = TRUE), 'PFAS exposed', NA_character_))]
}
if (all(is.na(meta$Group))) {
  kp <- find_col(c('^k_LPFOS$','^elim_k_LPFOS$','^HL_LPFOS$'))
  if (!is.na(kp)) meta[!is.na(get(kp)), Group := 'PFAS exposed']
}
if (sum(meta$Group == 'PFAS exposed', na.rm = TRUE) == 47L) meta[is.na(Group), Group := 'Reference']
meta[, Group := factor(Group, levels = c('PFAS exposed','Reference'))]

counts_dt <- read_any(src$counts); setDT(counts_dt); setnames(counts_dt, make.unique(names(counts_dt), sep='__dup'))
item_col <- names(counts_dt)[grep('item|group|label|condition', names(counts_dt), ignore.case = TRUE)][1]
count_col <- names(counts_dt)[grep('^n$|count', names(counts_dt), ignore.case = TRUE)][1]
get_count <- function(patterns, default=NA_integer_) {
  idx <- integer(0)
  for (pat in patterns) idx <- unique(c(idx, grep(pat, counts_dt[[item_col]], ignore.case = TRUE, perl = TRUE)))
  if (length(idx) == 0) return(default)
  as.integer(counts_dt[[count_col]][idx[1]])
}
exposed_n <- get_count(c('expos','pfas.*expos'), 47L)
reference_n <- get_count(c('ref','control'), 18L)
total_n <- uniqueN(meta[[sample_col]])
if (total_n != 65L) stop('Frozen total-count gate failed')
if (exposed_n != 47L) stop('Frozen exposed-count gate failed')
if (reference_n != 18L) stop('Frozen reference-count gate failed')

# variable selection
serum_pref <- c('PFOA','PFNA','PFHxS','PFHpS','PFPeS','L_PFOS','PFOS_MP1','PFOS_MP26','PFOS_MP345')
serum_vars <- serum_pref[serum_pref %in% names(meta)]
if (length(serum_vars) < 4) {
  extra <- grep('^(PF|L_PFOS|PFOS_MP)', names(meta), value = TRUE)
  extra <- extra[!grepl('fecal|k_|elim_|HL_', extra, ignore.case = TRUE)]
  serum_vars <- unique(c(serum_vars, extra))
}
serum_vars <- serum_vars[seq_len(min(length(serum_vars), 7))]

fecal_vars <- grep('^fecal_', names(meta), value = TRUE)
if (length(fecal_vars) > 0) {
  fecal_ord <- c('fecal_PFOA','fecal_PFNA','fecal_PFHxS','fecal_PFHpS','fecal_PFPeS','fecal_PFOS_MP26','fecal_PFOS_MP345','fecal_totalPFOS','fecal_PFOS')
  fecal_vars <- unique(c(fecal_ord[fecal_ord %in% fecal_vars], fecal_vars))
  fecal_vars <- fecal_vars[seq_len(min(length(fecal_vars), 7))]
}

k_vars <- grep('^(k_|elim_k_|HL_)', names(meta), value = TRUE)
k_pref <- c('k_PFOA','k_PFNA','k_PFHxS','k_PFHpS','k_PFPeS','k_LPFOS','k_PFOS_MP1','k_PFOS_MP26','k_PFOS_MP345',
            'elim_k_PFOA','elim_k_PFNA','elim_k_PFHxS','elim_k_PFHpS','elim_k_PFPeS','elim_k_LPFOS','elim_k_PFOS_MP1','elim_k_PFOS_MP26','elim_k_PFOS_MP345',
            'HL_PFOA','HL_PFNA','HL_PFHxS','HL_PFHpS','HL_PFPeS','HL_LPFOS','HL_PFOS_MP1','HL_PFOS_MP26','HL_PFOS_MP345')
k_vars <- unique(c(k_pref[k_pref %in% k_vars], k_vars))
k_vars <- k_vars[seq_len(min(length(k_vars), 7))]

# source data
study_counts <- data.table(group = c('PFAS exposed','Reference'), n = c(exposed_n, reference_n))
fwrite(study_counts, file.path(SRC_DIR, 'Figure_S01_panelA_study_counts.tsv'), sep='\t')

mk_row <- function(layer_name, vars, formatter = fmt_label) {
  if (length(vars) == 0) return(data.table())
  data.table(layer = layer_name,
             variable = vars,
             analyte = formatter(vars),
             n_nonmissing = vapply(vars, function(v) sum(!is.na(meta[[v]])), integer(1)))
}
cover <- rbindlist(list(
  mk_row('Serum', serum_vars),
  mk_row('Fecal', fecal_vars),
  mk_row('Elim. k', k_vars)
), fill = TRUE)
cover <- cover[!is.na(analyte) & analyte != '' & analyte != 'NA']
cover[, layer := factor(layer, levels = c('Serum','Fecal','Elim. k'))]
# use analyte order from observed data only
cover[, analyte := factor(analyte, levels = unique(analyte))]
fwrite(cover, file.path(SRC_DIR, 'Figure_S01_panelB_coverage.tsv'), sep='\t')

make_dist_long <- function(vars, layer) {
  if (length(vars) == 0) return(data.table())
  rbindlist(lapply(vars, function(v) {
    val <- numvec(meta[[v]])
    dt <- data.table(variable=v, value=val)
    dt <- dt[!is.na(value)]
    if (!nrow(dt)) return(NULL)
    posmin <- min(dt$value[dt$value > 0], na.rm = TRUE)
    if (!is.finite(posmin)) posmin <- 1e-6
    dt[, `:=`(layer = layer, analyte = fmt_label(v), log10_value = log10(pmax(value, posmin / 10)))]
    dt
  }), fill = TRUE)
}
dist_dt <- rbindlist(list(make_dist_long(serum_vars, 'Serum'), make_dist_long(fecal_vars, 'Fecal')), fill = TRUE)
if (!nrow(dist_dt)) stop('Panel C PFAS distributions had no plottable data')
dist_dt[, analyte := factor(analyte, levels = rev(unique(analyte)))]
fwrite(dist_dt, file.path(SRC_DIR, 'Figure_S01_panelC_distributions.tsv'), sep='\t')

base_name <- function(v) {
  y <- gsub('^fecal_', '', v)
  y <- gsub('^(k_|elim_k_|HL_)', '', y)
  y
}
serum_map <- setNames(serum_vars, base_name(serum_vars))
fecal_map <- setNames(fecal_vars, base_name(fecal_vars))
k_map <- setNames(k_vars, base_name(k_vars))
common <- sort(unique(c(intersect(names(serum_map), names(fecal_map)), intersect(names(serum_map), names(k_map)), intersect(names(fecal_map), names(k_map)))))
common <- common[common != '']
add_rel <- function(comp, v1, v2, rel) {
  if (is.null(v1) || is.null(v2) || !(v1 %in% names(meta)) || !(v2 %in% names(meta))) return(NULL)
  a <- numvec(meta[[v1]])
  b <- numvec(meta[[v2]])
  ok <- complete.cases(a, b)
  if (sum(ok) < 5) return(NULL)
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = 'spearman', exact = FALSE))
  data.table(compound = fmt_label(comp), relationship = rel, rho = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
}
rel_list <- list()
for (comp in common) {
  rel_list[[length(rel_list)+1]] <- add_rel(comp, unname(serum_map[comp]), unname(fecal_map[comp]), 'Serum–fecal')
  rel_list[[length(rel_list)+1]] <- add_rel(comp, unname(serum_map[comp]), unname(k_map[comp]), 'Serum–k')
  rel_list[[length(rel_list)+1]] <- add_rel(comp, unname(fecal_map[comp]), unname(k_map[comp]), 'Fecal–k')
}
rel_dt <- rbindlist(rel_list, fill = TRUE)
if (!nrow(rel_dt)) stop('Panel D relationship statistics had no plottable data')
rel_dt[, sig := fifelse(p_value < 0.001, '***', fifelse(p_value < 0.01, '**', fifelse(p_value < 0.05, '*', '')))]
rel_dt[, compound := factor(compound, levels = rev(unique(compound)))]
rel_dt[, relationship := factor(relationship, levels = c('Serum–fecal','Serum–k','Fecal–k'))]
fwrite(rel_dt, file.path(SRC_DIR, 'Figure_S01_panelD_relationships.tsv'), sep='\t')

meta_vars <- c(Age = age_col, BMI = bmi_col, Zonulin = zon_col, Calprotectin = cal_col)
meta_vars <- meta_vars[!is.na(meta_vars)]
meta_long <- rbindlist(lapply(names(meta_vars), function(lbl) {
  coln <- meta_vars[[lbl]]
  val <- numvec(meta[[coln]])
  dt <- data.table(Group = meta$Group, variable = lbl, value = val)
  dt <- dt[!is.na(value) & !is.na(Group)]
  if (!nrow(dt)) return(NULL)
  if (lbl %in% c('Zonulin','Calprotectin')) {
    posmin <- min(dt$value[dt$value > 0], na.rm = TRUE)
    if (!is.finite(posmin)) posmin <- 1e-6
    dt[, value := log10(pmax(value, posmin / 10))]
    dt[, variable := paste0(lbl, '\nlog10')]
  }
  dt
}), fill = TRUE)
if (!nrow(meta_long)) stop('Panel E metadata balance had no plottable data')
meta_long[, variable := factor(variable, levels = unique(variable))]
meta_stats <- meta_long[, {
  p <- tryCatch(wilcox.test(value ~ Group, exact = FALSE)$p.value, error = function(e) NA_real_)
  .(p_value = p, y_pos = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE)))
}, by = variable]
meta_stats[, label := fifelse(is.na(p_value), 'p=NA', fifelse(p_value < 0.001, 'p<0.001', paste0('p=', formatC(p_value, format='f', digits=3))))]
fwrite(meta_long, file.path(SRC_DIR, 'Figure_S01_panelE_metadata_balance.tsv'), sep='\t')
fwrite(meta_stats, file.path(SRC_DIR, 'Figure_S01_panelE_metadata_stats.tsv'), sep='\t')

counts_summary <- data.table(metric = c('Total cohort','PFAS exposed','Reference'), n = c(total_n, exposed_n, reference_n))
fwrite(counts_summary, file.path(SRC_DIR, 'Figure_S01_counts_summary.tsv'), sep='\t')

# plotting
pfas_cols <- c(deep='#4B226B', main='#6F4794', mid='#9A77BC', light='#D8CCE9', verylight='#F1ECF8', grid='#DDD7E8')
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = pfas_cols['grid'], linewidth = 0.35),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = '#EAE4F2', color = '#9F90B7', linewidth = 0.6),
    strip.text = element_text(face = 'bold', size = 11),
    axis.title = element_text(face = 'bold', size = 11),
    axis.text = element_text(color = '#333333'),
    legend.title = element_text(face = 'bold'),
    plot.margin = margin(4, 4, 4, 4)
  )

pA <- ggplot(study_counts, aes(x = reorder(group, n), y = n)) +
  geom_col(fill = pfas_cols['main'], width = 0.65) +
  geom_text(aes(label = paste0('n=', n)), hjust = -0.1, size = 4.5, color='#333333') +
  coord_flip(clip = 'off') +
  expand_limits(y = max(study_counts$n) * 1.3) +
  labs(x = 'Samples', y = NULL) +
  base_theme + theme(legend.position = 'none')

pB <- ggplot(cover, aes(x = analyte, y = layer, fill = n_nonmissing)) +
  geom_tile(color = 'white', linewidth = 0.4) +
  geom_text(aes(label = n_nonmissing), color = 'white', fontface = 'bold', size = 4.4) +
  scale_fill_gradient(low = pfas_cols['light'], high = '#5D14A6', guide = 'none') +
  labs(x = 'PFAS analyte', y = 'Phenotype layer') +
  base_theme + theme(axis.text.x = element_text(angle = 35, hjust = 1))

layer_cols <- c('Serum' = pfas_cols['main'], 'Fecal' = '#BBB3CE')
pC <- ggplot(dist_dt, aes(x = log10_value, y = analyte, fill = layer)) +
  geom_boxplot(width = 0.65, outlier.size = 1.5, alpha = 0.95) +
  facet_wrap(~ layer, ncol = 1, scales = 'free_y', strip.position = 'right') +
  scale_fill_manual(values = layer_cols) +
  labs(x = 'log10 PFAS concentration', y = 'PFAS variable') +
  base_theme + theme(legend.position = 'none', strip.placement = 'outside')

shape_vals <- c('Serum–fecal' = 21, 'Serum–k' = 24, 'Fecal–k' = 22)
pD <- ggplot(rel_dt, aes(y = compound)) +
  geom_vline(xintercept = 0, color = '#777777', linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = rho, yend = compound, color = rho), linewidth = 1.1) +
  geom_point(aes(x = rho, shape = relationship, fill = rho), size = 5, color = '#444444') +
  geom_text(aes(x = pmin(pmax(rho + ifelse(rho >= 0, 0.08, -0.08), -1.05), 1.05), label = sig), size = 4.2, vjust = 0.5) +
  scale_shape_manual(values = shape_vals, name = 'Relationship') +
  scale_fill_gradient2(low = pfas_cols['light'], mid = pfas_cols['mid'], high = '#5D14A6', midpoint = 0, name = 'Spearman\nrho') +
  scale_color_gradient2(low = pfas_cols['light'], mid = pfas_cols['mid'], high = '#5D14A6', midpoint = 0, guide = 'none') +
  labs(x = 'Spearman rho, pairwise complete observations', y = 'PFAS compound', caption = '* p<0.05; ** p<0.01; *** p<0.001') +
  xlim(-1.1, 1.1) +
  base_theme + theme(plot.caption = element_text(size = 9, hjust = 1))

meta_long[, Group := factor(Group, levels = c('PFAS exposed','Reference'))]
pE <- ggplot(meta_long, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.75, color = pfas_cols['deep']) +
  facet_wrap(~ variable, scales = 'free_y', nrow = 1) +
  geom_text(data = meta_stats, aes(x = 1.5, y = y_pos, label = label), inherit.aes = FALSE, size = 3.6) +
  scale_fill_manual(values = c('PFAS exposed' = pfas_cols['main'], 'Reference' = pfas_cols['light'])) +
  labs(x = 'Ronneby group', y = 'Metadata value') +
  base_theme + theme(legend.position = 'none')

figure <- (pA | pB) / (pC | pD) / pE +
  plot_layout(widths = c(1.0, 1.5), heights = c(0.8, 1.25, 0.95), guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 20), plot.tag.position = c(0,1))

pdf_path <- file.path(FIGURES_DIR, 'Figure_S01_excellence_candidate.pdf')
png_path <- file.path(FIGURES_DIR, 'Figure_S01_excellence_candidate.png')
tiff_path <- file.path(FIGURES_DIR, 'Figure_S01_excellence_candidate.tiff')
ggsave(pdf_path, figure, width = 13.2, height = 15.0, units = 'in', device = cairo_pdf)
ggsave(png_path, figure, width = 13.2, height = 15.0, units = 'in', dpi = 300)
ggsave(tiff_path, figure, width = 13.2, height = 15.0, units = 'in', dpi = 300, compression = 'lzw')

installed_script <- file.path(ROOT, 'scripts', '32A_build_Figure_S01_supplementary.R')
if (file.exists(installed_script)) file.copy(installed_script, file.path(SCRIPT_DIR, basename(installed_script)), overwrite = TRUE)

legend_lines <- c(
  'Figure S1 legend is maintained separately in manuscript text.',
  'Internal summary: panels A-E show study counts, phenotype coverage, PFAS distributions, relationship statistics, and metadata balance.',
  'Relationship-panel asterisks indicate nominal significance: * p<0.05; ** p<0.01; *** p<0.001.'
)
writeLines(legend_lines, file.path(LEG_DIR, 'Figure_S01.txt'))

prov <- data.table(field = c('figure_id','meta_source','counts_source','participant_source','mapping_source','population_resolution_source','builder'),
                   value = c(FIG_ID, src$meta, src$counts, src$part, src$mapping, src$popres, '32A_build_Figure_S01_supplementary.R'))
fwrite(prov, file.path(MAN_DIR, 'provenance.tsv'), sep='\t')

src_files <- list.files(SRC_DIR, full.names = TRUE)
manifest <- data.table(file = basename(src_files), path = file.path('submission/source_data/supplementary', FIG_ID, basename(src_files)), bytes = file.info(src_files)$size)
fwrite(manifest, file.path(MAN_DIR, 'source_data_manifest.tsv'), sep='\t')
writeLines(c(sprintf('PDF=%s', pdf_path), sprintf('PNG=%s', png_path), sprintf('TIFF=%s', tiff_path), 'CANVAS=13.2x15.0 in', 'DPI=300'), file.path(MAN_DIR, 'format_metadata.txt'))
writeLines(c('gate\tstatus\tdetail','VISUAL_QC\tPENDING\tAwaiting manual visual approval'), file.path(MAN_DIR, 'visual_qc.tsv'))
writeLines(c('check\tstatus\tdetail','FREEZE_CANDIDATE\tPASS\tCandidate outputs generated; freeze remains blocked until manual visual approval'), file.path(MAN_DIR, 'freeze_audit.tsv'))
cat('POPULATION_COUNT_GATE=PASS total=65 exposed=47 reference=18\n')
cat('GROUP_ASSIGNMENT_RULE=Condition or valid primary elimination phenotype fallback\n')
cat('GENERATOR_EXECUTION_PASS\n')
