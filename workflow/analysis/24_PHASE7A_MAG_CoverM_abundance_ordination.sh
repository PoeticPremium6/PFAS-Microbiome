#!/usr/bin/env bash

REV="/path/to/PFAS_mSystems_revision"
ORIG="/path/to/hpc_work/PFAS"
PH7="${REV}/06_MAG_resource_integration"

COVERM_DIR="${ORIG}/MAG_catalogue/coverm_MAG_abundance"
OUT="${PH7}/06_coverm_abundance"
FIG="${PH7}/07_figures"
QC="${PH7}/08_qc"
TAB="${REV}/12_tables"

mkdir -p "$OUT" "$FIG" "$QC" "$TAB"

echo "=== Phase 7A item 5: CoverM MAG abundance integration ==="

if [ ! -d "$COVERM_DIR" ]; then
  echo "ERROR missing CoverM dir: $COVERM_DIR"
  exit 2
fi

python3 - <<'PY'
from pathlib import Path
import csv
import re
import math
import statistics
from collections import defaultdict

REV = Path("/path/to/PFAS_mSystems_revision")
ORIG = Path("/path/to/hpc_work/PFAS")
PH7 = REV / "06_MAG_resource_integration"

COVERM_DIR = ORIG / "MAG_catalogue" / "coverm_MAG_abundance"
OUT = PH7 / "06_coverm_abundance"
FIG = PH7 / "07_figures"
QC = PH7 / "08_qc"
TAB = REV / "12_tables"

TABLE_S23 = TAB / "Table_S23_MAG_quality_summary_with_reheadered_paths.tsv"
if not TABLE_S23.exists():
    TABLE_S23 = TAB / "Table_S23_MAG_quality_summary.tsv"

TABLE_S24 = TAB / "Table_S24_MAG_GTDBTk_taxonomy_fastANI95_representatives.tsv"
TABLE_S24B = TAB / "Table_S24b_MAG_fastANI_95ANI_dereplication.tsv"
META = REV / "01_metadata_freeze" / "sample_metadata_revision_master.tsv"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

def clean_sample(x):
    m = re.search(r"[A-H][0-9]{2}", str(x))
    return m.group(0) if m else str(x)

def clean_mag(x):
    x = str(x or "").strip().strip('"')
    x = Path(x).name
    for ext in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
        if x.endswith(ext):
            x = x[:-len(ext)]
    return x

def fnum(x):
    try:
        if x is None or str(x).strip() == "":
            return 0.0
        return float(str(x).replace(",", ""))
    except Exception:
        return 0.0

def median(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return 0.0
    return statistics.median(vals)

# Quality table.
quality = {}
mag_order = []
if TABLE_S23.exists():
    with open(TABLE_S23, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            mag = clean_mag(r.get("MAG_ID", ""))
            if not mag:
                continue
            quality[mag] = r
            mag_order.append(mag)

# Derep mapping.
derep = {}
if TABLE_S24B.exists():
    with open(TABLE_S24B, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            mag = clean_mag(r.get("MAG_ID", ""))
            if mag:
                derep[mag] = r

# Representative taxonomy.
tax_by_rep = {}
if TABLE_S24.exists():
    with open(TABLE_S24, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            mag = clean_mag(r.get("MAG_ID", ""))
            if mag:
                tax_by_rep[mag] = r

# Metadata.
metadata = {}
if META.exists():
    with open(META, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sid = r.get("SampleID", "")
            if sid:
                metadata[sid] = r

coverm_files = sorted(COVERM_DIR.glob("*_MAG_ID_abundance.tsv"))

rel = defaultdict(dict)
mean_cov = defaultdict(dict)
covered_frac = defaultdict(dict)
samples = []

parse_warnings = []

for path in coverm_files:
    sample = clean_sample(path.name)
    samples.append(sample)

    with open(path, newline="", errors="ignore") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []

        genome_col = "Genome" if "Genome" in header else header[0]
        rel_cols = [h for h in header if "Relative Abundance" in h]
        mean_cols = [h for h in header if h.endswith(" Mean") or " Mean" in h]
        frac_cols = [h for h in header if "Covered Fraction" in h]

        if not rel_cols or not mean_cols or not frac_cols:
            parse_warnings.append(f"{path}\tHEADER_PARSE_WARNING\t{header}")
            continue

        rel_col = rel_cols[0]
        mean_col = mean_cols[0]
        frac_col = frac_cols[0]

        for r in reader:
            mag = clean_mag(r.get(genome_col, ""))
            if not mag:
                continue
            rel[mag][sample] = fnum(r.get(rel_col, 0))
            mean_cov[mag][sample] = fnum(r.get(mean_col, 0))
            covered_frac[mag][sample] = fnum(r.get(frac_col, 0))

samples = sorted(set(samples))

all_mags = []
seen = set()
for m in mag_order:
    if m not in seen:
        all_mags.append(m)
        seen.add(m)
for m in sorted(rel):
    if m not in seen:
        all_mags.append(m)
        seen.add(m)

def write_matrix(path, data):
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["MAG_ID"] + samples)
        for mag in all_mags:
            writer.writerow([mag] + [f"{data.get(mag, {}).get(s, 0.0):.10g}" for s in samples])

rel_matrix = OUT / "Ronneby_MAG_relative_abundance_percent.tsv"
mean_matrix = OUT / "Ronneby_MAG_mean_coverage.tsv"
frac_matrix = OUT / "Ronneby_MAG_covered_fraction.tsv"

write_matrix(rel_matrix, rel)
write_matrix(mean_matrix, mean_cov)
write_matrix(frac_matrix, covered_frac)

# Feature stats.
feature_stats = OUT / "Ronneby_MAG_abundance_feature_stats.tsv"
min_prev = max(1, math.ceil(0.10 * len(samples)))

feature_rows = []
for mag in all_mags:
    vals = [rel.get(mag, {}).get(s, 0.0) for s in samples]
    frac_vals = [covered_frac.get(mag, {}).get(s, 0.0) for s in samples]
    nonzero = [v for v in vals if v > 0]
    prev = len(nonzero)
    row = {
        "MAG_ID": mag,
        "samples_with_relative_abundance_gt0": prev,
        "prevalence_fraction": prev / len(samples) if samples else 0,
        "mean_relative_abundance_percent": statistics.mean(vals) if vals else 0,
        "median_relative_abundance_percent": statistics.median(vals) if vals else 0,
        "median_nonzero_relative_abundance_percent": median(nonzero),
        "max_relative_abundance_percent": max(vals) if vals else 0,
        "total_relative_abundance_percent_across_samples": sum(vals),
        "samples_with_covered_fraction_gt0": sum(1 for v in frac_vals if v > 0),
        "mean_covered_fraction": statistics.mean(frac_vals) if frac_vals else 0,
    }
    feature_rows.append(row)

with open(feature_stats, "w", newline="") as out:
    fields = [
        "MAG_ID", "samples_with_relative_abundance_gt0", "prevalence_fraction",
        "mean_relative_abundance_percent", "median_relative_abundance_percent",
        "median_nonzero_relative_abundance_percent", "max_relative_abundance_percent",
        "total_relative_abundance_percent_across_samples",
        "samples_with_covered_fraction_gt0", "mean_covered_fraction"
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in feature_rows:
        writer.writerow({k: (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k]) for k in fields})

# Prevalence-filtered matrix.
prev_mags = [r["MAG_ID"] for r in feature_rows if r["samples_with_relative_abundance_gt0"] >= min_prev]
prev_matrix = OUT / "Ronneby_MAG_relative_abundance_percent_prevalence10.tsv"
with open(prev_matrix, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["MAG_ID"] + samples)
    for mag in prev_mags:
        writer.writerow([mag] + [f"{rel.get(mag, {}).get(s, 0.0):.10g}" for s in samples])

# Sample summary metrics.
sample_summary = OUT / "Ronneby_MAG_sample_summary_metrics.tsv"
sample_rows = []

for s in samples:
    vals = [rel.get(mag, {}).get(s, 0.0) for mag in all_mags]
    frac_vals = [covered_frac.get(mag, {}).get(s, 0.0) for mag in all_mags]
    cov_vals = [mean_cov.get(mag, {}).get(s, 0.0) for mag in all_mags]

    total = sum(vals)
    p = [v / total for v in vals if total > 0 and v > 0]
    shannon = -sum(x * math.log(x) for x in p) if p else 0
    simpson = 1 - sum(x * x for x in p) if p else 0
    evenness = shannon / math.log(len(p)) if len(p) > 1 else 0
    sorted_vals = sorted(vals, reverse=True)

    m = metadata.get(s, {})

    sample_rows.append({
        "SampleID": s,
        "n_MAGs_relative_abundance_gt0": sum(1 for v in vals if v > 0),
        "n_MAGs_covered_fraction_gt0": sum(1 for v in frac_vals if v > 0),
        "sum_relative_abundance_percent": total,
        "shannon_MAG_abundance": shannon,
        "simpson_MAG_abundance": simpson,
        "pielou_evenness_MAG_abundance": evenness,
        "top1_relative_abundance_percent": sorted_vals[0] if sorted_vals else 0,
        "top5_cumulative_relative_abundance_percent": sum(sorted_vals[:5]),
        "mean_MAG_coverage_all_MAGs": statistics.mean(cov_vals) if cov_vals else 0,
        "median_MAG_coverage_all_MAGs": statistics.median(cov_vals) if cov_vals else 0,
        "ExposureGroup": m.get("ExposureGroup", ""),
        "Age": m.get("Age", ""),
        "Sex": m.get("Sex", ""),
        "BMI": m.get("BMI", ""),
        "SequencingBatch": m.get("SequencingBatch", ""),
    })

with open(sample_summary, "w", newline="") as out:
    fields = [
        "SampleID", "n_MAGs_relative_abundance_gt0", "n_MAGs_covered_fraction_gt0",
        "sum_relative_abundance_percent", "shannon_MAG_abundance", "simpson_MAG_abundance",
        "pielou_evenness_MAG_abundance", "top1_relative_abundance_percent",
        "top5_cumulative_relative_abundance_percent", "mean_MAG_coverage_all_MAGs",
        "median_MAG_coverage_all_MAGs", "ExposureGroup", "Age", "Sex", "BMI", "SequencingBatch"
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in sample_rows:
        outrow = {}
        for k in fields:
            v = r[k]
            outrow[k] = f"{v:.10g}" if isinstance(v, float) else v
        writer.writerow(outrow)

# Table S25: MAG abundance + quality + taxonomy + derep summary.
table_s25 = TAB / "Table_S25_MAG_abundance_taxonomy_quality_summary.tsv"

stats_by_mag = {r["MAG_ID"]: r for r in feature_rows}

fields = [
    "MAG_ID",
    "cluster_id", "representative_MAG_ID", "is_representative",
    "GTDB_classification", "GTDB_domain", "GTDB_phylum", "GTDB_class", "GTDB_order",
    "GTDB_family", "GTDB_genus", "GTDB_species",
    "completeness", "contamination", "genome_size_bp", "contigs", "n50_contigs", "gc", "predicted_genes",
    "samples_with_relative_abundance_gt0", "prevalence_fraction",
    "mean_relative_abundance_percent", "median_relative_abundance_percent",
    "median_nonzero_relative_abundance_percent", "max_relative_abundance_percent",
    "total_relative_abundance_percent_across_samples",
    "samples_with_covered_fraction_gt0", "mean_covered_fraction"
]

with open(table_s25, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for mag in all_mags:
        q = quality.get(mag, {})
        d = derep.get(mag, {})
        rep = clean_mag(d.get("representative_MAG_ID", mag if mag in tax_by_rep else ""))
        tax = tax_by_rep.get(rep, {})
        st = stats_by_mag.get(mag, {})

        row = {
            "MAG_ID": mag,
            "cluster_id": d.get("cluster_id", ""),
            "representative_MAG_ID": rep,
            "is_representative": d.get("is_representative", str(mag == rep) if rep else ""),
            "GTDB_classification": tax.get("classification", ""),
            "GTDB_domain": tax.get("domain", ""),
            "GTDB_phylum": tax.get("phylum", ""),
            "GTDB_class": tax.get("class", ""),
            "GTDB_order": tax.get("order", ""),
            "GTDB_family": tax.get("family", ""),
            "GTDB_genus": tax.get("genus", ""),
            "GTDB_species": tax.get("species", ""),
            "completeness": q.get("completeness", ""),
            "contamination": q.get("contamination", ""),
            "genome_size_bp": q.get("genome_size_bp", ""),
            "contigs": q.get("contigs", ""),
            "n50_contigs": q.get("n50_contigs", ""),
            "gc": q.get("gc", ""),
            "predicted_genes": q.get("predicted_genes", ""),
        }

        for k in [
            "samples_with_relative_abundance_gt0", "prevalence_fraction",
            "mean_relative_abundance_percent", "median_relative_abundance_percent",
            "median_nonzero_relative_abundance_percent", "max_relative_abundance_percent",
            "total_relative_abundance_percent_across_samples",
            "samples_with_covered_fraction_gt0", "mean_covered_fraction"
        ]:
            v = st.get(k, "")
            row[k] = f"{v:.10g}" if isinstance(v, float) else v

        writer.writerow(row)

# Top abundant source data.
top_abundant = FIG / "Figure_S17C_MAG_top20_abundance_source_data.tsv"
top20 = sorted(feature_rows, key=lambda r: r["mean_relative_abundance_percent"], reverse=True)[:20]
top20_mags = [r["MAG_ID"] for r in top20]

with open(top_abundant, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["MAG_ID", "SampleID", "relative_abundance_percent"])
    for mag in top20_mags:
        for s in samples:
            writer.writerow([mag, s, f"{rel.get(mag, {}).get(s, 0.0):.10g}"])

# Audit.
audit = QC / "PHASE7A_MAG_CoverM_abundance_integration_audit.tsv"
with open(audit, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["status", "COMPLETE"])
    writer.writerow(["coverm_dir", str(COVERM_DIR)])
    writer.writerow(["coverm_files", len(coverm_files)])
    writer.writerow(["samples", len(samples)])
    writer.writerow(["MAGs_in_Table_S23", len(mag_order)])
    writer.writerow(["MAGs_seen_in_CoverM", len(rel)])
    writer.writerow(["MAGs_in_output_matrix", len(all_mags)])
    writer.writerow(["minimum_prevalence_for_prevalence10", min_prev])
    writer.writerow(["MAGs_prevalence10", len(prev_mags)])
    writer.writerow(["parse_warnings", len(parse_warnings)])
    writer.writerow(["relative_abundance_matrix", str(rel_matrix)])
    writer.writerow(["mean_coverage_matrix", str(mean_matrix)])
    writer.writerow(["covered_fraction_matrix", str(frac_matrix)])
    writer.writerow(["prevalence10_matrix", str(prev_matrix)])
    writer.writerow(["feature_stats", str(feature_stats)])
    writer.writerow(["sample_summary", str(sample_summary)])
    writer.writerow(["Table_S25", str(table_s25)])
    writer.writerow(["Figure_S17C_top20_source_data", str(top_abundant)])

if parse_warnings:
    with open(QC / "PHASE7A_MAG_CoverM_parse_warnings.tsv", "w") as out:
        out.write("\n".join(parse_warnings) + "\n")

print(audit)
print(table_s25)
print(sample_summary)
print(prev_matrix)
PY

cat > "${REV}/scripts/24_PHASE7A_MAG_CoverM_figures.R" <<'RS'
rev <- "/path/to/PFAS_mSystems_revision"
ph7 <- file.path(rev, "06_MAG_resource_integration")

outdir <- file.path(ph7, "06_coverm_abundance")
figdir <- file.path(ph7, "07_figures")
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)

matfile <- file.path(outdir, "Ronneby_MAG_relative_abundance_percent_prevalence10.tsv")
samplefile <- file.path(outdir, "Ronneby_MAG_sample_summary_metrics.tsv")
topfile <- file.path(figdir, "Figure_S17C_MAG_top20_abundance_source_data.tsv")

pcoa_pdf <- file.path(figdir, "Figure_6D_MAG_abundance_Bray_PCoA.pdf")
summary_pdf <- file.path(figdir, "Figure_6E_MAG_abundance_summary_metrics.pdf")
heat_pdf <- file.path(figdir, "Figure_S17C_MAG_top20_abundance_heatmap.pdf")

score_out <- file.path(outdir, "Ronneby_MAG_abundance_Bray_PCoA_sample_scores.tsv")
ord_summary_out <- file.path(outdir, "Ronneby_MAG_abundance_ordination_summary.tsv")

purple_dark <- "#3B0F70"
purple_main <- "#6A00A8"
purple_mid <- "#8C2981"
purple_light <- "#B12A90"
lavender <- "#D7B5F7"
grey_dark <- "#333333"
grey_mid <- "#777777"

x <- read.table(matfile, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)
rownames(x) <- x[[1]]
mat <- as.matrix(x[, -1, drop=FALSE])
storage.mode(mat) <- "numeric"
mat <- t(mat)

sample_summary <- read.table(samplefile, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)
rownames(sample_summary) <- sample_summary$SampleID
sample_summary <- sample_summary[rownames(mat), , drop=FALSE]

bray <- function(m) {
  n <- nrow(m)
  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      denom <- sum(m[i, ] + m[j, ])
      d[i, j] <- ifelse(denom > 0, sum(abs(m[i, ] - m[j, ])) / denom, 0)
    }
  }
  rownames(d) <- rownames(m)
  colnames(d) <- rownames(m)
  as.dist(d)
}

bc <- bray(mat)
pc <- cmdscale(bc, k=2, eig=TRUE)
scores <- as.data.frame(pc$points)
colnames(scores) <- c("PCoA1", "PCoA2")
scores$SampleID <- rownames(scores)

eig <- pc$eig
pos <- eig[eig > 0]
p1 <- ifelse(length(pos) >= 1, pos[1] / sum(pos), NA)
p2 <- ifelse(length(pos) >= 2, pos[2] / sum(pos), NA)

scores <- merge(scores, sample_summary, by="SampleID", all.x=TRUE, sort=FALSE)

write.table(scores, score_out, sep="\t", quote=FALSE, row.names=FALSE)

write.table(
  data.frame(
    metric=c("n_samples", "n_MAGs_prevalence10", "Bray_PCoA1_positive_eigen_fraction", "Bray_PCoA2_positive_eigen_fraction"),
    value=c(nrow(mat), ncol(mat), p1, p2)
  ),
  ord_summary_out,
  sep="\t", quote=FALSE, row.names=FALSE
)

grp <- scores$ExposureGroup
if (all(is.na(grp)) || length(unique(na.omit(grp))) < 2) {
  grp <- rep("Ronneby", nrow(scores))
}
grp <- as.factor(grp)
pal <- c(lavender, purple_light, purple_main, purple_dark, grey_mid)
cols <- pal[((as.integer(grp) - 1) %% length(pal)) + 1]

pdf(pcoa_pdf, width=6.5, height=5.5)
par(mar=c(4.5,4.5,2,1), fg=grey_dark, col.axis=grey_dark, col.lab=grey_dark)
plot(
  scores$PCoA1, scores$PCoA2,
  pch=21, bg=cols, col=grey_dark, cex=1.3,
  xlab=paste0("PCoA1 (", round(100*p1, 1), "%)"),
  ylab=paste0("PCoA2 (", round(100*p2, 1), "%)"),
  main="MAG abundance Bray-Curtis PCoA"
)
legend("topright", legend=levels(grp), pt.bg=pal[seq_along(levels(grp))], pch=21, bty="n", cex=0.8)
dev.off()

pdf(summary_pdf, width=8, height=5)
par(mfrow=c(1,2), mar=c(5,4,2,1), fg=grey_dark, col.axis=grey_dark, col.lab=grey_dark)

if (length(unique(grp)) >= 2) {
  boxplot(
    scores$n_MAGs_relative_abundance_gt0 ~ grp,
    col=pal[seq_along(levels(grp))],
    border=grey_dark,
    las=2,
    ylab="MAGs detected",
    main="Detected MAG richness"
  )
  boxplot(
    scores$shannon_MAG_abundance ~ grp,
    col=pal[seq_along(levels(grp))],
    border=grey_dark,
    las=2,
    ylab="Shannon diversity",
    main="MAG abundance diversity"
  )
} else {
  hist(
    scores$n_MAGs_relative_abundance_gt0,
    breaks=20,
    col=purple_main,
    border="white",
    xlab="MAGs detected",
    main="Detected MAG richness"
  )
  hist(
    scores$shannon_MAG_abundance,
    breaks=20,
    col=purple_light,
    border="white",
    xlab="Shannon diversity",
    main="MAG abundance diversity"
  )
}
dev.off()

top <- read.table(topfile, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)
wide <- reshape(top, idvar="MAG_ID", timevar="SampleID", direction="wide")
rownames(wide) <- wide$MAG_ID
hm <- as.matrix(wide[, -1, drop=FALSE])
storage.mode(hm) <- "numeric"
colnames(hm) <- sub("^relative_abundance_percent\\.", "", colnames(hm))
hm <- log10(hm + 1)
hm <- hm[, order(colnames(hm)), drop=FALSE]

pdf(heat_pdf, width=9, height=6)
par(mar=c(7,10,2,2), fg=grey_dark, col.axis=grey_dark, col.lab=grey_dark)
image(
  x=seq_len(ncol(hm)),
  y=seq_len(nrow(hm)),
  z=t(hm[nrow(hm):1, ]),
  col=colorRampPalette(c("white", lavender, purple_light, purple_dark))(100),
  axes=FALSE,
  xlab="Sample",
  ylab="MAG",
  main="Top 20 MAGs by mean relative abundance"
)
axis(1, at=seq_len(ncol(hm)), labels=colnames(hm), las=2, cex.axis=0.5)
axis(2, at=seq_len(nrow(hm)), labels=rev(rownames(hm)), las=2, cex.axis=0.5)
box()
dev.off()

cat("DONE\n")
cat(pcoa_pdf, "\n")
cat(summary_pdf, "\n")
cat(heat_pdf, "\n")
RS

module purge || true
module load gcc12-env/12.1.0 || true
module load micromamba/1.4.2 || true
export MAMBA_ROOT_PREFIX=/path/to/hpc_work/micromamba
eval "$(micromamba shell hook --shell bash)"

micromamba run -n PFAS_PUBLICMETA Rscript "${REV}/scripts/24_PHASE7A_MAG_CoverM_figures.R"

echo
echo "=== CoverM abundance integration audit ==="
cat "${QC}/PHASE7A_MAG_CoverM_abundance_integration_audit.tsv"

echo
echo "=== Table S25 preview ==="
ls -lh "${TAB}/Table_S25_MAG_abundance_taxonomy_quality_summary.tsv"
head -5 "${TAB}/Table_S25_MAG_abundance_taxonomy_quality_summary.tsv"

echo
echo "=== sample summary preview ==="
head -10 "${OUT}/Ronneby_MAG_sample_summary_metrics.tsv"

echo
echo "=== ordination summary ==="
cat "${OUT}/Ronneby_MAG_abundance_ordination_summary.tsv"

echo
echo "=== figure outputs ==="
ls -lh \
  "${FIG}/Figure_6D_MAG_abundance_Bray_PCoA.pdf" \
  "${FIG}/Figure_6E_MAG_abundance_summary_metrics.pdf" \
  "${FIG}/Figure_S17C_MAG_top20_abundance_heatmap.pdf"
