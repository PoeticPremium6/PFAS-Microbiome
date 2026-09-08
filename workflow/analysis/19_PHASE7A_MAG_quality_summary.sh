#!/usr/bin/env bash

REV="/path/to/PFAS_mSystems_revision"
ORIG="/path/to/hpc_work/PFAS"
PH7="${REV}/06_MAG_resource_integration"

RAW_DIR="${ORIG}/MAG_catalogue/raw_bins"
FILT_DIR="${ORIG}/MAG_catalogue/filtered_80comp_10contam"
CHECKM_ALL="${ORIG}/MAG_catalogue/MAG_quality_checkm_all_chunks.tsv"
CHECKM_FILT="${ORIG}/MAG_catalogue/MAG_quality_80comp_10contam.tsv"

OUT="${PH7}/02_quality"
FIG="${PH7}/07_figures"
QC="${PH7}/08_qc"
TAB="${REV}/12_tables"

mkdir -p "$OUT" "$FIG" "$QC" "$TAB"

echo "=== Phase 7A item 1: MAG quality summary ==="

for F in "$RAW_DIR" "$FILT_DIR" "$CHECKM_ALL" "$CHECKM_FILT"; do
  if [ ! -e "$F" ]; then
    echo "ERROR missing: $F"
    exit 2
  fi
done

RAW_N=$(find "$RAW_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fna.gz" -o -name "*.fasta.gz" \) | wc -l)
FILT_N=$(find "$FILT_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fna.gz" -o -name "*.fasta.gz" \) | wc -l)

python3 - <<'PY'
from pathlib import Path
import csv
import statistics
import re

REV = Path("/path/to/PFAS_mSystems_revision")
ORIG = Path("/path/to/hpc_work/PFAS")
PH7 = REV / "06_MAG_resource_integration"

RAW_DIR = ORIG / "MAG_catalogue" / "raw_bins"
FILT_DIR = ORIG / "MAG_catalogue" / "filtered_80comp_10contam"
CHECKM_ALL = ORIG / "MAG_catalogue" / "MAG_quality_checkm_all_chunks.tsv"
CHECKM_FILT = ORIG / "MAG_catalogue" / "MAG_quality_80comp_10contam.tsv"

OUT = PH7 / "02_quality"
FIG = PH7 / "07_figures"
QC = PH7 / "08_qc"
TAB = REV / "12_tables"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

def clean(x):
    if x is None:
        return ""
    return str(x).strip().strip('"').strip()

def read_checkm(path):
    with open(path, newline="", errors="ignore") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = [clean(x) for x in next(reader)]
        rows = []
        for row in reader:
            if not row:
                continue
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            d = {header[i]: clean(row[i]) for i in range(len(header))}
            rows.append(d)
    return header, rows

def as_float(x):
    try:
        return float(str(x).replace(",", ""))
    except Exception:
        return None

def as_int(x):
    y = as_float(x)
    if y is None:
        return None
    return int(round(y))

def strip_fasta_ext(name):
    for ext in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name

def fasta_files(d):
    files = []
    for p in d.iterdir():
        if p.is_file() and re.search(r"\.(fa|fna|fasta)(\.gz)?$", p.name):
            files.append(p)
    return sorted(files)

all_header, all_rows = read_checkm(CHECKM_ALL)
filt_header, filt_rows = read_checkm(CHECKM_FILT)

raw_fastas = fasta_files(RAW_DIR)
filt_fastas = fasta_files(FILT_DIR)

filt_map = {strip_fasta_ext(p.name): str(p) for p in filt_fastas}

table_s23 = TAB / "Table_S23_MAG_quality_summary.tsv"
source_s16 = FIG / "Figure_S16_MAG_quality_source_data.tsv"
count_audit = QC / "PHASE7A_MAG_catalogue_count_audit.tsv"
summary_stats = QC / "PHASE7A_MAG_quality_summary_stats.tsv"

selected = [
    ("MAG_ID", "Bin Id"),
    ("marker_lineage", "Marker lineage"),
    ("completeness", "Completeness"),
    ("contamination", "Contamination"),
    ("strain_heterogeneity", "Strain heterogeneity"),
    ("genome_size_bp", "Genome size (bp)"),
    ("ambiguous_bases", "# ambiguous bases"),
    ("scaffolds", "# scaffolds"),
    ("contigs", "# contigs"),
    ("n50_scaffolds", "N50 (scaffolds)"),
    ("n50_contigs", "N50 (contigs)"),
    ("mean_scaffold_length_bp", "Mean scaffold length (bp)"),
    ("mean_contig_length_bp", "Mean contig length (bp)"),
    ("longest_scaffold_bp", "Longest scaffold (bp)"),
    ("longest_contig_bp", "Longest contig (bp)"),
    ("gc", "GC"),
    ("coding_density", "Coding density"),
    ("predicted_genes", "# predicted genes"),
]

numeric_out = {
    "completeness", "contamination", "strain_heterogeneity", "genome_size_bp",
    "ambiguous_bases", "scaffolds", "contigs", "n50_scaffolds", "n50_contigs",
    "mean_scaffold_length_bp", "mean_contig_length_bp", "longest_scaffold_bp",
    "longest_contig_bp", "gc", "coding_density", "predicted_genes"
}

with open(table_s23, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([x[0] for x in selected] + ["filtered_MAG_fasta_path", "in_filtered_MAG_directory"])
    for r in filt_rows:
        mag = clean(r.get("Bin Id", ""))
        row = []
        for outname, inname in selected:
            val = clean(r.get(inname, ""))
            if outname in numeric_out:
                num = as_float(val)
                val = "" if num is None else f"{num:.10g}"
            row.append(val)
        fasta = filt_map.get(mag, "")
        row.append(fasta)
        row.append("TRUE" if fasta else "FALSE")
        writer.writerow(row)

# Copy source data for Figure S16.
with open(table_s23, newline="") as inp, open(source_s16, "w", newline="") as out:
    reader = csv.DictReader(inp, delimiter="\t")
    fields = ["MAG_ID", "completeness", "contamination", "genome_size_bp", "contigs", "n50_contigs", "gc", "predicted_genes"]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in reader:
        writer.writerow({k: r.get(k, "") for k in fields})

def get_col(rows, col):
    vals = []
    for r in rows:
        x = as_float(r.get(col, ""))
        if x is not None:
            vals.append(x)
    return vals

def q(vals, p):
    if not vals:
        return ""
    vals = sorted(vals)
    idx = (len(vals)-1) * p
    lo = int(idx)
    hi = min(lo + 1, len(vals)-1)
    frac = idx - lo
    return vals[lo] * (1-frac) + vals[hi] * frac

metrics = {
    "Completeness": "Completeness",
    "Contamination": "Contamination",
    "Genome_size_bp": "Genome size (bp)",
    "Contigs": "# contigs",
    "N50_contigs": "N50 (contigs)",
    "GC": "GC",
    "Predicted_genes": "# predicted genes",
}

with open(summary_stats, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "n", "mean", "median", "q25", "q75", "min", "max"])
    for label, col in metrics.items():
        vals = get_col(filt_rows, col)
        writer.writerow([
            label,
            len(vals),
            f"{statistics.mean(vals):.6g}" if vals else "",
            f"{statistics.median(vals):.6g}" if vals else "",
            f"{q(vals, 0.25):.6g}" if vals else "",
            f"{q(vals, 0.75):.6g}" if vals else "",
            f"{min(vals):.6g}" if vals else "",
            f"{max(vals):.6g}" if vals else "",
        ])

missing_fasta = []
for r in filt_rows:
    mag = clean(r.get("Bin Id", ""))
    if mag not in filt_map:
        missing_fasta.append(mag)

with open(count_audit, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["raw_bins_fasta_count", len(raw_fastas)])
    writer.writerow(["filtered_MAG_fasta_count", len(filt_fastas)])
    writer.writerow(["checkm_all_rows", len(all_rows)])
    writer.writerow(["checkm_filtered_rows", len(filt_rows)])
    writer.writerow(["expected_raw_bins", 5925])
    writer.writerow(["expected_filtered_MAGs", 1973])
    writer.writerow(["raw_bins_count_matches_expected", str(len(raw_fastas) == 5925)])
    writer.writerow(["filtered_MAG_count_matches_expected", str(len(filt_fastas) == 1973)])
    writer.writerow(["checkm_filtered_rows_match_filtered_FASTA_count", str(len(filt_rows) == len(filt_fastas))])
    writer.writerow(["filtered_quality_rows_missing_FASTA", len(missing_fasta)])
    writer.writerow(["raw_bins_dir", str(RAW_DIR)])
    writer.writerow(["filtered_MAG_dir", str(FILT_DIR)])
    writer.writerow(["checkm_all_input", str(CHECKM_ALL)])
    writer.writerow(["checkm_filtered_input", str(CHECKM_FILT)])
    writer.writerow(["Table_S23", str(table_s23)])
    writer.writerow(["Figure_S16_source_data", str(source_s16)])

print("DONE_PY")
print(count_audit)
print(summary_stats)
print(table_s23)
print(source_s16)
PY

cat > "${REV}/scripts/19_PHASE7A_MAG_quality_figure.R" <<'RS'
rev <- "/path/to/PFAS_mSystems_revision"
ph7 <- file.path(rev, "06_MAG_resource_integration")
figdir <- file.path(ph7, "07_figures")
src <- file.path(figdir, "Figure_S16_MAG_quality_source_data.tsv")
outpdf <- file.path(figdir, "Figure_S16_MAG_quality_distributions.pdf")

x <- read.table(src, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)

purple_dark <- "#3B0F70"
purple_main <- "#6A00A8"
purple_light <- "#B12A90"
grey_dark <- "#333333"

pdf(outpdf, width=9, height=7)
par(mfrow=c(2,3), mar=c(4,4,2,1), col.axis=grey_dark, col.lab=grey_dark, fg=grey_dark)

hist(x$completeness, breaks=30, col=purple_main, border="white",
     main="MAG completeness", xlab="Completeness (%)")
abline(v=80, col=purple_dark, lwd=2, lty=2)

hist(x$contamination, breaks=30, col=purple_light, border="white",
     main="MAG contamination", xlab="Contamination (%)")
abline(v=10, col=purple_dark, lwd=2, lty=2)

hist(x$genome_size_bp / 1e6, breaks=30, col=purple_main, border="white",
     main="Genome size", xlab="Genome size (Mbp)")

hist(x$contigs, breaks=30, col=purple_light, border="white",
     main="Contig count", xlab="Contigs")

hist(x$n50_contigs / 1000, breaks=30, col=purple_main, border="white",
     main="Contig N50", xlab="N50 (kbp)")

hist(x$gc, breaks=30, col=purple_light, border="white",
     main="GC content", xlab="GC")
dev.off()

cat("DONE_R\n")
cat(outpdf, "\n")
RS

module purge || true
module load gcc12-env/12.1.0 || true
module load micromamba/1.4.2 || true
export MAMBA_ROOT_PREFIX=/path/to/hpc_work/micromamba
eval "$(micromamba shell hook --shell bash)"

if micromamba run -n PFAS_PUBLICMETA Rscript --version >/dev/null 2>&1; then
  micromamba run -n PFAS_PUBLICMETA Rscript "${REV}/scripts/19_PHASE7A_MAG_quality_figure.R"
else
  echo "WARNING: Rscript unavailable; Figure S16 PDF not generated."
fi

echo
echo "=== count audit ==="
cat "${QC}/PHASE7A_MAG_catalogue_count_audit.tsv"

echo
echo "=== quality summary stats ==="
cat "${QC}/PHASE7A_MAG_quality_summary_stats.tsv"

echo
echo "=== outputs ==="
ls -lh \
  "${TAB}/Table_S23_MAG_quality_summary.tsv" \
  "${FIG}/Figure_S16_MAG_quality_source_data.tsv" \
  "${FIG}/Figure_S16_MAG_quality_distributions.pdf" 2>/dev/null || true
