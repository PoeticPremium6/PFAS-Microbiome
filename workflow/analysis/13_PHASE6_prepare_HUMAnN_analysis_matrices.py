#!/usr/bin/env python3

import csv
import math
import heapq
from pathlib import Path

REV = Path("/path/to/PFAS_mSystems_revision")
PH6 = REV / "05_humann_functional_profiles"
MERGED = PH6 / "02_merged"
QC = PH6 / "03_qc"
AMAT = PH6 / "04_analysis_matrices"

QC.mkdir(parents=True, exist_ok=True)
AMAT.mkdir(parents=True, exist_ok=True)

SPECIAL = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED"}

def is_unstratified(x):
    return "|" not in x

def parse_vals(vals):
    out = []
    for v in vals:
        try:
            out.append(float(v))
        except Exception:
            out.append(0.0)
    return out

def write_row(writer, feat, vals):
    writer.writerow([feat] + [format(x, ".10g") for x in vals])

def process_matrix(input_path, feature_label, out_prefix, min_prev_fraction=0.10, top_n=500):
    unstrat_out = AMAT / f"{out_prefix}_unstratified_no_unmapped.tsv"
    filt_out = AMAT / f"{out_prefix}_unstratified_no_unmapped_prevalence10.tsv"
    stats_out = AMAT / f"{out_prefix}_unstratified_no_unmapped_feature_stats.tsv"
    top_out = AMAT / f"{out_prefix}_top_variable_features.tsv"

    top_heap = []
    total_rows = 0
    unstrat_rows = 0
    filt_rows = 0

    with open(input_path, newline="") as fh, \
         open(unstrat_out, "w", newline="") as uo, \
         open(filt_out, "w", newline="") as fo, \
         open(stats_out, "w", newline="") as so:

        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        n = len(samples)
        min_prev = max(1, math.ceil(n * min_prev_fraction))

        u_writer = csv.writer(uo, delimiter="\t")
        f_writer = csv.writer(fo, delimiter="\t")
        s_writer = csv.writer(so, delimiter="\t")

        u_writer.writerow([feature_label] + samples)
        f_writer.writerow([feature_label] + samples)
        s_writer.writerow([feature_label, "prevalence_n", "prevalence_fraction", "mean", "variance_log1p", "sum"])

        for row in reader:
            total_rows += 1
            if not row:
                continue
            feat = row[0].strip()
            if not feat:
                continue
            if feat in SPECIAL:
                continue
            if not is_unstratified(feat):
                continue

            vals = parse_vals(row[1:])
            if len(vals) < n:
                vals += [0.0] * (n - len(vals))
            vals = vals[:n]

            unstrat_rows += 1
            write_row(u_writer, feat, vals)

            prev = sum(1 for x in vals if x > 0)
            total = sum(vals)
            logs = [math.log1p(x) for x in vals]
            mean_log = sum(logs) / n
            var_log = sum((x - mean_log) ** 2 for x in logs) / (n - 1) if n > 1 else 0.0

            s_writer.writerow([feat, prev, format(prev / n, ".6g"), format(sum(vals) / n, ".10g"), format(var_log, ".10g"), format(total, ".10g")])

            if prev >= min_prev and total > 0:
                filt_rows += 1
                write_row(f_writer, feat, vals)

                item = (var_log, feat, prev, total)
                if len(top_heap) < top_n:
                    heapq.heappush(top_heap, item)
                else:
                    if item[0] > top_heap[0][0]:
                        heapq.heapreplace(top_heap, item)

    top_rows = sorted(top_heap, reverse=True)
    with open(top_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([feature_label, "variance_log1p", "prevalence_n", "sum"])
        for var_log, feat, prev, total in top_rows:
            writer.writerow([feat, format(var_log, ".10g"), prev, format(total, ".10g")])

    return {
        "input": str(input_path),
        "unstratified_no_unmapped": str(unstrat_out),
        "prevalence10": str(filt_out),
        "feature_stats": str(stats_out),
        "top_variable": str(top_out),
        "total_rows": total_rows,
        "unstratified_no_unmapped_rows": unstrat_rows,
        "prevalence10_rows": filt_rows,
        "samples": n,
        "min_prev": min_prev,
    }

results = []
results.append(process_matrix(
    MERGED / "Ronneby_HUMAnN_pathabundance_cpm.tsv",
    "Pathway",
    "Ronneby_HUMAnN_pathabundance_cpm",
    min_prev_fraction=0.10,
    top_n=500,
))
results.append(process_matrix(
    MERGED / "Ronneby_HUMAnN_genefamilies_cpm.tsv",
    "GeneFamily",
    "Ronneby_HUMAnN_genefamilies_cpm",
    min_prev_fraction=0.10,
    top_n=1000,
))

summary = QC / "PHASE6_HUMAnN_analysis_matrix_preparation_summary.tsv"
with open(summary, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["matrix", "metric", "value"])
    for r in results:
        name = Path(r["input"]).name
        for k, v in r.items():
            if k != "input":
                writer.writerow([name, k, v])

closeout = QC / "PHASE6_ITEM1_RUN_MERGE_QC_COMPLETE.txt"
with open(closeout, "w") as out:
    out.write("Phase 6 item 1 complete: HUMAnN run, merge, renormalization, and basic analysis-matrix preparation complete.\n")
    out.write("Inputs:\n")
    out.write(str(MERGED / "Ronneby_HUMAnN_genefamilies.tsv") + "\n")
    out.write(str(MERGED / "Ronneby_HUMAnN_genefamilies_cpm.tsv") + "\n")
    out.write(str(MERGED / "Ronneby_HUMAnN_pathabundance.tsv") + "\n")
    out.write(str(MERGED / "Ronneby_HUMAnN_pathabundance_cpm.tsv") + "\n")
    out.write(str(MERGED / "Ronneby_HUMAnN_pathcoverage.tsv") + "\n")
    out.write("\nQC summary:\n")
    out.write(str(summary) + "\n")

print("DONE")
print(summary)
print(closeout)
