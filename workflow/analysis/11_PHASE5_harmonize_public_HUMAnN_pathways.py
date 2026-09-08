#!/usr/bin/env python3

import csv
from pathlib import Path

REV = Path("/path/to/PFAS_mSystems_revision")
PH5 = REV / "04_public_dataset_harmonization"
PH6 = REV / "05_humann_functional_profiles"

ronneby_file = PH6 / "02_merged" / "Ronneby_HUMAnN_pathabundance_cpm.tsv"
nielsen_file = PH5 / "05_humann_pathways" / "Nielsen_pathway_abundance_matrix_only.tsv"
hmp_file = PH5 / "05_humann_pathways" / "HMP_pathway_abundance_matrix_only.tsv"

outdir = PH5 / "05_humann_pathways"
qcdir = PH5 / "06_qc"
outdir.mkdir(parents=True, exist_ok=True)
qcdir.mkdir(parents=True, exist_ok=True)

special_rows = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED"}

def clean_feature(x):
    x = x.strip()
    if not x:
        return ""
    return x

def is_unstratified(feature):
    return "|" not in feature

def read_matrix(path, cohort_prefix):
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        sample_names = [f"{cohort_prefix}_{s.strip()}" for s in header[1:]]
        data = {}
        for row in reader:
            if not row:
                continue
            feat = clean_feature(row[0])
            if not feat:
                continue
            if not is_unstratified(feat):
                continue
            vals = row[1:]
            if len(vals) < len(sample_names):
                vals = vals + ["0"] * (len(sample_names) - len(vals))
            vals = vals[:len(sample_names)]
            if feat not in data:
                data[feat] = [0.0] * len(sample_names)
            for i, v in enumerate(vals):
                try:
                    data[feat][i] += float(v)
                except Exception:
                    pass
    return sample_names, data

def write_matrix(path, features, sample_blocks):
    samples = []
    for _, sample_names, _ in sample_blocks:
        samples.extend(sample_names)
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Pathway"] + samples)
        for feat in features:
            row = [feat]
            for _, sample_names, data in sample_blocks:
                vals = data.get(feat)
                if vals is None:
                    row.extend(["0"] * len(sample_names))
                else:
                    row.extend([format(x, ".10g") for x in vals])
            writer.writerow(row)

print("Reading Ronneby...")
ron_samples, ron = read_matrix(ronneby_file, "Ronneby")
print("Reading Nielsen...")
nie_samples, nie = read_matrix(nielsen_file, "Nielsen")
print("Reading HMP...")
hmp_samples, hmp = read_matrix(hmp_file, "HMP")

sets = {
    "Ronneby": set(ron),
    "Nielsen": set(nie),
    "HMP": set(hmp),
}
union = sorted(set.union(*sets.values()))
shared = sorted(set.intersection(*sets.values()))
union_no_special = [x for x in union if x not in special_rows]
shared_no_special = [x for x in shared if x not in special_rows]

blocks = [
    ("Ronneby", ron_samples, ron),
    ("Nielsen", nie_samples, nie),
    ("HMP", hmp_samples, hmp),
]

outputs = {
    "Ronneby_Nielsen_HMP_HUMAnN_pathways_unstratified_union.tsv": union,
    "Ronneby_Nielsen_HMP_HUMAnN_pathways_unstratified_shared.tsv": shared,
    "Ronneby_Nielsen_HMP_HUMAnN_pathways_unstratified_union_no_unmapped.tsv": union_no_special,
    "Ronneby_Nielsen_HMP_HUMAnN_pathways_unstratified_shared_no_unmapped.tsv": shared_no_special,
}

for name, feats in outputs.items():
    outpath = outdir / name
    print(f"Writing {outpath} ({len(feats)} features)")
    write_matrix(outpath, feats, blocks)

summary = qcdir / "PHASE5_HUMAnN_pathway_harmonization_summary.tsv"
with open(summary, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["Ronneby_samples", len(ron_samples)])
    writer.writerow(["Nielsen_samples", len(nie_samples)])
    writer.writerow(["HMP_samples", len(hmp_samples)])
    writer.writerow(["total_samples", len(ron_samples) + len(nie_samples) + len(hmp_samples)])
    writer.writerow(["Ronneby_unstratified_pathways", len(ron)])
    writer.writerow(["Nielsen_unstratified_pathways", len(nie)])
    writer.writerow(["HMP_unstratified_pathways", len(hmp)])
    writer.writerow(["union_unstratified_pathways", len(union)])
    writer.writerow(["shared_unstratified_pathways", len(shared)])
    writer.writerow(["union_no_unmapped_pathways", len(union_no_special)])
    writer.writerow(["shared_no_unmapped_pathways", len(shared_no_special)])
    writer.writerow(["Ronneby_input", str(ronneby_file)])
    writer.writerow(["Nielsen_input", str(nielsen_file)])
    writer.writerow(["HMP_input", str(hmp_file)])

overlap = qcdir / "PHASE5_HUMAnN_pathway_feature_overlap.tsv"
with open(overlap, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["comparison", "n_features"])
    for a in sets:
        writer.writerow([a, len(sets[a])])
    writer.writerow(["Ronneby_intersect_Nielsen", len(sets["Ronneby"] & sets["Nielsen"])])
    writer.writerow(["Ronneby_intersect_HMP", len(sets["Ronneby"] & sets["HMP"])])
    writer.writerow(["Nielsen_intersect_HMP", len(sets["Nielsen"] & sets["HMP"])])
    writer.writerow(["Ronneby_intersect_Nielsen_intersect_HMP", len(shared)])
    writer.writerow(["Ronneby_union_Nielsen_union_HMP", len(union)])

print("DONE")
print(summary)
print(overlap)
