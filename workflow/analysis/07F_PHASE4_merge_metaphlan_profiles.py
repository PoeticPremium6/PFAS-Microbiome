from pathlib import Path
import csv
from collections import defaultdict

REV = Path("/path/to/PFAS_mSystems_revision")
PH4 = REV / "03_taxonomy_metaphlan"
PROFILE_DIR = PH4 / "01_profiles"
OUT_DIR = PH4 / "03_merged"
QC_DIR = PH4 / "04_qc"
TABLE_DIR = REV / "12_tables"

OUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

profile_files = sorted(PROFILE_DIR.glob("*_metaphlan_profile.tsv"))

def sample_id_from_profile(path):
    return path.name.replace("_metaphlan_profile.tsv", "")

def terminal_rank(clade):
    last = clade.split("|")[-1]
    if last.startswith("f__"):
        return "family", last
    if last.startswith("g__"):
        return "genus", last
    if last.startswith("s__"):
        return "species", last
    return None, None

rank_to_data = {
    "family": defaultdict(dict),
    "genus": defaultdict(dict),
    "species": defaultdict(dict),
}

long_rows = []
sample_status = []

for path in profile_files:
    sid = sample_id_from_profile(path)
    parsed_rows = 0
    kept_rows = 0

    if path.stat().st_size == 0:
        sample_status.append([sid, str(path), "EMPTY_PROFILE", 0, 0])
        continue

    with path.open(errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            clade = parts[0]
            try:
                relab = float(parts[2])
            except Exception:
                continue

            parsed_rows += 1
            rank, feature = terminal_rank(clade)
            if rank is None:
                continue

            rank_to_data[rank][feature][sid] = relab
            long_rows.append([sid, rank, feature, clade, relab])
            kept_rows += 1

    status = "COMPLETE" if parsed_rows > 0 and kept_rows > 0 else "NO_KEPT_RANK_ROWS"
    sample_status.append([sid, str(path), status, parsed_rows, kept_rows])

samples = sorted([sample_id_from_profile(p) for p in profile_files])

def write_matrix(rank, path):
    data = rank_to_data[rank]
    features = sorted(data.keys())

    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature_id"] + samples)
        for feat in features:
            w.writerow([feat] + [data[feat].get(sid, 0.0) for sid in samples])

def write_rank_sums(path):
    sums = defaultdict(float)
    for sid, rank, feature, clade, relab in long_rows:
        sums[(sid, rank)] += relab

    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["SampleID", "Rank", "RelativeAbundanceSum"])
        for sid in samples:
            for rank in ["family", "genus", "species"]:
                w.writerow([sid, rank, sums.get((sid, rank), 0.0)])

write_matrix("family", OUT_DIR / "Ronneby_MetaPhlAn_family.tsv")
write_matrix("genus", OUT_DIR / "Ronneby_MetaPhlAn_genus.tsv")
write_matrix("species", OUT_DIR / "Ronneby_MetaPhlAn_species.tsv")

write_matrix("family", TABLE_DIR / "Table_S6_MetaPhlAn_family.tsv")
write_matrix("genus", TABLE_DIR / "Table_S6_MetaPhlAn_genus.tsv")
write_matrix("species", TABLE_DIR / "Table_S6_MetaPhlAn_species.tsv")

with (OUT_DIR / "Ronneby_MetaPhlAn_profiles_long.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["SampleID", "Rank", "FeatureID", "FullClade", "RelativeAbundance"])
    w.writerows(long_rows)

with (QC_DIR / "PHASE4_metaphlan_merge_sample_status.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["SampleID", "Profile", "Status", "ParsedRows", "RankRowsKept"])
    w.writerows(sample_status)

write_rank_sums(QC_DIR / "PHASE4_metaphlan_rank_relative_abundance_sums.tsv")

summary_rows = [
    ["profile_files", len(profile_files)],
    ["samples", len(samples)],
    ["complete_profiles", sum(1 for r in sample_status if r[2] == "COMPLETE")],
    ["noncomplete_profiles", sum(1 for r in sample_status if r[2] != "COMPLETE")],
    ["family_features", len(rank_to_data["family"])],
    ["genus_features", len(rank_to_data["genus"])],
    ["species_features", len(rank_to_data["species"])],
    ["long_rows", len(long_rows)],
    ["family_table", str(OUT_DIR / "Ronneby_MetaPhlAn_family.tsv")],
    ["genus_table", str(OUT_DIR / "Ronneby_MetaPhlAn_genus.tsv")],
    ["species_table", str(OUT_DIR / "Ronneby_MetaPhlAn_species.tsv")],
    ["Table_S6_family", str(TABLE_DIR / "Table_S6_MetaPhlAn_family.tsv")],
    ["Table_S6_genus", str(TABLE_DIR / "Table_S6_MetaPhlAn_genus.tsv")],
    ["Table_S6_species", str(TABLE_DIR / "Table_S6_MetaPhlAn_species.tsv")],
]

with (QC_DIR / "PHASE4_metaphlan_merge_summary.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Field", "Value"])
    w.writerows(summary_rows)

if len(profile_files) != 65 or sum(1 for r in sample_status if r[2] == "COMPLETE") != 65:
    print("WARNING: expected 65 complete profiles.")
else:
    print("OK: 65 complete MetaPhlAn profiles merged.")

print(QC_DIR / "PHASE4_metaphlan_merge_summary.tsv")
