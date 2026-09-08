#!/usr/bin/env python3
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("Usage: 98A_build_MAG_sample_read_provenance.py ROOT")

ROOT = Path(sys.argv[1]).resolve()
PUB = ROOT / "submission/public_MAGs"
M1 = PUB / "manifests/MAG_catalogue_manifest.tsv"
OUT = PUB / "ncbi/sample_read_provenance"
OUT.mkdir(parents=True, exist_ok=True)

if not M1.is_file():
    raise RuntimeError(f"Missing M1 manifest: {M1}")

with M1.open("r", encoding="utf-8-sig", newline="") as fh:
    mags = list(csv.DictReader(fh, delimiter="\t"))

if len(mags) != 1973:
    raise RuntimeError(f"Expected 1973 MAGs; found {len(mags)}")

def source_id(mid):
    m = re.match(r"^([A-H]\d{2})_bin(?:\.|_|$)", mid)
    if not m:
        raise RuntimeError(f"Cannot parse source sample from MAG ID: {mid}")
    return m.group(1)

source_ids = sorted({source_id(r["mag_id"]) for r in mags})
if len(source_ids) != 65:
    raise RuntimeError(f"Expected 65 MAG source samples; found {len(source_ids)}: {source_ids}")

mag_counts = defaultdict(int)
eligible_counts = defaultdict(int)
for r in mags:
    sid = source_id(r["mag_id"])
    mag_counts[sid] += 1
    if r.get("ncbi_sequence_eligibility") == "ELIGIBLE":
        eligible_counts[sid] += 1

# Exact MAG -> source mapping.
with (OUT / "MAG_to_source_sample.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["mag_id","source_sample","ncbi_sequence_eligibility","public_filename"])
    for r in mags:
        w.writerow([
            r["mag_id"], source_id(r["mag_id"]),
            r.get("ncbi_sequence_eligibility",""),
            r.get("public_filename","")
        ])

# Target only the Ronneby provenance/read-processing parts of the repo.
targets = [
    ROOT / "01_metadata_freeze",
    ROOT / "02_read_qc_and_decontamination_stats",
    ROOT / "03_taxonomy_metaphlan/00_manifest",
    ROOT / "05_humann_functional_profiles/00_manifest",
    ROOT / "06_MAG_resource_integration",
]
top_files = [
    ROOT / "PFAS_Metadata.csv",
]

TEXT_EXTS = {".tsv",".csv",".txt",".md",".json",".yaml",".yml",".sh",".py"}
MAX_TEXT_BYTES = 50 * 1024 * 1024

RUN_RE = re.compile(r"\b(?:SRR|ERR|DRR)\d+\b", re.I)
SAMPLE_RE = re.compile(r"\b(?:SAMN|SAMEA|SAMD)\d+\b", re.I)
PROJECT_RE = re.compile(r"\b(?:PRJNA|PRJEB|PRJDB)\d+\b", re.I)
FASTQ_RE = re.compile(
    r"(?P<path>(?:/|\.{0,2}/)?[A-Za-z0-9_./+\-]+?\.(?:fastq|fq)(?:\.gz)?)",
    re.I
)

def exact_samples(text):
    return sorted(set(re.findall(r"(?<![A-Za-z0-9])([A-H]\d{2})(?![A-Za-z0-9])", text)))

def read_class(text):
    low = text.lower()
    if any(x in low for x in ("dehost", "host_removed", "host-removed", "decontam", "cleaned")):
        return "DEHOSTED"
    if "raw" in low or "login_" in low:
        return "RAW"
    return "UNCLASSIFIED"

def mate_class(text):
    low = text.lower()
    if re.search(r"(?:^|[_./-])r?1(?:[_./-]|$)", low) or "_1.fastq" in low or "_1.fq" in low:
        return "R1"
    if re.search(r"(?:^|[_./-])r?2(?:[_./-]|$)", low) or "_2.fastq" in low or "_2.fq" in low:
        return "R2"
    return "UNKNOWN"

def resolve_path(raw):
    p = Path(raw)
    candidates = []
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(ROOT / p)
        candidates.append(ROOT / "02_read_qc_and_decontamination_stats" / p)
    for c in candidates:
        if c.exists():
            return str(c), "YES", c.stat().st_size
    return raw, "NO", ""

evidence = []
files_scanned = 0

def scan_file(p):
    global files_scanned
    try:
        if not p.is_file() or p.stat().st_size <= 0 or p.stat().st_size > MAX_TEXT_BYTES:
            return
    except OSError:
        return
    if p.suffix.lower() not in TEXT_EXTS:
        return

    files_scanned += 1
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            for line_no, line in enumerate(fh, start=1):
                sids = exact_samples(line)
                if not sids:
                    continue
                fastqs = [m.group("path") for m in FASTQ_RE.finditer(line)]
                runs = sorted(set(x.upper() for x in RUN_RE.findall(line)))
                sams = sorted(set(x.upper() for x in SAMPLE_RE.findall(line)))
                prjs = sorted(set(x.upper() for x in PROJECT_RE.findall(line)))
                if not (fastqs or runs or sams or prjs or "fastq" in line.lower() or "read" in line.lower()):
                    continue

                snippet = re.sub(r"\s+", " ", line.strip())
                if len(snippet) > 700:
                    snippet = snippet[:697] + "..."

                for sid in sids:
                    if sid not in source_ids:
                        continue
                    if fastqs:
                        for fq in fastqs:
                            resolved, exists, size = resolve_path(fq)
                            evidence.append({
                                "source_sample": sid,
                                "evidence_file": str(p.relative_to(ROOT)),
                                "line": line_no,
                                "record_type": "FASTQ",
                                "read_class": read_class(line + " " + str(p)),
                                "mate": mate_class(fq),
                                "value": fq,
                                "resolved_path": resolved,
                                "exists": exists,
                                "bytes": size,
                                "sra_runs": ",".join(runs),
                                "biosamples": ",".join(sams),
                                "bioprojects": ",".join(prjs),
                                "snippet": snippet,
                            })
                    elif runs or sams or prjs:
                        evidence.append({
                            "source_sample": sid,
                            "evidence_file": str(p.relative_to(ROOT)),
                            "line": line_no,
                            "record_type": "ACCESSION",
                            "read_class": read_class(line + " " + str(p)),
                            "mate": "",
                            "value": "",
                            "resolved_path": "",
                            "exists": "",
                            "bytes": "",
                            "sra_runs": ",".join(runs),
                            "biosamples": ",".join(sams),
                            "bioprojects": ",".join(prjs),
                            "snippet": snippet,
                        })
    except Exception:
        return

for d in targets:
    if not d.exists():
        continue
    if d.is_file():
        scan_file(d)
        continue
    for dirpath, dirnames, filenames in os.walk(d):
        dirnames[:] = [
            x for x in dirnames
            if x not in {".git","__pycache__","logs","patch_backups"}
            and not x.startswith(".")
        ]
        for fn in filenames:
            scan_file(Path(dirpath) / fn)

for p in top_files:
    scan_file(p)

efields = [
    "source_sample","evidence_file","line","record_type","read_class","mate",
    "value","resolved_path","exists","bytes","sra_runs","biosamples",
    "bioprojects","snippet"
]
with (OUT / "targeted_read_accession_evidence.tsv").open(
    "w", encoding="utf-8", newline=""
) as fh:
    w = csv.DictWriter(fh, fieldnames=efields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(evidence)

# Aggregate per source sample.
by_sample = defaultdict(list)
for e in evidence:
    by_sample[e["source_sample"]].append(e)

summary_rows = []
for sid in source_ids:
    evs = by_sample.get(sid, [])

    def paths(cls, mate):
        return sorted({
            e["resolved_path"] for e in evs
            if e["record_type"] == "FASTQ"
            and e["read_class"] == cls
            and e["mate"] == mate
            and e["resolved_path"]
        })

    def existing_paths(cls, mate):
        return sorted({
            e["resolved_path"] for e in evs
            if e["record_type"] == "FASTQ"
            and e["read_class"] == cls
            and e["mate"] == mate
            and e["exists"] == "YES"
        })

    raw_r1 = paths("RAW","R1")
    raw_r2 = paths("RAW","R2")
    dh_r1 = paths("DEHOSTED","R1")
    dh_r2 = paths("DEHOSTED","R2")

    eraw_r1 = existing_paths("RAW","R1")
    eraw_r2 = existing_paths("RAW","R2")
    edh_r1 = existing_paths("DEHOSTED","R1")
    edh_r2 = existing_paths("DEHOSTED","R2")

    runs = sorted({
        x for e in evs
        for x in e["sra_runs"].split(",") if x
    })
    sams = sorted({
        x for e in evs
        for x in e["biosamples"].split(",") if x
    })
    prjs = sorted({
        x for e in evs
        for x in e["bioprojects"].split(",") if x
    })

    if edh_r1 and edh_r2:
        read_status = "PAIRED_DEHOSTED_LOCAL"
    elif eraw_r1 and eraw_r2:
        read_status = "PAIRED_RAW_LOCAL_ONLY"
    elif raw_r1 or raw_r2 or dh_r1 or dh_r2:
        read_status = "READ_PATHS_FOUND_NOT_LOCALLY_COMPLETE"
    else:
        read_status = "NO_FASTQ_PATH_FOUND"

    if runs:
        accession_status = "SRA_RUN_PRESENT"
    elif sams:
        accession_status = "BIOSAMPLE_PRESENT_NO_RUN"
    else:
        accession_status = "NO_PUBLIC_SOURCE_ACCESSION"

    summary_rows.append({
        "source_sample": sid,
        "all_MAGs": mag_counts[sid],
        "ncbi_eligible_MAGs": eligible_counts[sid],
        "read_status": read_status,
        "accession_status": accession_status,
        "raw_R1_candidates": " | ".join(raw_r1),
        "raw_R2_candidates": " | ".join(raw_r2),
        "dehosted_R1_candidates": " | ".join(dh_r1),
        "dehosted_R2_candidates": " | ".join(dh_r2),
        "existing_raw_R1": " | ".join(eraw_r1),
        "existing_raw_R2": " | ".join(eraw_r2),
        "existing_dehosted_R1": " | ".join(edh_r1),
        "existing_dehosted_R2": " | ".join(edh_r2),
        "SRA_runs": ",".join(runs),
        "physical_BioSamples": ",".join(sams),
        "BioProjects": ",".join(prjs),
        "evidence_records": len(evs),
        "evidence_files": " | ".join(sorted({e["evidence_file"] for e in evs})),
    })

sfields = list(summary_rows[0].keys())
with (OUT / "sample_read_provenance.tsv").open(
    "w", encoding="utf-8", newline=""
) as fh:
    w = csv.DictWriter(fh, fieldnames=sfields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary_rows)

# Minimal public-deposition working sheet. Do not copy phenotype or personal
# metadata into this layer.
with (OUT / "source_sample_deposition_working.tsv").open(
    "w", encoding="utf-8", newline=""
) as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([
        "source_sample","all_MAGs","ncbi_eligible_MAGs","read_status",
        "preferred_public_read_source_review","SRA_run","physical_BioSample",
        "BioProject","collection_date","geo_loc_name","isolation_source",
        "human_read_screening_confirmed","consent_public_release_confirmed",
        "status"
    ])
    for r in summary_rows:
        preferred = (
            "REVIEW_DEHOSTED_PAIR"
            if r["read_status"] == "PAIRED_DEHOSTED_LOCAL"
            else "REVIEW_RAW_PAIR"
            if r["read_status"] == "PAIRED_RAW_LOCAL_ONLY"
            else "RESOLVE_READ_SOURCE"
        )
        w.writerow([
            r["source_sample"], r["all_MAGs"], r["ncbi_eligible_MAGs"],
            r["read_status"], preferred, r["SRA_runs"],
            r["physical_BioSamples"], r["BioProjects"],
            "", "", "human gut metagenome", "", "",
            "NEEDS_DEPOSITION_METADATA_REVIEW"
        ])

counts = defaultdict(int)
acc_counts = defaultdict(int)
for r in summary_rows:
    counts[r["read_status"]] += 1
    acc_counts[r["accession_status"]] += 1

(OUT / "phaseM3.status").write_text(
    "\n".join([
        "EXIT_CODE=0",
        "STATUS=READY_FOR_TARGETED_READ_PROVENANCE_REVIEW",
        "PHASE=M3_TARGETED_MAG_SAMPLE_READ_PROVENANCE",
        "SOURCE_SAMPLES=65",
        "AUTHORITATIVE_MAGS=1973",
        "NCBI_ELIGIBLE_MAGS=1383",
        f"TARGET_FILES_SCANNED={files_scanned}",
        f"EVIDENCE_RECORDS={len(evidence)}",
        f"PAIRED_DEHOSTED_LOCAL={counts['PAIRED_DEHOSTED_LOCAL']}",
        f"PAIRED_RAW_LOCAL_ONLY={counts['PAIRED_RAW_LOCAL_ONLY']}",
        f"READ_PATHS_FOUND_NOT_LOCALLY_COMPLETE={counts['READ_PATHS_FOUND_NOT_LOCALLY_COMPLETE']}",
        f"NO_FASTQ_PATH_FOUND={counts['NO_FASTQ_PATH_FOUND']}",
        f"SAMPLES_WITH_SRA_RUN={acc_counts['SRA_RUN_PRESENT']}",
        f"SAMPLES_WITH_BIOSAMPLE_NO_RUN={acc_counts['BIOSAMPLE_PRESENT_NO_RUN']}",
        f"SAMPLES_WITH_NO_PUBLIC_SOURCE_ACCESSION={acc_counts['NO_PUBLIC_SOURCE_ACCESSION']}",
        "RAW_READS_COPIED=NO",
        "DEHOSTED_READS_COPIED=NO",
        "FASTA_HEADERS_MODIFIED=NO",
        "T3_T3R3_MODIFIED=NO",
        "NEXT=SELECT_PUBLIC_READ_SOURCE_AND_COMPLETE_SOURCE_DEPOSITION_METADATA",
        ""
    ]),
    encoding="utf-8"
)

print("PHASE_M3=PASS")
print("STATUS=READY_FOR_TARGETED_READ_PROVENANCE_REVIEW")
print("SOURCE_SAMPLES=65")
print("AUTHORITATIVE_MAGS=1973")
print("NCBI_ELIGIBLE_MAGS=1383")
print(f"TARGET_FILES_SCANNED={files_scanned}")
print(f"EVIDENCE_RECORDS={len(evidence)}")
for key in [
    "PAIRED_DEHOSTED_LOCAL",
    "PAIRED_RAW_LOCAL_ONLY",
    "READ_PATHS_FOUND_NOT_LOCALLY_COMPLETE",
    "NO_FASTQ_PATH_FOUND",
]:
    print(f"{key}={counts[key]}")
print(f"SAMPLES_WITH_SRA_RUN={acc_counts['SRA_RUN_PRESENT']}")
print(f"SAMPLES_WITH_BIOSAMPLE_NO_RUN={acc_counts['BIOSAMPLE_PRESENT_NO_RUN']}")
print(f"SAMPLES_WITH_NO_PUBLIC_SOURCE_ACCESSION={acc_counts['NO_PUBLIC_SOURCE_ACCESSION']}")
