#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("Usage: 96A_build_public_MAG_catalogue.py ROOT")

ROOT = Path(sys.argv[1]).resolve()
SUB = ROOT / "submission"
OUT = SUB / "public_MAGs"
MAN = OUT / "manifests"
META = OUT / "metadata"
CAT = OUT / "catalogue/all_MAGs"
NCBI = OUT / "ncbi"
BATCH = NCBI / "batch_plans"

EXPECTED_MAG_ROWS = 1973

for p in (MAN, META, CAT, NCBI, BATCH):
    p.mkdir(parents=True, exist_ok=True)

def sha256(path):
    h = hashlib.sha256()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def file_sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def norm_col(name):
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")

def pick_col(headers, exact=(), contains=()):
    nmap = {norm_col(h): h for h in headers}
    for x in exact:
        if norm_col(x) in nmap:
            return nmap[norm_col(x)]
    for h in headers:
        nh = norm_col(h)
        if all(token in nh for token in contains):
            return h
    return None

def parse_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in {"na","nan","none","null","missing"}:
        return None
    s = s.replace("%","").replace(",","")
    try:
        return float(s)
    except ValueError:
        return None

def strip_fasta_suffix(name):
    low = name.lower()
    if low.endswith(".gz"):
        name = name[:-3]
        low = low[:-3]
    for ext in (".fasta", ".fna", ".fa", ".fas"):
        if low.endswith(ext):
            return name[:-len(ext)]
    return Path(name).stem

def safe_id(s):
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())
    x = re.sub(r"_+", "_", x).strip("_")
    return x

def fasta_open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8", errors="replace")

def fasta_stats(path):
    seq_count = 0
    total_bp = 0
    min_bp = None
    max_bp = 0
    current_len = 0
    headers = set()
    duplicate_headers = 0
    invalid_sequence_chars = 0
    current_header = None

    def finish_seq():
        nonlocal seq_count, total_bp, min_bp, max_bp, current_len
        if current_header is not None:
            seq_count += 1
            total_bp += current_len
            min_bp = current_len if min_bp is None else min(min_bp, current_len)
            max_bp = max(max_bp, current_len)

    with fasta_open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                finish_seq()
                current_len = 0
                current_header = line[1:].strip().split()[0] if line[1:].strip() else ""
                if current_header in headers:
                    duplicate_headers += 1
                headers.add(current_header)
            else:
                seq = re.sub(r"\s+", "", line).upper()
                if seq and re.search(r"[^ACGTRYSWKMBDHVN.-]", seq):
                    invalid_sequence_chars += 1
                current_len += len(seq.replace("-", "").replace(".", ""))
        finish_seq()

    return {
        "contigs": seq_count,
        "total_bp": total_bp,
        "min_contig_bp": 0 if min_bp is None else min_bp,
        "max_contig_bp": max_bp,
        "duplicate_contig_ids": duplicate_headers,
        "invalid_sequence_lines": invalid_sequence_chars,
    }

# ---------------------------------------------------------------------
# Authoritative MAG-quality table from T1E TSV mirrors.
# ---------------------------------------------------------------------
tsv_dir = SUB / "tables/supplementary/section_tsv_T1E"
quality_candidates = sorted(tsv_dir.glob("Table_S07__*MAG*quality*.tsv"))
if not quality_candidates:
    quality_candidates = sorted(tsv_dir.glob("Table_S07__*.tsv"))

quality_file = None
for p in quality_candidates:
    if "quality" in p.name.lower():
        quality_file = p
        break
if quality_file is None:
    raise RuntimeError(f"Could not locate Table S07 MAG-quality TSV in {tsv_dir}")

with quality_file.open("r", encoding="utf-8-sig", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    headers = reader.fieldnames or []

if len(rows) != EXPECTED_MAG_ROWS:
    raise RuntimeError(
        f"Expected {EXPECTED_MAG_ROWS} authoritative MAG rows; found {len(rows)} in {quality_file}"
    )

id_col = pick_col(
    headers,
    exact=("mag_id","bin_id","genome_id","genome","bin","name"),
)
if id_col is None:
    # Use the first column only if it is unique and nonempty.
    first = headers[0] if headers else None
    if first and len({r.get(first,"") for r in rows}) == len(rows):
        id_col = first
    else:
        raise RuntimeError(f"Could not identify MAG ID column. Headers={headers}")

comp_col = pick_col(headers, exact=("completeness","checkm_completeness","checkm2_completeness"), contains=("complet",))
contam_col = pick_col(headers, exact=("contamination","checkm_contamination","checkm2_contamination"), contains=("contam",))
tax_col = pick_col(headers, exact=("gtdb_taxonomy","classification","taxonomy","gtdb_lineage"), contains=("tax",))
size_col = pick_col(headers, exact=("genome_size","genome_size_bp","size_bp","total_length","length_bp"))

mag_rows = {}
for r in rows:
    mid = str(r.get(id_col,"")).strip()
    if not mid:
        raise RuntimeError("Blank MAG ID in authoritative quality table")
    if mid in mag_rows:
        raise RuntimeError(f"Duplicate MAG ID in authoritative quality table: {mid}")
    mag_rows[mid] = r

# ---------------------------------------------------------------------
# Discover candidate MAG FASTAs.
# ---------------------------------------------------------------------
fasta_exts = (".fa",".fna",".fasta",".fas",".fa.gz",".fna.gz",".fasta.gz",".fas.gz")
skip_parts = {
    ".git","__pycache__","submission","patch_backups","recovery",
    "node_modules","software","references","reference","bundles",
}
skip_contains = (".patch_stage", ".debug", "final_evidence_bundle")

by_stem = defaultdict(list)
files_scanned = 0

for dirpath, dirnames, filenames in os.walk(ROOT):
    rel_parts = Path(dirpath).relative_to(ROOT).parts

    kept = []
    for d in dirnames:
        dl = d.lower()
        if dl in skip_parts:
            continue
        if any(tok in dl for tok in skip_contains):
            continue
        kept.append(d)
    dirnames[:] = kept

    for fn in filenames:
        low = fn.lower()
        if not any(low.endswith(ext) for ext in fasta_exts):
            continue
        p = Path(dirpath) / fn
        files_scanned += 1
        stem = strip_fasta_suffix(fn)
        by_stem[stem].append(p)
        by_stem[stem.lower()].append(p)

inventory = []
selected = {}
missing = []
ambiguous = []

for mid in mag_rows:
    candidates = []
    for key in (mid, mid.lower(), safe_id(mid), safe_id(mid).lower()):
        candidates.extend(by_stem.get(key, []))
    # de-duplicate paths
    candidates = sorted(set(candidates))

    if not candidates:
        missing.append(mid)
        inventory.append([mid, "MISSING", 0, "", "", ""])
        continue

    if len(candidates) == 1:
        selected[mid] = candidates[0]
        inventory.append([mid, "UNIQUE", 1, str(candidates[0].relative_to(ROOT)), file_sha256(candidates[0]), ""])
        continue

    # Multiple copies are acceptable only when byte-identical.
    hash_to_paths = defaultdict(list)
    for p in candidates:
        hash_to_paths[file_sha256(p)].append(p)

    if len(hash_to_paths) == 1:
        # Prefer path containing MAG/bin terminology and shortest relative path.
        chosen = sorted(
            candidates,
            key=lambda p: (
                0 if any(tok in str(p).lower() for tok in ("mag","bin","metabat","drep")) else 1,
                len(str(p.relative_to(ROOT))),
            )
        )[0]
        selected[mid] = chosen
        inventory.append([
            mid, "MULTIPLE_IDENTICAL", len(candidates),
            str(chosen.relative_to(ROOT)), file_sha256(chosen),
            " | ".join(str(p.relative_to(ROOT)) for p in candidates),
        ])
    else:
        ambiguous.append(mid)
        inventory.append([
            mid, "AMBIGUOUS_DISTINCT", len(candidates), "", "",
            " | ".join(str(p.relative_to(ROOT)) for p in candidates),
        ])

with (MAN / "MAG_source_inventory.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["mag_id","status","candidate_count","selected_source","selected_sha256","all_candidates"])
    w.writerows(inventory)

# Hard stop before staging if provenance is not exact.
if missing or ambiguous or len(selected) != EXPECTED_MAG_ROWS:
    (MAN / "phaseM1.status").write_text(
        "\n".join([
            "EXIT_CODE=1",
            "STATUS=DISCOVERY_INCOMPLETE",
            "PHASE=M1_PUBLIC_MAG_CATALOGUE_AND_NCBI_READINESS",
            f"AUTHORITATIVE_MAGS={EXPECTED_MAG_ROWS}",
            f"FASTA_FILES_SCANNED={files_scanned}",
            f"UNIQUE_OR_IDENTICAL_MATCHES={len(selected)}",
            f"MISSING_MAGS={len(missing)}",
            f"AMBIGUOUS_MAGS={len(ambiguous)}",
            "STAGING_RUN=NO",
            "NEXT=RESOLVE_MAG_FASTA_SOURCE_MAPPING",
            "",
        ]),
        encoding="utf-8",
    )
    print("PHASE_M1_DISCOVERY=INCOMPLETE")
    print(f"AUTHORITATIVE_MAGS={EXPECTED_MAG_ROWS}")
    print(f"FASTA_FILES_SCANNED={files_scanned}")
    print(f"UNIQUE_OR_IDENTICAL_MATCHES={len(selected)}")
    print(f"MISSING_MAGS={len(missing)}")
    print(f"AMBIGUOUS_MAGS={len(ambiguous)}")
    raise SystemExit(2)

# ---------------------------------------------------------------------
# Stage all authoritative MAGs without altering sequences.
# ---------------------------------------------------------------------
if CAT.exists():
    shutil.rmtree(CAT)
CAT.mkdir(parents=True, exist_ok=True)

catalogue_rows = []
eligible = []
ineligible = []
public_ids = set()

for i, mid in enumerate(sorted(selected), start=1):
    source = selected[mid]
    out_id = safe_id(mid)
    if not out_id:
        out_id = f"RonnebyMAG{i:05d}"
    if out_id in public_ids:
        raise RuntimeError(f"Sanitized MAG filename collision: {mid} -> {out_id}")
    public_ids.add(out_id)

    dst = CAT / f"{out_id}.fna"

    # Re-materialize compressed inputs; otherwise hard-link/copy unchanged text.
    if source.suffix.lower() == ".gz":
        with fasta_open(source) as src, dst.open("w", encoding="utf-8", newline="\n") as out:
            shutil.copyfileobj(src, out)
        copy_mode = "DECOMPRESSED_COPY"
    else:
        try:
            os.link(source, dst)
            copy_mode = "HARDLINK"
        except OSError:
            shutil.copy2(source, dst)
            copy_mode = "COPY"

    stats = fasta_stats(dst)
    row = mag_rows[mid]

    completeness = parse_float(row.get(comp_col)) if comp_col else None
    contamination = parse_float(row.get(contam_col)) if contam_col else None
    taxonomy = str(row.get(tax_col,"")).strip() if tax_col else ""

    # Prefer computed size from sequence itself.
    total_bp = stats["total_bp"]

    ncbi_ok = (
        completeness is not None
        and completeness >= 90.0
        and total_bp >= 100000
        and stats["min_contig_bp"] >= 200
        and stats["duplicate_contig_ids"] == 0
        and stats["invalid_sequence_lines"] == 0
    )

    reason = []
    if completeness is None:
        reason.append("missing_checkm_completeness")
    elif completeness < 90:
        reason.append("checkm_lt_90")
    if total_bp < 100000:
        reason.append("assembly_lt_100kb")
    if stats["min_contig_bp"] < 200:
        reason.append("contig_lt_200bp")
    if stats["duplicate_contig_ids"] > 0:
        reason.append("duplicate_contig_ids")
    if stats["invalid_sequence_lines"] > 0:
        reason.append("invalid_sequence_characters")

    proposed_isolate = f"RonnebyMAG{i:05d}"

    out = {
        "mag_id": mid,
        "public_filename": dst.name,
        "proposed_isolate": proposed_isolate,
        "source_path": str(source.relative_to(ROOT)),
        "copy_mode": copy_mode,
        "sha256": file_sha256(dst),
        "contigs": stats["contigs"],
        "total_bp": total_bp,
        "min_contig_bp": stats["min_contig_bp"],
        "max_contig_bp": stats["max_contig_bp"],
        "duplicate_contig_ids": stats["duplicate_contig_ids"],
        "invalid_sequence_lines": stats["invalid_sequence_lines"],
        "checkm_or_checkm2_completeness": "" if completeness is None else completeness,
        "checkm_or_checkm2_contamination": "" if contamination is None else contamination,
        "gtdb_lineage": taxonomy,
        "ncbi_sequence_eligibility": "ELIGIBLE" if ncbi_ok else "NOT_ELIGIBLE",
        "ncbi_sequence_eligibility_reason": "PASS" if ncbi_ok else ";".join(reason),
    }
    catalogue_rows.append(out)
    (eligible if ncbi_ok else ineligible).append(out)

fields = list(catalogue_rows[0].keys())
with (MAN / "MAG_catalogue_manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(catalogue_rows)

with (NCBI / "NCBI_MAG_sequence_eligibility.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(catalogue_rows)

# NCBI taxonomy-review working table.
with (NCBI / "NCBI_taxonomy_review_request.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["proposed_isolate","original_mag_id","GTDB_lineage_unmodified","NCBI_organism_name_to_be_confirmed"])
    for r in eligible:
        w.writerow([r["proposed_isolate"], r["mag_id"], r["gtdb_lineage"], ""])

# NCBI BioSample working metadata table. This is intentionally not claimed to
# be the portal template until current official MIMAG template is downloaded.
with (NCBI / "NCBI_MIMAG_BioSample_working_metadata.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([
        "proposed_isolate","original_mag_id","sample_name","organism",
        "sample_type","isolation_source","geo_loc_name","collection_date",
        "derived_from_SRA_run_or_physical_BioSample","BioProject","BioSample",
        "status",
    ])
    for r in eligible:
        w.writerow([
            r["proposed_isolate"], r["mag_id"], r["proposed_isolate"], "",
            "metagenomic assembly", "human gut metagenome", "", "",
            "", "", "", "NEEDS_ACCESSIONS_AND_NCBI_TAXON_NAME",
        ])

# Genome batch plan: <=400 assemblies per batch.
with (NCBI / "NCBI_genome_batch_plan.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["batch_id","batch_row","proposed_isolate","original_mag_id","filename","BioProject","BioSample","status"])
    for idx, r in enumerate(eligible, start=1):
        batch_no = (idx - 1) // 400 + 1
        row_no = (idx - 1) % 400 + 1
        w.writerow([
            f"batch_{batch_no:03d}", row_no, r["proposed_isolate"],
            r["mag_id"], r["public_filename"], "", "",
            "NEEDS_BIOPROJECT_BIOSAMPLE_AND_HEADER_ACCESSIONS",
        ])

# Split explicit batch manifests.
batch_groups = defaultdict(list)
for idx, r in enumerate(eligible, start=1):
    batch_no = (idx - 1) // 400 + 1
    batch_groups[batch_no].append(r)

for batch_no, batch_rows in batch_groups.items():
    p = BATCH / f"batch_{batch_no:03d}.tsv"
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["row","proposed_isolate","original_mag_id","filename","checkm_completeness","total_bp","gtdb_lineage"])
        for j, r in enumerate(batch_rows, start=1):
            w.writerow([
                j, r["proposed_isolate"], r["mag_id"], r["public_filename"],
                r["checkm_or_checkm2_completeness"], r["total_bp"], r["gtdb_lineage"],
            ])

# Checksums for the complete public catalogue.
with (MAN / "SHA256SUMS_MAG_FASTA").open("w", encoding="utf-8") as fh:
    for r in catalogue_rows:
        fh.write(f'{r["sha256"]}  catalogue/all_MAGs/{r["public_filename"]}\n')

# Human-readable README.
(OUT / "README_MAG_PUBLIC_RELEASE.md").write_text(
    f"""# PFAS mSystems public MAG release package

This directory contains the authoritative genome-resolved catalogue used for
the manuscript and is separate from the frozen manuscript/figure submission.

## Catalogue

- Authoritative MAGs: {len(catalogue_rows)}
- FASTA files staged: {len(catalogue_rows)}
- Current NCBI sequence-eligible MAGs: {len(eligible)}
- Current NCBI sequence-ineligible MAGs: {len(ineligible)}

`catalogue/all_MAGs/` preserves every authoritative MAG regardless of current
NCBI MAG acceptance thresholds.

## NCBI readiness

The NCBI eligibility table applies sequence-level gates only:
- CheckM or CheckM2 completeness >= 90%
- total assembly size >= 100,000 nt
- every contig >= 200 nt
- no duplicate contig identifiers
- no invalid sequence characters detected

Passing these checks does NOT mean a MAG can yet be uploaded. NCBI additionally
requires a BioProject, one organism-specific MIMAG BioSample per MAG, an NCBI
organism name confirmed from the unmodified GTDB lineage, and SRA run
accessions (or an appropriate physical mixed-sample BioSample) linking the MAG
to its source data.

The working BioSample and taxonomy-review tables in `ncbi/` intentionally leave
those accession fields blank until authoritative public accessions are known.

Genome submission batches are planned at no more than 400 MAGs each.

## Provenance

`manifests/MAG_source_inventory.tsv` records the repository source FASTA chosen
for every MAG. Multiple source copies were accepted only when byte-identical.
No source MAG file was modified.

`manifests/MAG_catalogue_manifest.tsv` records sequence statistics, quality,
taxonomy context, final public filename and SHA256 for every MAG.
""",
    encoding="utf-8",
)

status = (
    "READY_FOR_PUBLIC_CATALOGUE_QC"
    if len(catalogue_rows) == EXPECTED_MAG_ROWS
    else "FAIL"
)

(MAN / "phaseM1.status").write_text(
    "\n".join([
        "EXIT_CODE=0",
        f"STATUS={status}",
        "PHASE=M1_PUBLIC_MAG_CATALOGUE_AND_NCBI_READINESS",
        f"AUTHORITATIVE_MAGS={EXPECTED_MAG_ROWS}",
        f"FASTAS_STAGED={len(catalogue_rows)}",
        f"NCBI_SEQUENCE_ELIGIBLE={len(eligible)}",
        f"NCBI_SEQUENCE_INELIGIBLE={len(ineligible)}",
        f"NCBI_BATCHES_PLANNED={len(batch_groups)}",
        "SOURCE_MAG_FILES_MODIFIED=NO",
        "T3_T3R3_FREEZES_MODIFIED=NO",
        "NCBI_ACCESSION_METADATA_COMPLETE=NO",
        "NEXT=REVIEW_MAG_CATALOGUE_AND_RESOLVE_SOURCE_ACCESSIONS",
        "",
    ]),
    encoding="utf-8",
)

print("PHASE_M1=PASS")
print(f"AUTHORITATIVE_MAGS={EXPECTED_MAG_ROWS}")
print(f"FASTAS_STAGED={len(catalogue_rows)}")
print(f"NCBI_SEQUENCE_ELIGIBLE={len(eligible)}")
print(f"NCBI_SEQUENCE_INELIGIBLE={len(ineligible)}")
print(f"NCBI_BATCHES_PLANNED={len(batch_groups)}")
print("STATUS=READY_FOR_PUBLIC_CATALOGUE_QC")
