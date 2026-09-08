from pathlib import Path
import csv
import re

REV = Path("/path/to/PFAS_mSystems_revision")
PH4 = REV / "03_taxonomy_metaphlan"
PH6 = REV / "05_humann_functional_profiles"
PH6_MANIFEST = PH6 / "00_manifest" / "Ronneby_HUMAnN_manifest.tsv"
PH4_MANIFEST = PH4 / "00_manifest" / "Ronneby_MetaPhlAn_manifest.tsv"
PROFILE_DIR = PH4 / "01_profiles"

PH6_MANIFEST.parent.mkdir(parents=True, exist_ok=True)

def norm(x):
    return re.sub(r"[^a-z0-9]", "", x.lower())

def pick_col(header, patterns):
    hnorm = [norm(h) for h in header]
    for p in patterns:
        pn = norm(p)
        for i, h in enumerate(hnorm):
            if pn == h or pn in h:
                return i
    return None

rows = []
if PH4_MANIFEST.exists():
    with PH4_MANIFEST.open(errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sid_i = pick_col(header, ["SampleID", "sample", "sample_id"])
        r1_i = pick_col(header, ["R1", "Read1", "fastq1", "fq1"])
        r2_i = pick_col(header, ["R2", "Read2", "fastq2", "fq2"])

        if sid_i is None:
            sid_i = 0
        if r1_i is None or r2_i is None:
            raise SystemExit(f"Could not detect R1/R2 columns in {PH4_MANIFEST}: {header}")

        for r in reader:
            if not r or len(r) <= max(sid_i, r1_i, r2_i):
                continue
            sid = r[sid_i].strip()
            r1 = r[r1_i].strip()
            r2 = r[r2_i].strip()
            if sid:
                profile = PROFILE_DIR / f"{sid}_metaphlan_profile.tsv"
                status = "READY" if Path(r1).exists() and Path(r2).exists() and profile.exists() and profile.stat().st_size > 0 else "MISSING_INPUT"
                rows.append([sid, r1, r2, str(profile), status])
else:
    raise SystemExit(f"Missing required Phase 4 manifest: {PH4_MANIFEST}")

rows = sorted(rows, key=lambda x: x[0])

with PH6_MANIFEST.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["SampleID", "R1", "R2", "MetaPhlAnProfile", "Status"])
    w.writerows(rows)

ready = sum(1 for r in rows if r[4] == "READY")
print(f"manifest={PH6_MANIFEST}")
print(f"samples={len(rows)}")
print(f"ready={ready}")
if ready != 65:
    print("WARNING: expected 65 READY samples.")
