# PFAS–Microbiome processed-data deposit plan

## Purpose and scope

Create a Zenodo **Dataset** record separate from the GitHub code **Software** archive. Eligible candidates may include public/cleared *aggregate or derived* MetaPhlAn, HUMAnN, AGORA2, statistical summary and figure-source tables. Merely being processed, host-depleted, coded, or de-identified does **not** guarantee public-release permission.

**Do not upload:** raw or host-filtered reads, human variants, BAM/CRAM/FASTQ files, participant-level PFAS/clinical metadata, participant linkage keys, MAG assembly FASTAs awaiting NCBI release, or any sample-level file not individually approved for disclosure.

## Inventory — review only, never automatically package

From the repo root:

```bash
python3 scripts/zenodo_inventory.py --out zenodo_candidate_inventory.tsv
```

The TSV contains relative paths, file sizes, SHA-256 hashes, header names (when CSV/TSV), review flags, and a default `PENDING` value for every file. Keep the working inventory private until headers/filenames have been checked and approvals recorded. **The script does not approve, copy, package, upload or publish any data.**

## Dataset record contents (only after sign-off)

1. `README_DATA.md`: study context; analysis units; generation pipeline; limitations; permitted use.
2. `DATA_DICTIONARY.tsv`: each released file's columns, definitions, units, transformation, missingness and interpretation.
3. Approved summary/source tables, with explicit version and corresponding figure/table references.
4. `SHA256SUMS.txt`: cryptographic checksums of released files.
5. Citation and provenance: manuscript status, GitHub commit/release, BioProject PRJNA1522707 and references to external databases and model versions.

## Suggested Zenodo metadata

- **Title:** PFAS–microbiome systems analysis of the Ronneby cohort: approved processed analysis tables
- **Resource type:** Dataset
- **Creators:** Confirm the authorised author/contributor list before publishing; do not assume the software author list is identical.
- **Related identifiers:** GitHub software release DOI (when assigned), BioProject PRJNA1522707 and manuscript DOI (when assigned).
- **Description:** Derived taxonomic, functional and metabolic-model analysis outputs supporting the manuscript; no raw sequencing data or restricted individual-level clinical/PFAS data; state the exact tables made available.
- **License:** Choose only after institutional/coauthor approval for the specific released tables.
- **DOI:** A draft can reserve a DOI before publication, but it is registered only upon publishing.

## Verification

Keep the Zenodo deposit as a **draft** until files, ethics/privacy review, manuscript authorisation, checksums, data dictionary and citations are complete. Inspect the uploaded contents and metadata before publishing. Do not duplicate NCBI genome FASTAs in this dataset deposit.
