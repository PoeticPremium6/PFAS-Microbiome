# PFAS–Microbiome: public release gate (v1.2.0)

**Status: pre-release, not approved for archival until every relevant gate is signed off.**

## 1. Integrity and traceability

- [ ] `git status` is clean; changes reviewed by a contributor with authority to publish.
- [ ] `VERSION`, README and `CITATION.cff` agree; no fabricated DOI or release date.
- [ ] The public results and figure source tables correspond to the final manuscript figures.
- [ ] Each public analysis script declares required inputs and outputs; published software dependencies and external databases are documented.
- [ ] A second environment or independent reviewer reproduces the publication-facing figures from permitted source tables; record versions, expected outputs and tolerances.
- [ ] The repository clearly distinguishes reproducibility from restricted-source reanalysis (raw cohort data not public).
- [ ] Generated figures, figure identifiers and caption claims match the submitted manuscript.

## 2. Human-participant data governance — **mandatory manual review**

- [ ] Permission from relevant cohort/institutional data controllers for this specific public deposit is documented.
- [ ] No raw or host-depleted sequencing reads, read-level alignments, human variants, participant identifiers, restricted metadata or individual PFAS/clinical measurements have been committed.
- [ ] Every *sample-level* microbial abundance/functional matrix is separately assessed for membership/linkage/re-identification risk (absence of names is insufficient).
- [ ] Any analysis tables are reviewed for quasi-identifiers, outliers and small cells and aggregation is considered when appropriate.
- [ ] Human-decontamination of the 321 MAG assemblies is independently confirmed before NCBI public release; this is separate from the repository scan.
- [ ] `git log --all --name-only` and git history are reviewed; deleting a file in a later commit does not delete it from previous public history.
- [ ] No credentials, tokens or private filesystem paths; enable GitHub security features when available.
- [ ] `python3 scripts/release_qa.py` exits successfully; its warnings have been reviewed, not merely ignored.

## 3. Authorship, citations and licenses

- [ ] The manuscript/software author list, contributor approval and ownership have been confirmed.
- [ ] Choose a software license only with all appropriate institutional/coauthor rights confirmed; do not add a license automatically.
- [ ] Data re-use permissions/license are decided separately from the software license.
- [ ] Third-party dataset, code, model and database licensing/attribution checked.
- [ ] DOI and accession links are added **only** when real identifiers are assigned.

## 4. Release, Zenodo and manuscript

- [ ] Zenodo GitHub integration enabled **before** creating the GitHub release.
- [ ] Create GitHub release `v1.2.0` only when checks above complete.
- [ ] Confirm software archive was actually ingested by Zenodo; capture version DOI and concept DOI.
- [ ] Zenodo *Dataset* record for explicitly approved processed tables created separately, with README, data dictionary, checksums and disclosure review.
- [ ] NCBI BioProject PRJNA1522707, BioSample submission SUB16581910 and approved MAG assembly accessions reconciled; do not invent accessions.
- [ ] Manuscript Data Availability and code links match the final published artifacts and any privacy limitations.

*Automated checks are preliminary safeguards, not permission to disclose or proof that results are reproducible.*
