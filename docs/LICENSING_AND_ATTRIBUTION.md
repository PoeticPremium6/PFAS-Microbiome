# Licensing, provenance and reusable outputs

This project has **separate software and research-data licences**:

| Output | Proposed public licence | Location |
|---|---|---|
| Original analysis and figure-generating software | GPL-3.0-only | `LICENSE` |
| Original processed summary tables approved for release | CC BY 4.0 | `DATA_LICENSE.md` |
| Original documentation and original display material where rights permit | CC BY 4.0 | `DATA_LICENSE.md` |
| Third-party tools, reference data, models and annotations | Their original licences | Upstream sources |
| Restricted original sequencing and participant-level material | Not made public | Not included |

## Required user attribution

Reuse of original processed tables under CC BY 4.0 requires proper credit, reference to the licence and an indication of changes. Please cite the paper and the specific Zenodo dataset version DOI when available. Code distributed under GPL-3.0-only requires compliance with GNU GPL terms, including relevant source-distribution obligations. Citation requests are distinct from licence obligations.

**No separate permission is required for commercial use permitted by those licenses.** These terms do not grant the right to misuse confidential source data or redistribute third-party references.

## Boundaries and provenance

1. `results/machine_readable/` primarily contains analysis summaries and derived results. Neither its directory name nor a successful automated QA run constitutes ethical clearance.
2. For `Table_S07__MAG_abundance_context.tsv`, `Table_S05__Beta_diversity.tsv` and public-context tables, retain the file-specific privacy/redistribution decision in the deposit records.
3. `results/machine_readable/Table_S04__HUMAnN_associations.tsv` contains aggregate-style association columns; see `scripts/check_humann_summary.py` to check tabular integrity without copying records into logs.
4. The 321 MAG assemblies are intended for NCBI, under BioProject PRJNA1522707, once that independent submission process is complete. Do not put unapproved genome assemblies in the GitHub repository or Zenodo dataset.
5. The sequence reads and individual-level PFAS and clinical metadata are withheld from public release as described in the manuscript. Repository history needs a separate review.
6. The public-data deposit must include only files individually cleared for public redistribution. `scripts/zenodo_inventory.py` is an inventory, not an approval mechanism.

## Release citation work

- Ensure `CITATION.cff` is accurate and lists the actual software authors.
- Once issued, record a **software version DOI** and a separate **processed-data dataset DOI**. Do not invent either DOI.
- Add the MAG `SAMN`/assembly accessions when NCBI review and assembly submission finish.
- Record manuscript methods, software/database versions, checksums and the exact source commit.

Copyright remains with the relevant owners. Confirm that any relevant University of Glasgow and collaborating institution arrangements, author agreements and original-source terms permit the offered licences **before publishing this patch**.
