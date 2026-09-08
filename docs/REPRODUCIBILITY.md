# Reproducibility

## Scope

This repository is a public reproducibility companion to a human-cohort study.
It is not intended to reconstruct restricted participant-level inputs from
GitHub alone.

## Code

`workflow/analysis/` contains the frozen public analysis/provenance snapshot.
`workflow/final_figures/` contains current journal-facing figure builders that
were available in the project at publication-release time.

Additional post-v1.0 public analyses, when available, are placed under
`workflow/additional_analysis/`.

## Environments

The original analyses used multiple R/Python/HPC environments. Environment
snapshots are retained as provenance, not as a fully portable lockfile. Local
conda-build paths are removed from the public package snapshot.

## Integrity

`manifests/public_release_manifest.tsv` and root `SHA256SUMS` provide release
checksums. `manifests/syntax_audit.tsv` records Python and shell syntax checks.

## Scientific scope

Publication finalization changes packaging and documentation only. It does not
recalculate the frozen statistical results.
