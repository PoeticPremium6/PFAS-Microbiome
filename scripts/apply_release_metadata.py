#!/usr/bin/env python3
"""Safely align publication version metadata while retaining existing CFF authors.

VERSION is the local source of truth. The helper adds a clearly marked README
release section when one does not already exist. This script does not publish.
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXPECTED_VERSION = '1.2.0'
MARKER = '<!-- RELEASE_STANDARDS_v1.2 -->'


def main() -> None:
    readme_path = ROOT / 'README.md'
    cff_path = ROOT / 'CITATION.cff'
    version_path = ROOT / 'VERSION'
    for path in (readme_path, cff_path, version_path):
        if not path.is_file():
            raise SystemExit(f'Missing {path}. Apply in the PFAS-Microbiome repository checkout.')

    raw_version = version_path.read_text(encoding='utf-8').strip()
    if not re.fullmatch(r'v?\d+\.\d+\.\d+', raw_version):
        raise SystemExit(f'Unexpected VERSION format: {raw_version!r}')
    version = raw_version.removeprefix('v')
    if version != EXPECTED_VERSION:
        raise SystemExit(f'Expected VERSION={EXPECTED_VERSION}, found {version}; inspect before changing metadata.')

    cff = cff_path.read_text(encoding='utf-8')
    old_versions = re.findall(r'^version:\s*["\']?(\d+\.\d+\.\d+)', cff, flags=re.M)
    if len(old_versions) != 1:
        raise SystemExit('CITATION.cff must contain exactly one top-level version field; manual review required.')
    old_version = old_versions[0]
    if old_version not in ('1.1.0', version):
        raise SystemExit(f'Unexpected CFF version {old_version}; manual review required.')
    updated_cff, count = re.subn(
        r'^version:\s*["\']?\d+\.\d+\.\d+["\']?\s*$',
        f'version: "{version}"', cff, count=1, flags=re.M,
    )
    if count != 1:
        raise SystemExit('Could not update CFF version reliably; no files changed.')
    if old_version != version:
        # A historical release date should not be reused for an unreleased version.
        updated_cff = re.sub(r'^date-released:.*(?:\r?\n|$)', '', updated_cff, flags=re.M)

    readme = readme_path.read_text(encoding='utf-8')
    labels = re.findall(r'\bRelease:\s*v?(\d+\.\d+\.\d+)', readme, flags=re.I)
    if any(label != version for label in labels):
        raise SystemExit(f'Existing README release label {labels} conflicts with VERSION={version}. No changes made.')
    if MARKER not in readme and not labels:
        readme = readme.rstrip() + f'''

{MARKER}
## Publication release and public data

Release: v{version} (pre-release metadata; no GitHub release published yet)

Run `python3 scripts/release_qa.py` (Linux) or `py -3 scripts/release_qa.py` (Windows) before release. Automated checks do **not** certify privacy, scientific reproducibility, permissions, or rights to publish; see `docs/PUBLIC_RELEASE_CHECKLIST.md` and `docs/ZENODO_DATA_GUIDE.md`.

Public analysis code and approved processed data should have separate archival records. Add the Zenodo dataset DOI and NCBI assembly accessions after assignment. Original human-associated raw sequencing reads and restricted participant-level data are not publicly released.
'''
    elif MARKER not in readme:
        # Respect existing metadata; do not duplicate release sections.
        print(f'Existing README release label detected for v{version}; preserving README layout.')
    elif not labels:
        # Recover safely if an older patch added the release section but not
        # a machine-readable Release: label. Only change the marked section.
        readme = readme.replace(
            MARKER,
            f'{MARKER}\nRelease: v{version} (pre-release metadata; no GitHub release published yet)',
            1,
        )

    # Write only after every preflight condition above has passed.
    version_path.write_text(version + '\n', encoding='utf-8')
    if updated_cff != cff:
        cff_path.write_text(updated_cff, encoding='utf-8')
    if readme_path.read_text(encoding='utf-8') != readme:
        readme_path.write_text(readme, encoding='utf-8')
    print(f'PASS: README.md, VERSION and CITATION.cff aligned to v{version}.')
    print('No DOI, GitHub release date, or NCBI accession was invented.')


if __name__ == '__main__':
    main()
