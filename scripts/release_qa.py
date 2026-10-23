#!/usr/bin/env python3
"""Publication preflight: metadata consistency and *heuristic* disclosure checks.

Not a privacy certification, licensing review, or scientific reproduction test.
Uses Python 3 standard library only; does not inspect restricted raw data.
"""
from __future__ import annotations
import argparse
import pathlib
import re
import subprocess
import sys

RAW_ENDINGS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bam', '.cram',
               '.sam', '.fasta', '.fasta.gz', '.fa', '.fa.gz', '.fna', '.fna.gz',
               '.vcf', '.vcf.gz', '.bcf', '.sra')
TEXT_ENDINGS = ('.py', '.r', '.rmd', '.sh', '.md', '.txt', '.yml', '.yaml',
                '.json', '.cff', '.tsv', '.csv', '.toml', '.ini', '.R')
SECRET_PATTERNS = [
    ('possible GitHub token', re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{25,}|github_pat_[A-Za-z0-9_]{35,})\b')),
    ('possible AWS access key', re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b')),
    ('private key delimiter', re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----')),
]
PERSONAL_PATH = re.compile(r'(?:/home/[A-Za-z0-9_.-]+/|[A-Za-z]:\\\\?Users\\[A-Za-z0-9_. -]+\\)')
SAMPLE_HEADERS = ('participant', 'patient_id', 'subject_id', 'sample_id',
                  'individual_id', 'person_id', 'clinical', 'date_of_birth',
                  'birth_date', 'postcode', 'postal_code', 'latitude', 'longitude',
                  'pfas_serum', 'pfas_concentration')


def git_lines(root: pathlib.Path, args: list[str]) -> list[str] | None:
    try:
        res = subprocess.run(['git', '-C', str(root), *args], check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             encoding='utf-8', errors='replace')
        return [x for x in res.stdout.splitlines() if x.strip()]
    except (OSError, subprocess.CalledProcessError):
        return None


def tracked_paths(root: pathlib.Path) -> list[pathlib.Path]:
    # Include staged and untracked (non-ignored) files locally so a freshly
    # extracted patch is checked before the first commit, not only afterward.
    names = git_lines(root, ['ls-files', '--cached', '--others', '--exclude-standard'])
    if names is None:
        return [p.relative_to(root) for p in root.rglob('*') if p.is_file()
                and '.git' not in p.relative_to(root).parts]
    return sorted({pathlib.Path(name) for name in names})


def parse_version(root: pathlib.Path):
    version = (root / 'VERSION').read_text(encoding='utf-8').strip() if (root / 'VERSION').exists() else ''
    cff = (root / 'CITATION.cff').read_text(encoding='utf-8') if (root / 'CITATION.cff').exists() else ''
    readme = (root / 'README.md').read_text(encoding='utf-8') if (root / 'README.md').exists() else ''
    cff_match = re.search(r'^version:\s*["\']?([0-9]+\.[0-9]+\.[0-9]+)', cff, re.M)
    readme_match = re.search(r'\bRelease:\s*v?([0-9]+\.[0-9]+\.[0-9]+)', readme, re.I)
    return version.removeprefix('v'), cff_match.group(1) if cff_match else '', readme_match.group(1) if readme_match else ''


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    errors: list[str] = []
    warnings: list[str] = []
    version, citation, readme = parse_version(root)
    if not version or version != citation or version != readme:
        errors.append(f'Metadata mismatch: VERSION={version!r} CITATION.cff={citation!r} README.md={readme!r}')
    for expected in ['README.md', 'CITATION.cff', 'data/ACCESSIONS.tsv',
                     'docs/REPRODUCIBILITY.md', 'workflow/README.md',
                     'workflow/final_figures', 'results/machine_readable']:
        if not (root / expected).exists():
            errors.append(f'Missing expected publication artifact: {expected}')
    if not (root / 'LICENSE').exists() and not (root / 'LICENSE.md').exists():
        warnings.append('No LICENSE file. Confirm contributors/institutions allow the chosen software license.')

    paths = tracked_paths(root)
    for relative in paths:
        path = root / relative
        if not path.is_file():
            continue
        label = relative.as_posix()
        lower = label.lower()
        if lower.endswith(RAW_ENDINGS):
            errors.append(f'Potential raw/genomic sequence artifact committed: {label}')
        if path.stat().st_size > 50 * 1024 * 1024:
            warnings.append(f'Large tracked file >50 MiB (review release size/content): {label}')
        if lower.endswith(TEXT_ENDINGS) or relative.name in ('VERSION', '.gitignore', '.gitattributes'):
            if path.stat().st_size > 8 * 1024 * 1024:
                warnings.append(f'Text file too large for lightweight content scan: {label}')
                continue
            raw = path.read_text(encoding='utf-8', errors='replace')
            for name, pattern in SECRET_PATTERNS:
                if pattern.search(raw):
                    errors.append(f'Potential {name} in {label}; inspect privately')
            if PERSONAL_PATH.search(raw):
                warnings.append(f'Absolute personal filesystem path may be present: {label}')
            if relative.suffix.lower() in {'.csv', '.tsv'}:
                head = raw.splitlines()[0].lower() if raw else ''
                found = [h for h in SAMPLE_HEADERS if h in head]
                if found:
                    warnings.append(f'Manual participant/sample-level disclosure review required: {label} ({", ".join(found)})')

    history = git_lines(root, ['log', '--all', '--pretty=format:', '--name-only'])
    if history:
        previously_named = sorted({p for p in history if p.lower().endswith(RAW_ENDINGS)})
        if previously_named:
            warnings.append(f'{len(previously_named)} possible restricted sequence file path(s) appear in Git history. Review history separately before public release.')

    print(f'PFAS publication preflight: {root}')
    print(f'Versioned + untracked non-ignored files examined: {len(paths)}')
    for item in sorted(set(errors)):
        print('ERROR: ' + item)
    for item in sorted(set(warnings)):
        print('REVIEW: ' + item)
    print('NOTE: This scan cannot establish privacy safety, validate results, or certify Git history.')
    print(f'STATUS: {"FAIL" if errors else "PASS WITH MANUAL REVIEW REQUIRED"}')
    return 1 if errors else 0

if __name__ == '__main__':
    sys.exit(main())
