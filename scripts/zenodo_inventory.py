#!/usr/bin/env python3
"""Inventory *candidate* public files; DOES NOT copy, package or publish data."""
from __future__ import annotations
import argparse
import csv
import hashlib
import pathlib

SUSPECT = ('sample', 'subject', 'participant', 'patient', 'clinical',
           'serum', 'person', 'postcode', 'latitude', 'longitude', 'birth')
SUFFIXES = {'.csv', '.tsv', '.txt', '.json', '.parquet', '.xlsx', '.rds', '.rdata'}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def headers(path):
    if path.suffix.lower() not in {'.csv', '.tsv'}:
        return []
    try:
        with path.open(encoding='utf-8-sig', errors='replace', newline='') as f:
            return next(csv.reader(f, delimiter='\t' if path.suffix.lower() == '.tsv' else ','), [])[:200]
    except (OSError, StopIteration, csv.Error):
        return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=pathlib.Path, default=pathlib.Path('zenodo_candidate_inventory.tsv'))
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / 'results' / 'machine_readable'
    if not source.is_dir():
        parser.error(f'Missing source directory: {source}')
    output = args.out.resolve()
    if output.is_relative_to(source):
        parser.error('Write the inventory outside the source directory')
    with output.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['relative_path', 'bytes', 'sha256', 'columns', 'review_flags', 'release_approval'])
        count = 0
        for path in sorted(source.rglob('*')):
            if not path.is_file() or path.suffix.lower() not in SUFFIXES:
                continue
            fields = headers(path)
            flags = []
            if any(any(s in field.lower() for s in SUSPECT) for field in fields):
                flags.append('sensitive_or_sample_level_column_review')
            if path.suffix.lower() not in {'.csv', '.tsv'}:
                flags.append('manual_content_review')
            if not flags:
                flags.append('manual_disclosure_review_still_required')
            writer.writerow([path.relative_to(root), path.stat().st_size, sha256(path),
                             '|'.join(fields), '|'.join(flags), 'PENDING'])
            count += 1
    print(f'Wrote {count} candidate records to {output}')
    print('NOT APPROVED FOR PUBLICATION: review every file, consent, restrictions and disclosure risks.')

if __name__ == '__main__':
    main()
