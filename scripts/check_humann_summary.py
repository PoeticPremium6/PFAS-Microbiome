#!/usr/bin/env python3
"""Stream-scan the large HUMAnN association TSV without printing study data.

Checks structure and a subset of metadata-related warning patterns. It is not
an anonymisation audit, scientific validation or clearance for public release.
"""
from __future__ import annotations
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / 'results' / 'machine_readable' / 'Table_S04__HUMAnN_associations.tsv'
EXPECTED = ['Target', 'Feature Type', 'Feature', 'Phenotype', 'Model', 'N',
            'Beta Or Group Effect', 'P value', 'FDR q value', 'Interpretation']
SENSITIVE_COLUMNS = ('patient', 'participant', 'sample_id', 'subject_id',
                     'individual_id', 'person_id', 'birth_date',
                     'postcode', 'latitude', 'longitude')

def main() -> int:
    if not TABLE.is_file():
        print('ERROR: HUMAnN table not found'); return 1
    count = 0
    inconsistent = 0
    invalid_n = 0
    with TABLE.open(encoding='utf-8-sig', newline='', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, [])
        suspicious = [name for name in header if any(x in name.casefold() for x in SENSITIVE_COLUMNS)]
        if header != EXPECTED or suspicious:
            print(f'ERROR: Unexpected column schema; names={header!r}, sensitive={suspicious!r}')
            return 1
        for row in reader:
            if not row or not any(cell.strip() for cell in row):
                continue
            count += 1
            if len(row) != len(EXPECTED):
                inconsistent += 1
            elif row[5].strip() and not row[5].strip().isdigit():
                invalid_n += 1
    print('HUMAnN association table structural audit')
    print(f'Rows examined: {count:,}')
    print(f'Columns: {len(header)}; unexpected-width rows: {inconsistent}; non-integer N values: {invalid_n}')
    print('NOTE: This detects structural anomalies only; it does not establish disclosure safety.')
    if inconsistent:
        print('STATUS: FAIL (row structure)'); return 1
    if invalid_n:
        print('STATUS: REVIEW (N values might be valid missing-value markers)'); return 0
    print('STATUS: PASS (structure only)'); return 0

if __name__ == '__main__':
    sys.exit(main())
