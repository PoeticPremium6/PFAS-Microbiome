#!/usr/bin/env python3
"""Add an idempotent licence summary to PFAS-Microbiome README."""
from pathlib import Path

root = Path(__file__).resolve().parents[1]
readme = root / 'README.md'
if not (root / 'CITATION.cff').is_file() or not (root / 'VERSION').is_file():
    raise SystemExit('ERROR: Not the PFAS-Microbiome publication checkout.')
license_file = root / 'LICENSE'
data_file = root / 'DATA_LICENSE.md'
if not license_file.is_file() or not data_file.is_file():
    raise SystemExit('ERROR: LICENSE or DATA_LICENSE.md is missing.')
if 'GNU GENERAL PUBLIC LICENSE' not in license_file.read_text(encoding='utf-8'):
    raise SystemExit('ERROR: LICENSE is not recognisable as the GNU GPL text.')
if 'CC BY 4.0' not in data_file.read_text(encoding='utf-8'):
    raise SystemExit('ERROR: DATA_LICENSE.md has unexpected content.')
if not readme.is_file():
    raise SystemExit('ERROR: README.md is missing.')
original = readme.read_text(encoding='utf-8')
marker = '<!-- PFAS_DUAL_LICENSE_v1.2.0 -->'
if marker in original:
    print('PASS: Licence summary already present in README (no duplicate added).')
else:
    section = '''<!-- PFAS_DUAL_LICENSE_v1.2.0 -->
## Licensing and reuse

- **Original analysis/figure-generation code:** GNU GPL v3.0 only (`GPL-3.0-only`); see [LICENSE](LICENSE).
- **Original, public-release-approved processed research tables:** CC BY 4.0; see [DATA_LICENSE.md](DATA_LICENSE.md).
- **Third-party resources and restricted human-participant material:** not relicensed by this repository. Refer to [licensing and attribution notes](docs/LICENSING_AND_ATTRIBUTION.md).

Commercial reuse is permitted under the relevant licence. Please credit the creators, cite the associated manuscript and cite the archived software/data DOI once available. The public software and processed-data deposits have separate licences and DOIs.
'''
    readme.write_text(original.rstrip() + '\n\n' + section, encoding='utf-8')
    print('PASS: Appended separate code/data licence scope to README.')
print('NOTE: Review institutional, contributor and third-party permissions before public commit.')
