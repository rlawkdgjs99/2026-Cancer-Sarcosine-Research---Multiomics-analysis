"""Independent CSV reader checks for the six-KO input export. No statistical fit."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import math
import os
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
base = next(root.glob('*/HGMT_CRC_WGS-Healthy_vs_Cancer'))
analysis = base / 'Healthy_vs_Cancer_4_CRC_cohorts_integrated'
audit = analysis / 'results_integrated/CRC_WGS_6KO_input_export_26.09.07'
out = next(base.glob('*/CRC_WGS_6KO_4cohort_analysis_inputs_26.09.07'))
usage = out.parent
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        records = list(reader)
        assert all(None not in r and None not in r.values() for r in records)
        return reader.fieldnames, records

manifest = json.loads((audit / 'export_manifest.json').read_text())
for item in manifest['sources'] + manifest['outputs']:
    assert sha(root / item['path']) == item['sha256'], item['path']
_, source = rows(analysis / 'results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv')
header, data = rows(out / 'CRC_WGS_6KO_per_sample.csv')
mh, metadata = rows(out / 'CRC_WGS_6KO_sample_metadata.csv')
kos = [x[0] for x in manifest['ko_names']]
assert header == ['Run.ID', 'Cohort', 'Group'] + kos
assert len(data) == len(source) == len(metadata) == 1647
assert len(mh) == 27
assert len({x['Run.ID'] for x in data}) == 1647
assert len({x['Run ID'] for x in metadata}) == 1647
for row, src, meta in zip(data, source, metadata):
    assert row == {k: src[k] for k in header}  # Exact numeric text, not a tolerance test.
    assert row['Run.ID'] == meta['Run ID']
    assert row['Cohort'] == meta['Cohort']
    assert row['Group'] == meta['Analysis Group']
    assert all(math.isfinite(float(row[k])) and float(row[k]) >= 0 for k in kos)
assert sha(out / 'CRC_WGS_6KO_sample_metadata.csv') == sha(usage / 'CRC_WGS_4cohort_pooled_26.09.01/CRC_WGS_4cohort_pooled_metadata.csv')
assert len(list(out.glob('*.csv'))) + len(list((out / 'Forest_plot').glob('*.csv'))) == 8
for p in (out / 'Forest_plot').glob('*.csv'):
    assert sha(p) == sha(usage / 'CRC_WGS_7KO_4cohort_forest_source_data_26.09.01' / p.name)
    h, rr = rows(p)
    assert len(h) == 10 and len(rr) == 4
observed = collections.Counter((r['Cohort'],r['Group']) for r in data)
assert observed == {('PRJEB6070','Healthy'):476,('PRJEB6070','Cancer'):590,
                    ('PRJEB10878','Healthy'):54,('PRJEB10878','Cancer'):74,
                    ('PRJEB27928','Healthy'):120,('PRJEB27928','Cancer'):140,
                    ('PRJNA429097','Healthy'):95,('PRJNA429097','Cancer'):98}
for ko, cohort in [('K00301','PRJEB10878'),('K00306','PRJEB27928'),('K00315','PRJEB27928')]:
    assert all(float(r[ko]) == 0 for r in data if r['Cohort'] == cohort)
for p in [*out.rglob('*'), audit / 'KO_data_preview.png']:
    if p.is_file():
        assert not p.name.startswith('.')
        # APFS may report UF_TRACKED (0x40); it is not a hidden-file flag.
        assert not (getattr(p.stat(), 'st_flags', 0) & 0x8000), p
        xa = subprocess.run(['/usr/bin/xattr','-px','com.apple.FinderInfo',str(p)], capture_output=True, text=True)
        if xa.returncode:
            assert 'No such xattr' in xa.stderr, xa.stderr
            info = b''
        else:
            info = bytes.fromhex(xa.stdout)
        assert not (len(info) >= 10 and int.from_bytes(info[8:10], 'big') & 0x4000), p
report = {
    'status':'PASS', 'python':sys.version, 'rows':1647,
    'KO_cells_exact_string_match':1647*6, 'KO_columns':kos,
    'metadata_byte_identical':True, 'sample_ids_and_row_order_match':True,
    'cohort_group_labels_match':True, 'six_forest_CSVs_byte_identical':True,
    'inputs_and_outputs_match_manifest':True, 'file_visibility_verified':True,
    'counts':{f'{c}/{g}':n for (c,g),n in sorted(observed.items())},
    'zeros_by_KO':{k:sum(float(r[k]) == 0 for r in data) for k in kos},
    'missing_KO_values':0, 'statistics_recomputed':False,
}
(audit / 'verification_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
