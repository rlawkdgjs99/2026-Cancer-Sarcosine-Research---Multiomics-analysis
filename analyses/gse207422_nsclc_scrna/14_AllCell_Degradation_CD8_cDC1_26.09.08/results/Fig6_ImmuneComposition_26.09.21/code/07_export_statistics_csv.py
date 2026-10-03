"""Export existing patient data and original model inputs; no new tests."""
from pathlib import Path
import csv
import hashlib
import math

base = Path(__file__).resolve().parents[1]
tables = base / 'tables'
sources = [tables / n for n in ['Patients_ImmuneCell_percent.csv', 'patient_values.csv']]
hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

wide, original = map(read, sources)
cell_types = list(dict.fromkeys(r['cell_type'] for r in original))
lookup = {(r['Patient'], r['cell_type']): r for r in original}
assert len(lookup) == len(original) == 156
assert len(wide) == len({r['Patient'] for r in wide}) == 12
assert sum(r['group'] == 'Low' for r in wide) == 5
assert sum(r['group'] == 'High' for r in wide) == 7
assert {r['histology'] for r in wide} == {'Adeno', 'Squamous'}

meta = ['Patient', 'Sample', 'group', 'High_01', 'histology',
        'Squamous_01', 'degradation_mean_z']
fields = meta + [c + ' (%)' for c in cell_types] + [c + ' (log_ratio)' for c in cell_types]
out = []
for r in wide:
    result = {k: r[k] for k in meta if k in r}
    result['High_01'] = int(r['group'] == 'High')
    result['Squamous_01'] = int(r['histology'] == 'Squamous')
    for ct in cell_types:
        src = lookup[(r['Patient'], ct)]
        assert all(src[k] == r[k] for k in ['Patient', 'Sample', 'group', 'histology', 'degradation_mean_z'])
        assert r[ct + ' (%)'] == src['percent']
        pct = float(src['percent'])
        assert 0 < pct < 100
        assert math.isclose(math.log(pct / (100 - pct)), float(src['log_ratio']), abs_tol=1e-12)
        result[ct + ' (%)'] = src['percent']
        result[ct + ' (log_ratio)'] = src['log_ratio']
    out.append(result)

target = tables / 'Patients_ImmuneCell_PRISM_statistics.csv'
with target.open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(out)

saved = read(target)
assert len(saved) == 12 and len(saved[0]) == 33
for actual, expected in zip(saved, out):
    assert all(actual[k] == str(expected[k]) for k in fields)
assert all(hashlib.sha256(p.read_bytes()).hexdigest() == h for p, h in hashes.items())
print(f'Saved: {target}')
print('Verified: 12 patients (Low5/High7), 33 columns, 156 percentages and 156 original log-ratios.')
print('Original files unchanged. No statistical tests rerun; no native PRISM validation performed.')
