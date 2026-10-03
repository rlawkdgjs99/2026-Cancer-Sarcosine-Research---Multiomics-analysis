#!/usr/bin/env python3
"""Independent numerical and file checks for the Figure 6 e/g/i re-display.

No GSEA, group redefinition or inferential modification. Recompute original full-
family BH only as a diagnostic, never across the displayed twelve results.
Run after 06_scRNA_GSEA_heatmap.R. Paths resolve relative to this script.
"""
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import re
import struct
import sys

import numpy as np
from PIL import Image, ImageCms
from pypdf import PdfReader
import pypdfium2 as pdfium

np.random.seed(260917)  # Verification is deterministic and uses no resampling.
out = Path(__file__).resolve().parents[1]
root = next(p for p in out.parents if (p / 'PROJECT_HANDOFF.md').exists())
qa = out / 'results/qa'
receipt = json.loads((qa / 'SC_build.json').read_text())
rows = list(csv.DictReader((out / 'results/tables/Fig6_scRNA_GSEA_heatmap.csv').open()))
checks = {}


def check(name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1)
    return result


check('Twelve distinct results', len(rows) == 12 and
      len({(r['pathway'], r['compartment']) for r in rows}) == 12)
for item in receipt['inputs']:
    check('Input unchanged: ' + item['path'], sha(root / item['path']) == item['sha256'])

for column, compartment in enumerate(['Whole_tumour', 'CD8', 'cDC']):
    subset = [r for r in rows if r['compartment'] == compartment]
    check(compartment + ' four rows', len(subset) == 4)
    src = list(csv.DictReader((root / subset[0]['source_table']).open()))
    if compartment == 'CD8':
        src = [r for r in src if r['cell_type'] == 'CD8']
    lookup = {r['pathway']: r for r in src}
    full_input = next(item for item in receipt['inputs']
                      if item['path'].endswith('02_GSEA_all_scopes.csv') and
                      ('15_AllCell' if column == 0 else '14_AllCell' if column == 1 else '17_cDC') in item['path'])
    full = [r for r in csv.DictReader((root / full_input['path']).open()) if r['scope'] == 'post_group']
    check(compartment + ' original BH family', len(full) == int(subset[0]['bh_family_n']))
    qfield = subset[0]['source_q_field']
    check(compartment + ' full-family BH reproduced', np.allclose(
        bh([r['pval'] for r in full]), [float(r[qfield]) for r in full], rtol=1e-12, atol=1e-14))
    for row in subset:
        original = lookup[row['pathway']]
        prefix = compartment + '/' + row['pathway']
        for field, source_field in [('NES', 'NES'), ('BH_q', qfield), ('n', 'n'), ('High', 'High'), ('Low', 'Low')]:
            check(prefix + '/' + field, float(row[field]) == float(original[source_field]))
        check(prefix + '/asterisk', ('*' in row['nes_label']) == (float(row['BH_q']) < .05))
        check(prefix + '/sign', row['nes_label'].startswith('+' if float(row['NES']) > 0 else '-'))

png = out / 'results/figures/SC.png'
pdf = out / 'results/figures/SC.pdf'
image = Image.open(png)
pixels = np.asarray(image.convert('RGB')).copy()
# Add explicit sRGB/600-dpi transport metadata without altering any RGB pixel.
profile = ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
temp = png.with_name('SC_encoding_tmp.png')
image.convert('RGB').save(temp, dpi=(600, 600), icc_profile=profile, compress_level=6)
image.close()
check('Encoding preserves all RGB pixels', np.array_equal(pixels, np.asarray(Image.open(temp))))
temp.replace(png)
im = Image.open(png)
check('PNG opaque RGB', im.mode == 'RGB')
check('PNG dimensions', im.size == (6240, 4020))
check('PNG dpi', all(abs(v-600) < .02 for v in im.info['dpi']))
check('PNG sRGB profile', 'srgb' in ImageCms.getProfileDescription(
      ImageCms.ImageCmsProfile(__import__('io').BytesIO(im.info['icc_profile']))).lower())
ihdr = png.read_bytes()[16:29]
_, _, depth, colour_type, _, _, interlace = struct.unpack('>IIBBBBB', ihdr)
check('PNG 8-bit RGB, noninterlaced', (depth, colour_type, interlace) == (8, 2, 0))
for j, compartment in enumerate(['Whole_tumour', 'CD8', 'cDC']):
    subset = [r for r in rows if r['compartment'] == compartment]
    for i, r in enumerate(subset):
        x = round((receipt['plot']['x_centres'][j] + .065) * im.width)
        y = round((1 - receipt['plot']['y_centres'][i] - .038) * im.height)
        expected = tuple(int(r['fill'][k:k+2], 16) for k in (1, 3, 5))
        check(f'Cell fill {j}/{i}', im.getpixel((x, y)) == expected)
        check(f'Signed colour {j}/{i}', (expected[0] > expected[2]) == (float(r['NES']) > 0))
        def luminance(hex_colour):
            v = np.array([int(hex_colour[k:k+2], 16)/255 for k in (1,3,5)])
            v = np.where(v <= .04045, v/12.92, ((v+.055)/1.055)**2.4)
            return v @ np.array([.2126, .7152, .0722])
        lo, hi = sorted([luminance(r['fill']), luminance(r['text_colour'])])
        check(f'Text contrast >= 4.5 {j}/{i}', (hi+.05)/(lo+.05) >= 4.5)

reader = PdfReader(pdf)
check('PDF single page', len(reader.pages) == 1)
page = reader.pages[0]
check('PDF vector, no raster image', len(page.images) == 0)
check('PDF size', abs(float(page.mediabox.width)-10.4*72) < 1 and
      abs(float(page.mediabox.height)-6.7*72) < 1)
text = page.extract_text()
compact = re.sub(r'\s+', '', text).replace('−', '-')
for row in rows:
    check('PDF NES ' + row['compartment'] + row['pathway'], row['nes_label'] in compact)
    check('PDF q ' + row['compartment'] + row['pathway'], row['q_label'].replace(' ', '') in compact)
check('PDF unicode labels', all(s in text for s in ['κ', 'γ', 'CD8']))
check('GSEA only', 'CAMERA' not in text and 'score_beta' not in text)
doc = pdfium.PdfDocument(str(pdf))
pg = doc[0]
tp = pg.get_textpage()
w, h = pg.get_size()
inside = True
for i in range(tp.count_chars()):
    a,b,c,d = tp.get_charbox(i)
    if c > a and d > b:
        inside &= a >= -1 and b >= -1 and c <= w+1 and d <= h+1
check('No text outside PDF page', inside)
pg.render(scale=2).to_pil().save(qa/'SC_pdf_preview.png')
tp.close(); pg.close(); doc.close()
versions = {p: importlib.metadata.version(p) for p in ['numpy','Pillow','pypdf','pypdfium2']}
artifacts = [png,pdf,out/'results/tables/Fig6_scRNA_GSEA_heatmap.csv',
             out/'scripts/06_scRNA_GSEA_heatmap.R',Path(__file__).resolve()]
result = {'checks':checks,'count':len(checks),'passed':all(checks.values()),
          'python':sys.version,'packages':versions,
          'outputs':{str(p.relative_to(root)):sha(p) for p in artifacts},
          'limits':'No author PPT edits or native Finder drag test. Visual review recorded separately.'}
(qa/'SC_verification.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(f'PASS: {len(checks)} independent checks; RGB PNG and vector PDF verified.')
