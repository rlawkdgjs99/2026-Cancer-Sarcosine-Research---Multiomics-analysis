#!/usr/bin/env python3
"""Verify the paired display against the frozen 12 original GSEA contrasts.

No new tests or group-wise expression estimates. Run after script 08.
Only the new SC6 PNG is re-encoded, preserving every RGB pixel.
"""
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import io
import json
import re
import struct
import sys

import numpy as np
from PIL import Image, ImageCms
from pypdf import PdfReader
import pypdfium2 as pdfium

np.random.seed(260917)  # No resampling; deterministic verification.
out = Path(__file__).resolve().parents[1]
root = next(p for p in out.parents if (p / 'PROJECT_HANDOFF.md').exists())
qa = out / 'results/qa'
receipt = json.loads((qa / 'SC6_build.json').read_text())
original = list(csv.DictReader((root / receipt['source_csv']).open()))
rows = list(csv.DictReader((root / receipt['display_csv']).open()))
lookup = {(r['compartment'], r['pathway']): r for r in original}
checks = {}


def check(name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


check('12 original contrasts / 24 display cells', len(lookup) == 12 and len(rows) == 24)
check('24 unique positions', len({(r['column_index'], r['row_index']) for r in rows}) == 24)
check('Both orientations for each original contrast', all(
    sum((r['compartment'], r['pathway']) == key for r in rows) == 2 for key in lookup))
for item in receipt['original_inputs'] + receipt['preserved']:
    check('Input preserved ' + item['path'], sha(root / item['path']) == item['sha256'])
for r in rows:
    prefix = r['compartment'] + '/' + r['pathway'] + '/' + r['contrast']
    src = lookup[r['compartment'], r['pathway']]
    sign = int(r['multiplier'])
    check(prefix + '/original NES', float(r['original_NES']) == float(src['NES']))
    check(prefix + '/display NES', float(r['NES']) == sign * float(src['NES']))
    check(prefix + '/correct orientation',
          (r['contrast'], r['numerator'], r['denominator'], sign) in
          [('High vs Low', 'High', 'Low', 1), ('Low vs High', 'Low', 'High', -1)])
    check(prefix + '/paired column order',
          (int(r['column_index']) % 2 == 1) == (sign == 1))
    for field in ['BH_q', 'n', 'High', 'Low', 'cells', 'bh_family_n']:
        check(prefix + '/' + field, float(r[field]) == float(src[field]))
    for field in ['scope', 'eligible_patients', 'bh_family', 'source_q_field', 'source_table']:
        check(prefix + '/' + field, r[field] == src[field])
    check(prefix + '/asterisk', ('*' in r['nes_label']) == (float(r['BH_q']) < .05))
    check(prefix + '/NES label', r['nes_label'].rstrip('*') == f"{float(r['NES']):+.2f}")
for key in lookup:
    a, b = [r for r in rows if (r['compartment'], r['pathway']) == key]
    check(str(key) + '/mirror pair', float(a['NES']) == -float(b['NES']) and a['BH_q'] == b['BH_q'])

png = out / 'results/figures/SC6.png'
pdf = out / 'results/figures/SC6.pdf'
image = Image.open(png)
pixels = np.asarray(image.convert('RGB')).copy()
profile = ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
temp = png.with_name('SC6_encoding_tmp.png')
image.convert('RGB').save(temp, dpi=(600, 600), icc_profile=profile, compress_level=6)
image.close()
check('Encoding preserves every RGB pixel', np.array_equal(pixels, np.asarray(Image.open(temp))))
temp.replace(png)
im = Image.open(png)
plot = receipt['plot']
check('PNG opaque RGB', im.mode == 'RGB')
check('PNG dimensions', im.size == (7920, 4020))
check('PNG 600 dpi', all(abs(v - 600) < .02 for v in im.info['dpi']))
check('PNG sRGB profile', 'srgb' in ImageCms.getProfileDescription(
    ImageCms.ImageCmsProfile(io.BytesIO(im.info['icc_profile']))).lower())
_, _, depth, colour_type, _, _, interlace = struct.unpack('>IIBBBBB', png.read_bytes()[16:29])
check('PNG 8-bit RGB, noninterlaced', (depth, colour_type, interlace) == (8, 2, 0))


def luminance(hex_colour):
    v = np.array([int(hex_colour[k:k+2], 16) / 255 for k in (1, 3, 5)])
    v = np.where(v <= .04045, v / 12.92, ((v + .055) / 1.055)**2.4)
    return v @ np.array([.2126, .7152, .0722])


for r in rows:
    j, i = int(r['column_index']) - 1, int(r['row_index']) - 1
    x = round((plot['x_centres'][j] + .36 * plot['cell_width']) * im.width)
    y = round((1 - plot['y_centres'][i] - .35 * plot['cell_height']) * im.height)
    expected = tuple(int(r['fill'][k:k+2], 16) for k in (1, 3, 5))
    check(f'Cell fill {j}/{i}', im.getpixel((x, y)) == expected)
    check(f'Signed colour for stated contrast {j}/{i}',
          (expected[0] > expected[2]) == (float(r['NES']) > 0))
    lo, hi = sorted([luminance(r['fill']), luminance(r['text_colour'])])
    check(f'Text contrast >= 4.5 {j}/{i}', (hi + .05) / (lo + .05) >= 4.5)

reader = PdfReader(pdf)
check('PDF single page', len(reader.pages) == 1)
page = reader.pages[0]
check('PDF vector without raster image', len(page.images) == 0)
check('PDF page dimensions', abs(float(page.mediabox.width) - 13.2 * 72) < 1 and
      abs(float(page.mediabox.height) - 6.7 * 72) < 1)
text = page.extract_text()
compact = re.sub(r'\s+', '', text).replace('−', '-')
for r in rows:
    key = r['compartment'] + r['pathway'] + r['contrast']
    check('PDF NES ' + key, r['nes_label'] in compact)
    check('PDF q ' + key, r['q_label'].replace(' ', '') in compact)
check('PDF biological labels', all(v in text for v in ['κ', 'γ', 'CD8']))
check('PDF direction labels', all(v in text for v in ['High vs Low', 'Low vs High']))
check('PDF mirrored contrast disclosure', 'mirrored NES, identical q' in text)
check('PDF correct colour interpretation', all(v in text for v in
      ['First group enriched', 'Second group enriched', 'NES for the stated contrast']))
check('GSEA only', 'CAMERA' not in text and 'score_beta' not in text)
doc = pdfium.PdfDocument(str(pdf))
pg = doc[0]
tp = pg.get_textpage()
w, h = pg.get_size()
inside = True
for i in range(tp.count_chars()):
    a, b, c, d = tp.get_charbox(i)
    if c > a and d > b:
        inside &= a >= -1 and b >= -1 and c <= w + 1 and d <= h + 1
check('No text outside PDF page', inside)
pg.render(scale=2).to_pil().save(qa / 'SC6_pdf_preview.png')
tp.close(); pg.close(); doc.close()
for item in receipt['preserved']:
    check('Original remains intact after verification ' + item['path'],
          sha(root / item['path']) == item['sha256'])
artifacts = [png, pdf, root / receipt['display_csv'], out / 'scripts/08_scRNA_two_direction_heatmap.R',
             Path(__file__).resolve()]
result = {
    'checks': checks, 'count': len(checks), 'passed': all(checks.values()),
    'python': sys.version,
    'packages': {p: importlib.metadata.version(p) for p in ['numpy', 'Pillow', 'pypdf', 'pypdfium2']},
    'outputs': {str(p.relative_to(root)): sha(p) for p in artifacts},
    'interpretation': '24 display cells represent 12 original contrasts, with mirrored NES and identical BH q.',
    'limits': 'Visual review recorded separately. No native Finder-to-PPT drag test or author PPT edit.'}
(qa / 'SC6_verification.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
print(f'PASS: {len(checks)} independent checks; SC6 RGB PNG and vector PDF verified.')
