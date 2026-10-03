"""Verify original GSEA values, deterministic selection and lossless PNG encoding."""
from pathlib import Path
import csv, hashlib, json, math, struct
from collections import defaultdict
from PIL import Image, ImageCms
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
TAB, QA, FIG = (ROOT / 'results' / s for s in ('tables', 'qa', 'figures'))
def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
checks = []
def check(name, value):
    assert value, name
    checks.append(name)
roots = {'Whole_tumour': ROOT,
         'CD8': BASE / '14_AllCell_Degradation_CD8_cDC1_26.09.08',
         'cDC': BASE / '17_cDC_AllCell_Degradation_26.09.16'}
sources = {}
for ct, root in roots.items():
    rows = [r for r in read(root / 'results/tables/02_GSEA_all_scopes.csv')
            if r['scope'] == 'post_group' and (ct != 'CD8' or r['cell_type'] == 'CD8')]
    sources[ct] = {r['pathway']: r for r in rows}
audit = read(TAB / 'Fig6_HR_all_results_selection_audit.csv')
selected = read(TAB / 'Fig6_HR_selected_plotdata.csv')
sets, universes = defaultdict(set), defaultdict(set)
for r in read(QA / 'HR_membership.csv'):
    sets[r['pathway']].add(r['gene_symbol'])
for r in read(QA / 'HR_ranked_genes.csv'):
    universes[r['compartment']].add(r['gene'])
prior_universes = defaultdict(set)
for r in read(QA / 'GSEA_curve_frozen_ranks.csv'):
    prior_universes[r['compartment']].add(r['gene'])
check('ranked universes equal prior independent frozen-rank export', universes == prior_universes)
for ct in roots:
    original = sources[ct]
    qcol = 'q_global_cDC' if ct == 'cDC' else 'q_global'
    candidates = [r for r in original.values() if r['collection'] in ('Hallmark','Reactome')
                  and float(r['NES']) > 0 and float(r[qcol]) < .05]
    candidates.sort(key=lambda r: (float(r[qcol]), -abs(float(r['NES'])), r['pathway']))
    chosen, decisions = [], {}
    for r in candidates:
        key = r['pathway']; genes = sets[key] & universes[ct]
        ovs = [len(genes & (sets[p] & universes[ct])) / len(genes | (sets[p] & universes[ct])) for p in chosen]
        maximum = max(ovs, default=0)
        state = 'DISPLAY_LIMIT_8' if len(chosen) >= 8 else 'GENE_OVERLAP_GE_0.50' if maximum >= .5 else 'SELECTED'
        decisions[key] = (state, maximum)
        if state == 'SELECTED': chosen.append(key)
    actual = sorted([r for r in selected if r['compartment']==ct],key=lambda r:int(r['display_order']))
    check(f'{ct}: independent exact selected order', [r['pathway'] for r in actual] == chosen)
    for r in [x for x in audit if x['compartment'] == ct]:
        source = original[r['pathway']]
        for field in ('NES','ES','leading_edge_n','n','High','Low'):
            check(f'{ct}/{r["pathway"]}/{field}', math.isclose(float(r[field]),float(source[field]),rel_tol=5e-14,abs_tol=0))
        check(f'{ct}/{r["pathway"]}/BH_q', math.isclose(float(r['BH_q']),float(source[qcol]),rel_tol=5e-14,abs_tol=0))
        check(f'{ct}/{r["pathway"]}/gene_set_size', int(r['gene_set_size'])==int(source['size']))
        check(f'{ct}/{r["pathway"]}/LE_genes', r['leading_edge']==source['leading_edge'])
        check(f'{ct}/{r["pathway"]}/measured_genes', set(r['measured_genes'].split(';')) == sets[r['pathway']] & universes[ct])
        if r['pathway'] in decisions:
            dec, mx = decisions[r['pathway']]
            check(f'{ct}/{r["pathway"]}/selection', r['display_status']==dec and math.isclose(float(r['max_Jaccard_to_selected']),mx,rel_tol=1e-12,abs_tol=1e-14))
        else:
            check(f'{ct}/{r["pathway"]}/not eligible', r['display_status'] != 'SELECTED')
    for r in actual:
        a = next(x for x in audit if x['compartment']==ct and x['pathway']==r['pathway'])
        check(f'{ct}/{r["pathway"]}/plotted row matches audit', all(r[k]==a[k] for k in a))
        check(f'{ct}/{r["pathway"]}/label collection', r['display_label'].endswith('[H]' if r['collection']=='Hallmark' else '[R]'))

profile_path = Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc')
profile = profile_path.read_bytes() if profile_path.exists() else ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
artifacts = []
for stem in ('HW','HT','HD','HR'):
    path = FIG / f'{stem}.png'
    with Image.open(path) as im:
        rgb = im.convert('RGB'); before_pixels = hashlib.sha256(rgb.tobytes()).hexdigest()
        tmp = path.with_suffix('.checked.png')
        rgb.save(tmp,format='PNG',dpi=(600,600),icc_profile=profile,compress_level=6)
    with Image.open(tmp) as im:
        check(f'{stem}: opaque RGB pixels preserved', im.mode=='RGB' and hashlib.sha256(im.tobytes()).hexdigest()==before_pixels)
        check(f'{stem}: 600 dpi and sRGB', abs(im.info['dpi'][0]-600)<.01 and bool(im.info.get('icc_profile')))
        width,height=im.size
        check(f'{stem}: dimensions', (width,height)==((12240,3540) if stem=='HR' else (4080,3540)))
        preview=im.copy();preview.thumbnail((2200,1100));preview.save(QA / f'HR_preview_{stem}.png')
    tmp.replace(path)
    header = path.read_bytes()[:33]
    check(f'{stem}: PNG header RGB8', header[:8]==b'\x89PNG\r\n\x1a\n' and header[24:26]==bytes((8,2)))
    pdf=FIG/f'{stem}.pdf'; reader=PdfReader(pdf)
    check(f'{stem}: PDF one page',len(reader.pages)==1)
    text=reader.pages[0].extract_text()
    check(f'{stem}: PDF content labels', 'Hallmark' in text and 'Reactome' in text and 'Degradation-High' in text and 'Leading-edge' in text)
    check(f'{stem}: no raw p or forbidden legacy method', 'CAMERA' not in text and 'raw p' not in text.lower())
    (QA/f'HR_pdf_text_{stem}.txt').write_text(text)
    artifacts.append({'stem':stem,'width':width,'height':height,'png_sha256':sha(path),'pdf_sha256':sha(pdf),'pixel_sha256':before_pixels})
for r in read(QA/'HR_input_hashes.csv'):
    check(f'protected input {r["path"]}',sha(Path(r['path']))==r['sha256'])
summary={'status':'PASS','number_of_checks':len(checks),'selected_counts':dict((ct,sum(r['compartment']==ct for r in selected)) for ct in roots),'all_GSEA_values_unchanged':True,'selection_independently_reproduced':True,'protected_inputs_unchanged':True,'artifacts':artifacts,'native_powerpoint_drag_tested':False,'visual_review':'Recorded separately in HR_visual_review.md'}
(QA/'HR_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
