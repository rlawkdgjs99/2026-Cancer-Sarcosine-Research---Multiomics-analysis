from pathlib import Path
import csv,gzip,json,hashlib,shutil
import numpy as np
from PIL import Image,ImageCms
from pypdf import PdfReader
out=Path(__file__).resolve().parent.parent
b=out.parents[2]
root=next(p for p in out.parents if (p/'PROJECT_HANDOFF.md').exists())
def rows(p):
    with (gzip.open(p,'rt',encoding='utf-8-sig') if p.suffix=='.gz' else p.open(encoding='utf-8-sig')) as f: return list(csv.DictReader(f))
original=rows(b/'results/tables/01_cell_paper_style_scores.csv.gz')
source=rows(out/'tables/source_cells.csv.gz'); plotted=rows(out/'tables/plotted_cells.csv.gz')
assert len(original)==len(source)==92053 and len(plotted)==90512
fields=['umap_1','umap_2','production_module_raw','degradation_module_raw','production_module_shifted','degradation_module_shifted']
maxdiff=0.
for a,c in zip(original,source):
    for k in ['cell_id','Patient','Sample','Resource','final_lineage','analysis_eligible']: assert a[k]==c[k],k
    for k in fields: maxdiff=max(maxdiff,abs(float(a[k])-float(c[k])))
assert maxdiff<1e-12
selected=[x for x in source if x['analysis_eligible'].lower()=='true']
assert plotted==selected
assert len({x['cell_id'] for x in plotted})==90512
assert len({x['Patient'] for x in plotted})==15
caps=rows(out/'tables/display_specification.csv')
for row in caps:
    key=row['axis'].lower()+'_module_shifted'
    cap=float(row['display_cap'])
    q=float(np.quantile([float(x[key]) for x in original],.995))
    assert abs(cap-q)<1e-12
    assert sum(float(x[key])>cap for x in plotted)==int(row['above_cap'])
for item in json.loads((out/'qa/input_manifest.json').read_text()):
    assert hashlib.sha256((root/item['path']).read_bytes()).hexdigest()==item['sha256']
profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
figchecks=[]
for p in sorted((out/'figures/PPT_insert').glob('*.png')):
    im=Image.open(p);im.load();before=im.convert('RGB').tobytes();dpi=im.info.get('dpi')
    assert im.mode=='RGB' and all(abs(v-600)<.1 for v in dpi)
    im.save(p,format='PNG',dpi=(600,600),icc_profile=profile)
    re=Image.open(p);re.load();assert re.tobytes()==before and re.mode=='RGB' and re.info.get('icc_profile')
    data=p.read_bytes();assert data[:8]==b'\x89PNG\r\n\x1a\n' and data[24]==8 and data[25]==2
    pdf=out/'figures'/(p.stem+'.pdf');pages=PdfReader(pdf).pages;assert len(pages)==1
    txt=pages[0].extract_text()
    if p.stem in ['P','DP','M']:assert 'GNMT + DMGDH' in txt
    if p.stem in ['D','DP','M']:assert 'SARDH + PIPOX' in txt
    if p.stem!='L':assert 'Shifted module score' in txt
    figchecks.append({'file':p.name,'pixels':re.size,'mode':re.mode,'dpi':re.info['dpi'],'pixel_preserved_sRGB_export':True,'sha256':hashlib.sha256(data).hexdigest()})
delivery=out/'PPT_ready';delivery.mkdir(exist_ok=True)
for p in list((out/'figures/PPT_insert').glob('*.png'))+list((out/'figures').glob('*.pdf')):
    target=delivery/p.name;shutil.copy2(p,target);assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(target.read_bytes()).digest()
report={'status':'PASS','source_csv_vs_RDS_export_max_abs_difference':maxdiff,'exact_plot_roster':True,'original_cells':92053,'plotted_cells':90512,'patients':15,'input_hashes_unchanged':True,'caps_independently_recomputed':True,'figures':figchecks,'delivery':str(delivery),'native_PowerPoint_drag_tested':False}
(out/'qa/verification.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps(report,ensure_ascii=False,indent=2))
