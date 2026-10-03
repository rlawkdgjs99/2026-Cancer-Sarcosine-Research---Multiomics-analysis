#!/usr/bin/env python3
"""Independently check output membership, displayed colored cells, and PDF labels."""
from pathlib import Path
import csv, hashlib, json, unicodedata, shutil
from PIL import Image
import pypdfium2 as pdfium
BASE=Path(__file__).resolve().parent
CO=['PRJNA751792','PRJNA1023797','PRJEB22863']
PAIRS={'a':['Production-associated','NR-enriched'], 'b':['Degradation-associated','R-enriched']}
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig')))
source_path=BASE/'csv/Venn_3cohorts_membership.csv'
source=read(source_path)
source_lookup={(r['Category'],r['Species_full']):r for r in source}
output_path=BASE/'csv/Species_4criteria_cohort_membership_plot.csv'
rows=read(output_path)
assert len(rows)==102
assert len({(r['Panel'],r['Species_full'],r['Category'],r['Cohort']) for r in rows})==len(rows)
summary={}
for panel,cats in PAIRS.items():
 expected={r['Species_full'] for r in source if r['Category'] in cats and any(r[c]=='TRUE' for c in CO)}
 actual={r['Species_full'] for r in rows if r['Panel']==panel}
 assert expected==actual
 summary[panel]={'n_species':len(actual),'n_squares':len([r for r in rows if r['Panel']==panel]),'criteria':{}}
 for cat in cats:
  counts={co:sum(r[co]=='TRUE' for r in source if r['Category']==cat) for co in CO}
  summary[panel]['criteria'][cat]=counts
for r in rows:
 src=source_lookup[(r['Category'],r['Species_full'])]
 assert r['Category']==PAIRS[r['Panel']][int(r['Criterion_order'])-1]
 assert r['Cohort']==CO[int(r['Cohort_order'])-1]
 assert r['Species']==src['Species'] and r['Species_label']==src['Species'].replace('_',' ')
 assert r['Meets_criterion']==src[r['Cohort']]
 assert r['Square_fill']==('#EB9438' if src[r['Cohort']]=='TRUE' else '#FFFFFF')
# Independently reconstruct the complete Venn memberships from frozen statistics.
import math
TEMP=Path('/private/tmp/nsclc_panel_deck_review_260909')
TEMP.mkdir(parents=True,exist_ok=True)
def num(s):
 try:return float(s)
 except (TypeError,ValueError):return math.nan
categories=[x for pair in PAIRS.values() for x in pair]
rebuilt={(cat,c):set() for cat in categories for c in CO}
scorepath=BASE/'csv/Species_pathway_score_Spearman_all_tests.csv'
dapath=BASE/'csv/Species_DA_all_eligible.csv'
for r in read(scorepath):
 if r['Analysis'] in CO and num(r['rho'])>.3 and num(r['BH_q'])<.05:
  rebuilt[(r['Score']+'-associated',r['Analysis'])].add(r['Species_full'])
for r in read(dapath):
 if r['Analysis'] in CO and num(r['BH_q_global'])<.05:
  if num(r['log2FC_NR_vs_R'])>1:rebuilt[('NR-enriched',r['Analysis'])].add(r['Species_full'])
  elif num(r['log2FC_NR_vs_R'])< -1:rebuilt[('R-enriched',r['Analysis'])].add(r['Species_full'])
for r in source:
 for c in CO:assert (r[c]=='TRUE') == (r['Species_full'] in rebuilt[(r['Category'],c)])
input_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_path,scorepath,dapath]}
protected=json.loads((BASE/'qa/cohort_membership_protected_before.json').read_text())
for filename,digest in protected.items():
 assert hashlib.sha256((BASE/filename).read_bytes()).hexdigest()==digest
layout=[('Main_species_4criteria_cohort_membership',['a','b'],10.50),
        ('Main_species_Production_NR_cohort_membership',['a'],4.34),
        ('Main_species_Degradation_R_cohort_membership',['b'],7.64)]
image_checks=0; pdf_checks=0; files={}; renders=[]
for name,panels,height in layout:
 png=BASE/'figures'/f'{name}.png'; pdfpath=BASE/'figures'/f'{name}.pdf'
 im=Image.open(png).convert('RGB');assert im.size==(2805,round(height*300))
 pdf=pdfium.PdfDocument(pdfpath);assert len(pdf)==1
 page=pdf[0];widthpt,heightpt=page.get_size();assert abs(widthpt-9.35*72)<1 and abs(heightpt-height*72)<1
 textpage=page.get_textpage();txt=textpage.get_text_range()
 for panel in panels:
  for species in {r['Species_label'] for r in rows if r['Panel']==panel}:assert species in txt, species
  for cohort in CO:assert cohort in txt
 for ci in range(textpage.count_chars()):
  box=textpage.get_charbox(ci)
  assert box[0]>=-.5 and box[1]>=-.5 and box[2]<=widthpt+.5 and box[3]<=heightpt+.5,box
 pim=page.render(scale=300/72).to_pil().convert('RGB')
 ytop=.88
 for panel in panels:
  for r in [r for r in rows if r['Panel']==panel]:
   x=[4.39,6.77][int(r['Criterion_order'])-1]+[.34,.96,1.58][int(r['Cohort_order'])-1]
   y=ytop+1.42+(int(r['Row_order'])-.5)*.30
   expect=(235,148,56) if r['Meets_criterion']=='TRUE' else (255,255,255)
   for pic in [im,pim]:
    px=pic.getpixel((round(x*300),round(y*300)))
    assert max(abs(a-b) for a,b in zip(px,expect))<=2,(name,r,px,expect)
   image_checks+=1;pdf_checks+=1
  ytop+=1.42+summary[panel]['n_species']*.30+.06+.48
 render=TEMP/f'{name}_PDF.png';pim.save(render);renders.append(str(render))
 for path in [png,pdfpath]:files[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
 textpage.close();page.close();pdf.close()
# Keep the minimal plot CSV beside the previously delivered NSCLC plotting CSVs.
NS=BASE.parents[2]
used=next(p for p in NS.iterdir() if unicodedata.normalize('NFC',p.name)=='사용데이터_모음')
delivery=used/'NSCLC_7KO_CRC_style_panels_26.09.09'/output_path.name
shutil.copyfile(output_path,delivery)
assert delivery.read_bytes()==output_path.read_bytes()
for p in [output_path,BASE/'plot_cohort_membership.R',Path(__file__)]:files[str(p.relative_to(BASE))]=hashlib.sha256(p.read_bytes()).hexdigest()
report={'source_sha256':input_sha256,'selection':'Union of either paired criterion across >=1 cohort; no >=2-cohort gate','summary':summary,'csv_cells_verified':102,'PNG_cell_centers_verified':image_checks,'PDF_cell_centers_verified':pdf_checks,'PDF_text_and_bounds':'all three one-page PDFs passed','prior_files_unchanged':len(protected),'plot_csv_delivery':str(delivery),'output_sha256':files}
(BASE/'qa/cohort_membership_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps({'summary':summary,'PNG_checks':image_checks,'PDF_checks':pdf_checks,'protected_unchanged':len(protected),'CSV_copy':str(delivery),'renders':renders},ensure_ascii=False,indent=2))
