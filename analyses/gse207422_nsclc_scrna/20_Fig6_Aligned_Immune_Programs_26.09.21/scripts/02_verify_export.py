from pathlib import Path
import csv, hashlib, json, subprocess
from PIL import Image, ImageCms
from pypdf import PdfReader
out=Path(__file__).resolve().parents[1]
b=out.parent
source=b/'15_AllCell_Degradation_Pathway_Overview_26.09.14/results/tables/Fig6_scRNA_GSEA_heatmap.csv'
original=list(csv.DictReader(source.open(encoding='utf-8-sig')))
rows=list(csv.DictReader((out/'tables/plot_data.csv').open()))
assert len(rows)==8
for r in rows:
    ref=[v for v in original if (v['pathway'],v['compartment'])==(r['pathway'],r['compartment'])]
    assert len(ref)==1 and ref[0]==r
assert sum(float(r['BH_q'])>=.05 for r in rows)==1
assert {r['bh_family_n'] for r in rows}=={'4958','5959','3797'}
im=Image.open(out/'figures/A.png').convert('RGB')
im.save(out/'figures/A.png',dpi=(600,600),icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes())
check=Image.open(out/'figures/A.png');check.verify()
check=Image.open(out/'figures/A.png')
assert check.mode=='RGB' and check.size==(8700,3000) and check.info.get('icc_profile')
preview=check.copy();preview.thumbnail((2100,1000));preview.save(out/'qa/preview.png')
pdf=PdfReader(out/'figures/A.pdf');assert len(pdf.pages)==1
text=pdf.pages[0].extract_text()
for label in ['2.96','2.56','2.23','3.25','1.73','1.63','1.61','2.65','0.068']:
    assert label in text,label
assert 'greaterorequal' not in text
subprocess.run(['/opt/homebrew/bin/pdftoppm','-scale-to','2100','-singlefile','-png',str(out/'figures/A.pdf'),str(out/'qa/pdf_preview')],check=True)
json.dump({'displayed_rows':8,'source_values_preserved':True,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'full_BH_families':[4958,5959,3797],'nonsignificant_rows':1,'PNG_mode':check.mode,'PNG_size':check.size,'PNG_sRGB':True,'PDF_pages':1,'NES_annotations_verified':True,'files':{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [out/'figures/A.png',out/'figures/A.pdf',out/'tables/plot_data.csv']}},(out/'qa/verification.json').open('w'),indent=2)
print(text)
print('All source, annotation, PDF and PNG checks passed.')
