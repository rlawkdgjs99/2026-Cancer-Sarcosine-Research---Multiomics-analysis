from pathlib import Path
import csv, hashlib, json, subprocess
from PIL import Image, ImageCms
from pypdf import PdfReader
out=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((out/'tables/plot_data.csv').open()))
source=out.parent/'15_AllCell_Degradation_Pathway_Overview_26.09.14/results/tables/Fig6_scRNA_GSEA_heatmap.csv'
original=list(csv.DictReader(source.open(encoding='utf-8-sig')))
for row in rows:
    assert row in original
specs={'A2':(14.5,6.4,None),'W':(8.8,4.6,'Whole_tumour'),'D':(8.8,2.65,'cDC'),'T':(8.8,3.95,'CD8')}
checks={}
for name,(w,h,comp) in specs.items():
    path=out/'figures'/f'{name}.png'
    im=Image.open(path).convert('RGB')
    im.save(path,dpi=(600,600),icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes())
    assert im.size==(round(w*600),round(h*600))
    im.thumbnail((2000,1100));im.save(out/'qa'/f'{name}_preview.png')
    pdfpath=path.with_suffix('.pdf');pdf=PdfReader(pdfpath)
    assert len(pdf.pages)==1
    text=pdf.pages[0].extract_text()
    selected=[r for r in rows if comp is None or r['compartment']==comp]
    for r in selected:
        assert f"{float(r['NES']):.2f}" in text
        q=float(r['BH_q'])
        if q>=.001: assert f"{q:.3f}" in text
    assert 'greaterorequal' not in text and '\x00' not in text
    if name in ['A2','T']: assert 'CD8' in text and '+' in text
    (out/'qa'/f'{name}_pdf_text.txt').write_text(text)
    subprocess.run(['/opt/homebrew/bin/pdftoppm','-scale-to','2000','-singlefile','-png',str(pdfpath),str(out/'qa'/f'{name}_pdf_preview')],check=True)
    checks[name]={'rows':len(selected),'PNG_RGB_sRGB_600dpi':True,'PDF_one_page':True,'NES_labels_checked':True,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
json.dump({'source_rows_unchanged':True,'exports':checks},(out/'qa/large_verification.json').open('w'),indent=2)
print('Four PNG/PDF pairs exported and checked; source data unchanged.')
