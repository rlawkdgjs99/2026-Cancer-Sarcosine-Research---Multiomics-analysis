from pathlib import Path
from PIL import Image,ImageCms
import subprocess,json
b=Path(__file__).resolve().parents[1]
f=b/'figures/H.png'
with Image.open(f) as im:
    im.convert('RGB').save(f,dpi=(600,600),icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes())
with Image.open(f) as im:
    assert im.mode=='RGB' and im.size==(7500,4740) and im.info.get('icc_profile')
    assert abs(im.info['dpi'][0]-600)<.1
subprocess.run(['pdftoppm','-scale-to','1800','-singlefile','-png',str(b/'figures/H.pdf'),str(b/'qa/H_pdf_preview')],check=True,capture_output=True)
r=subprocess.run(['pdftotext','-layout',str(b/'figures/H.pdf'),'-'],capture_output=True,text=True,check=True)
(b/'qa/pdf_text.txt').write_text(r.stdout)
for text in ['Immune-cell composition','Low (n = 5)','High (n = 7)','13 tests','0.775','0.506','0.043','P02','P15']:
    assert text in r.stdout,text
print('PNG format and PDF annotations verified')
