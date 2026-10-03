from pathlib import Path
import csv,hashlib,json,math,os,shutil,stat,subprocess,xml.etree.ElementTree as ET,tempfile
import numpy as np
from PIL import Image,ImageCms
base=Path(__file__).resolve().parent.parent
qa=base/'results/qa'; tables=base/'results/tables'; checks=[]
key=qa/'scratch_path.txt'
if key.exists(): tmp=Path(key.read_text().strip())
else:
 tmp=Path(tempfile.mkdtemp(prefix='sarcosine_cdc17_qa_',dir='/private/tmp'));key.write_text(str(tmp))
def ck(name,passed,detail=''):
 checks.append(dict(check=name,pass_check=bool(passed),detail=str(detail)))
 if not passed:raise AssertionError((name,detail))
def read(name):
 with (tables/name).open() as f:return list(csv.DictReader(f))
allg=read('02_GSEA_all_scopes.csv')
for scope in sorted(set(r['scope'] for r in allg)):
 rows=[r for r in allg if r['scope']==scope];p=np.array([float(r['pval']) for r in rows]);order=np.argsort(p);q=np.zeros(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
 ck(scope+' NumPy BH',np.max(np.abs(q-np.array([float(r['q_global_cDC']) for r in rows])))<1e-12)
primary={r['pathway']:r for r in allg if r['scope']=='post_group'}
for name in ['plotdata_GO_display.csv','plotdata_Exact_display.csv']:
 for r in read(name):
  src=primary[r['pathway']]
  for k in ['NES','q_global_cDC','leading_edge_n']:ck(name+' '+r['pathway']+' '+k,float(r[k])==float(src[k]))
profile=Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc').read_bytes(); artifacts=[]
specs=[('cDC_GO_Overview','G',(7440,4260)),('cDC_Exact_Immune_Pathways','I',(5280,3420))]
for stem,short,dims in specs:
 p=base/'results/figures'/f'{stem}.png';im=Image.open(p).convert('RGB');before=hashlib.sha256(im.tobytes()).hexdigest();im.save(p,format='PNG',dpi=(600,600),icc_profile=profile)
 t=p.parent/'PPT_insert'/f'{short}.png';shutil.copyfile(p,t)
 for f in [p,t]:
  subprocess.run(['/usr/bin/xattr','-c',str(f)],check=True)
  subprocess.run([shutil.which('chflags'),'0',str(f)],check=True)
  with Image.open(f) as img:
   ck(f.name+' RGB/dimensions',img.mode=='RGB' and img.size==dims,(img.mode,img.size));ck(f.name+' exact_pixels',hashlib.sha256(img.tobytes()).hexdigest()==before);ck(f.name+' dpi',all(abs(v-600)<.1 for v in img.info['dpi']));ck(f.name+' ICC',img.info.get('icc_profile')==profile)
  header=f.read_bytes()[:33];ck(f.name+' noninterlaced_8bit',header[24]==8 and header[28]==0)
  sip=subprocess.run(['/usr/bin/sips','-g','pixelWidth','-g','pixelHeight','-g','profile',str(f)],capture_output=True,text=True,check=True).stdout
  ck(f.name+' CoreGraphics_decode',str(dims[0]) in sip and str(dims[1]) in sip,sip.strip());ck(f.name+' visible_no_xattrs',not subprocess.run(['/usr/bin/xattr',str(f)],capture_output=True,text=True,check=True).stdout.strip() and (f.stat().st_flags & stat.UF_HIDDEN)==0)
 ck(short+' byte_copy',p.read_bytes()==t.read_bytes())
 svg=p.with_suffix('.svg');root=ET.parse(svg).getroot();ck(stem+' SVG_valid',root.tag.endswith('svg'))
 pdf=p.with_suffix('.pdf');info=subprocess.run(['/opt/homebrew/bin/pdfinfo',str(pdf)],capture_output=True,text=True,check=True).stdout;ck(stem+' PDF_one_page','Pages:           1' in info)
 fonts=subprocess.run(['/opt/homebrew/bin/pdffonts',str(pdf)],capture_output=True,text=True,check=True).stdout;(qa/(stem+'_pdf_fonts.txt')).write_text(fonts)
 subprocess.run(['/opt/homebrew/bin/pdftoppm','-scale-to','2000','-png','-singlefile',str(pdf),str(tmp/(stem+'_pdf'))],check=True)
 subprocess.run(['/opt/homebrew/bin/pdftotext','-bbox-layout',str(pdf),str(tmp/(stem+'_bbox.html'))],check=True)
 doc=ET.parse(tmp/(stem+'_bbox.html'));ns={'h':'http://www.w3.org/1999/xhtml'};pg=doc.find('.//h:page',ns);pw,ph=float(pg.attrib['width']),float(pg.attrib['height']);words=pg.findall('.//h:word',ns);bad=[w.text for w in words if float(w.attrib['xMin'])<-.5 or float(w.attrib['yMin'])<-.5 or float(w.attrib['xMax'])>pw+.5 or float(w.attrib['yMax'])>ph+.5];ck(stem+' PDF_text_within_page',len(bad)==0,bad)
 artifacts.append(dict(stem=stem,transport=str(t),pixel_sha256=before,sha256=hashlib.sha256(t.read_bytes()).hexdigest(),width=dims[0],height=dims[1]))
(qa/'artifact_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));(qa/'transport_manifest.json').write_text(json.dumps(artifacts,indent=2))
print('NUMPY/ARTIFACT CHECKS PASS',len(checks));print(json.dumps(artifacts,indent=2))
