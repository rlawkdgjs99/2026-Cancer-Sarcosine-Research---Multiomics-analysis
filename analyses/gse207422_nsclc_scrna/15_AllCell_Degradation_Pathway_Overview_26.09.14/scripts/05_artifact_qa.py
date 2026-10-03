from pathlib import Path
import csv, hashlib, json, math, os, shutil, stat, subprocess, xml.etree.ElementTree as ET
import numpy as np
from PIL import Image, ImageCms
base=Path(__file__).resolve().parent.parent
tmp=Path('/private/tmp/sarcosine_fig8e_s4mgvcl_');tmp.mkdir(exist_ok=True)
tables=base/'results/tables'; qa=base/'results/qa'; checks=[]
def read(name):
 with (tables/name).open() as f:return list(csv.DictReader(f))
def ck(name,passed,detail=''):
 checks.append(dict(check=name,pass_check=bool(passed),detail=str(detail)))
 if not passed:raise AssertionError((name,detail))
def twop(t,df):
 t=abs(t);n=8192;x=np.linspace(0,t,n+1);a=math.gamma((df+1)/2)/(math.sqrt(df*math.pi)*math.gamma(df/2));v=a*(1+x*x/df)**(-(df+1)/2);inte=(t/n)/3*(v[0]+v[-1]+4*v[1:-1:2].sum()+2*v[2:-1:2].sum());return 1-2*inte
crit={}
def critical(df):
 if df not in crit:
  lo,hi=0.,20.
  for _ in range(55):
   mid=(lo+hi)/2
   if twop(mid,df)>.05:lo=mid
   else:hi=mid
  crit[df]=(lo+hi)/2
 return crit[df]
sv=read('04_patient_score_values.csv');ss=read('04_exact_scores_HC3.csv')
for sc in sorted({r['scope'] for r in ss}):
 design=read('verification_design_'+sc+'.csv');cols=[k for k in design[0] if k!='Patient'];X=np.array([[float(r[k]) for k in cols] for r in design]);j=cols.index('exposure_z' if sc.endswith('continuous') else 'groupHigh'); inv=np.linalg.inv(X.T@X);h=np.sum((X@inv)*X,axis=1);df=len(X)-len(cols)
 for result in [r for r in ss if r['scope']==sc]:
  id=result['pathway'];vm={r['Patient']:float(r['score']) for r in sv if r['scope']==sc and r['pathway']==id};y=np.array([vm[r['Patient']] for r in design]);bt=inv@X.T@y;res=y-X@bt;u=res/(1-h);cov=inv@((X*u[:,None]).T@(X*u[:,None]))@inv;se=math.sqrt(cov[j,j]);p=twop(bt[j]/se,df);qcrit=critical(df)
  vals={'beta':bt[j],'se':se,'p':p,'ci_low':bt[j]-qcrit*se,'ci_high':bt[j]+qcrit*se}
  for k,value in vals.items():ck(sc+' '+id+' independent_numpy_'+k,abs(value-float(result[k]))<1e-9,abs(value-float(result[k])))
gsea={(r['scope'],r['pathway']):r for r in read('02_GSEA_all_scopes.csv')};plot=read('plotdata_Fig8e_GO_overview.csv')
ck('GO_display_16_rows',len(plot)==16)
for r in plot:
 q=gsea[('post_group',r['pathway'])]
 for key in ['NES','q_global','leading_edge_n']:ck('plot_source '+r['pathway']+' '+key,float(r[key])==float(q[key]))
 ck('direction '+r['pathway'],(float(r['NES'])>0)==(r['direction']=='High'))
 ck('significance '+r['pathway'],float(r['q_global'])<.05)
profile=Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc').read_bytes(); artifacts=[]
for stem,short,dims in [('Fig8e_AllCell_Degradation_GO_Overview','PPT_INSERT_Fig8e_GO',(7440,4260)),('Fig8e_Exact_Fig7_Pathway_Checks','PPT_INSERT_Fig8e_Pathways',(6360,5100))]:
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
