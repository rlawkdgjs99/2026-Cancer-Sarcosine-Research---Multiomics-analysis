#!/usr/bin/env python3
"""Independently reconstruct proportions, stratified permutation tests, and PNG transport."""
from pathlib import Path
import json,hashlib,itertools,subprocess,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
import pypdfium2 as pdfium
base=Path(__file__).resolve().parents[1];tab=base/'results/tables';fig=base/'results/figures'
tmp=Path('/private/tmp/fig8_cellfrac_qa');tmp.mkdir(exist_ok=True)
plan=json.loads((tab/'Fig8_CellFrac_plan.json').read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=[]
def check(name,condition):
    assert condition,name
    checks.append(name)
for name,x in plan['inputs'].items():check('input unchanged '+name,sha(x['path'])==x['sha256'])
md=pd.read_csv(plan['inputs']['patients']['path']);post=md[md.treatment=='Post'].sort_values('Patient').reset_index(drop=True)
check('patient roster',len(post)==12 and post.group.value_counts().to_dict()=={'High':7,'Low':5})
check('original labels',list(np.where(md.degradation_mean_z>md.cutoff,'High','Low'))==list(md.group))
cells=pd.read_csv(plan['inputs']['cells']['path'],usecols=['cell_id','Patient','Sample','final_lineage'])
cells=cells[(cells.final_lineage!='Excluded residual doublet')&cells.Patient.isin(post.Patient)]
check('78192 unique retained cells',len(cells)==78192 and not cells.cell_id.duplicated().any())
check('patient totals',cells.groupby('Patient').size().to_dict()==post.set_index('Patient').all_cells.to_dict())
comp=pd.read_csv(plan['inputs']['composition']['path']);actual=cells.groupby(['Patient','final_lineage']).size()
expected=comp[comp.Patient.isin(post.Patient)].set_index(['Patient','final_lineage']).cells
check('all lineage counts',actual.to_dict()==expected.to_dict())
check('target totals',int((cells.final_lineage=='CD8 T cell').sum())==16229 and int((cells.final_lineage=='Conventional DC').sum())==901)
pl=pd.read_csv(tab/'Fig8_CellFrac_patients.csv');stats=pd.read_csv(tab/'Fig8_CellFrac_statistics.csv');summ=pd.read_csv(tab/'Fig8_CellFrac_summary.csv');pt=pd.read_csv(tab/'Fig8_CellFrac_permutations.csv')
ids=list(post.Patient);h=(post.group=='High').astype(float).to_numpy();hist=(post.histology=='Squamous').astype(float).to_numpy()
strata=[np.flatnonzero(hist==v) for v in [0,1]]
all_h=[]
for a in itertools.combinations(strata[0],int(h[strata[0]].sum())):
 for b in itertools.combinations(strata[1],int(h[strata[1]].sum())):
  u=np.zeros(12);u[list(a)+list(b)]=1;all_h.append(u)
all_h=np.asarray(all_h);check('300 unique stratified assignments',all_h.shape==(300,12) and len(np.unique(all_h,axis=0))==300)
def fit(y,g):
    X=np.column_stack([np.ones(12),hist,g]);coef=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@coef;cov=np.linalg.inv(X.T@X)*(resid@resid/9)
    return coef[2],coef[2]/np.sqrt(cov[2,2])
def bh(p):
    p=np.asarray(p);o=np.argsort(p);q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return q
all_num=[]
for scope in ['all_singlets','immune_singlets']:
 raw=cells if scope=='all_singlets' else cells[~cells.final_lineage.isin(['Epithelial','CAF'])]
 den=raw.groupby('Patient').size().reindex(ids).to_numpy()
 if scope=='immune_singlets':check('immune total',den.sum()==68093)
 ps=[]
 for cell in ['CD8 T cell','Conventional DC']:
  num=raw[raw.final_lineage==cell].groupby('Patient').size().reindex(ids,fill_value=0).to_numpy()
  pct=100*num/den;y=np.log(num/(den-num));b,t=fit(y,h);ts=np.array([fit(y,u)[1] for u in all_h]);p=np.mean(np.abs(ts)>=abs(t)-1e-12);ps.append(p)
  r=stats[(stats.scope==scope)&(stats.cell_type==cell)].iloc[0]
  d=pl[(pl.scope==scope)&(pl.cell_type==cell)].set_index('Patient').loc[ids]
  check(scope+' '+cell+' plot numerator/denominator',np.array_equal(num,d.target_cells) and np.array_equal(den,d.denominator_cells))
  np.testing.assert_allclose(pct,d.percent,rtol=0,atol=1e-11);np.testing.assert_allclose(y,d.log_ratio,rtol=0,atol=1e-12)
  np.testing.assert_allclose([b,t,p],[r.log_ratio_beta_High_minus_Low,r.t_statistic,r.raw_P],rtol=0,atol=1e-11)
  perm=pt[(pt.scope==scope)&(pt.cell_type==cell)]
  np.testing.assert_allclose(ts,perm.t_statistic,rtol=0,atol=1e-11)
  check(scope+' '+cell+' independent full permutation',int(r.extreme_assignments)==int(round(p*300)))
  for group in ['Low','High']:
   v=pct[(post.group==group).to_numpy()];s=summ[(summ.scope==scope)&(summ.cell_type==cell)&(summ.group==group)].iloc[0]
   np.testing.assert_allclose([v.mean(),v.std(ddof=1),np.median(v)],[s.mean_percent,s.SD_percent,s.median_percent],rtol=0,atol=1e-11)
  if scope=='all_singlets':
   prism=pd.read_csv(tab/('Fig8_'+('CD8' if cell=='CD8 T cell' else 'cDC')+'_Prism.csv'))
   for group in ['Low','High']:
    vals=prism[[group+'_Patient',group+'_percent']].dropna().set_index(group+'_Patient').iloc[:,0]
    for patient,value in vals.items():check('Prism '+cell+' '+patient,abs(value-d.loc[patient,'percent'])<1e-11)
 np.testing.assert_allclose(bh(ps),stats[stats.scope==scope].BH_q,rtol=0,atol=1e-12)
 check(scope+' independent BH2',True)
# Lossless scientific-figure PNG packaging. No drawing/content editing here.
icc=Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc').read_bytes()
for f in [fig/'Fig8_CD8_cDC_Proportions.png',fig/'C8.png',fig/'DC.png']:
 with Image.open(f) as im:
  before=np.array(im.convert('RGB'));im.convert('RGB').save(f,dpi=(600,600),icc_profile=icc,compress_level=6)
 with Image.open(f) as im:check(f.name+' reencoding pixel identity',np.array_equal(before,np.array(im)))
shutil.copyfile(fig/'Fig8_CD8_cDC_Proportions.png',fig/'CD.png')
outputs={}
for stem in ['CD','C8','DC']:
 f=fig/(stem+'.png');pdf=fig/(stem+'.pdf')
 for path in [f,pdf]:
  subprocess.run(['/usr/bin/xattr','-c',str(path)],check=True);subprocess.run(['/usr/bin/chflags','nohidden',str(path)],check=True)
  check(path.name+' xattrs clear',not subprocess.check_output(['/usr/bin/xattr',str(path)]).strip())
  check(path.name+' visible',not(path.stat().st_flags & 0x8000))
 with Image.open(f) as im:
  im.load();pixels=np.asarray(im);size=im.size
  check(stem+' RGB8 DPI ICC',im.mode=='RGB' and all(abs(v-600)<.02 for v in im.info['dpi']) and im.info['icc_profile']==icc)
  check(stem+' PNG structure',f.read_bytes()[24:26]==bytes([8,2]) and f.read_bytes()[28]==0)
  check(stem+' expected dimensions',size==((4080,2490) if stem=='CD' else (2070,2100)))
 raw=tmp/(stem+'.rgb')
 subprocess.run(['/usr/local/bin/Rscript','-e',"a<-commandArgs(TRUE);x<-png::readPNG(a[1]);writeBin(as.raw(round(as.vector(aperm(x,c(3,2,1)))*255)),a[2])",str(f),str(raw)],check=True,cwd=str(tmp))
 check(stem+' Pillow/libpng identity',raw.read_bytes()==pixels.tobytes())
 reader=PdfReader(pdf);check(stem+' one PDF page',len(reader.pages)==1);txt=reader.pages[0].extract_text();(tmp/(stem+'_pdf.txt')).write_text(txt)
 for word in (['0.537','0.467'] if stem=='CD' else ['0.537'] if stem=='C8' else ['0.467']):check(stem+' PDF q '+word,word in txt)
 check(stem+' PDF no raw P','raw P' not in txt)
 doc=pdfium.PdfDocument(str(pdf));doc[0].render(scale=2).to_pil().save(tmp/(stem+'_pdf.png'))
 outputs[stem]={'png':str(f),'pdf':str(pdf),'size':size,'sha256':sha(f),'file':subprocess.check_output(['/usr/bin/file',str(f)],text=True).strip()}
check('scientific and short transport identical',sha(fig/'Fig8_CD8_cDC_Proportions.png')==sha(fig/'CD.png'))
old=Path('/private/tmp/fig8_cellfrac_protected.json')
if old.exists():
 for p,h0 in json.loads(old.read_text()).items():check('protected '+Path(p).name,sha(p)==h0)
qa={'checks_passed':len(checks),'checks':checks,'inputs':plan['inputs'],'plan_sha256':sha(tab/'Fig8_CellFrac_plan.json'),'outputs':outputs,'scratch_PPTX':'pending','visual_review':'pending','native_drag_tested':False,'stats':stats.to_dict(orient='records')}
(base/'results/qa/Fig8_CellFrac_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2)+'\n')
print('Independent count/statistics/PNG/PDF structural checks passed:',len(checks))
