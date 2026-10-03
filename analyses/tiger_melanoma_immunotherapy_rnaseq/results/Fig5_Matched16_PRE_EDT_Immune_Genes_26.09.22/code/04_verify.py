from pathlib import Path
import csv,json,math,hashlib,collections
import numpy as np
from openpyxl import load_workbook
p=Path(__file__).resolve().parent.parent
def read(name):
 with (p/name).open() as f:return list(csv.DictReader(f))
def ranks(x):
 x=np.asarray(x);return np.array([1+np.sum(x<v)+(np.sum(x==v)-1)/2 for v in x],dtype=float)
def bh(ps):
 a=np.array(ps);order=np.argsort(a);b=np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];z=np.empty(len(a));z[order]=np.minimum(b,1);return z
checks=0
def eq(a,b,tol=1e-10):
 global checks
 assert math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol),(a,b)
 checks+=1
def mw(x,y):
 x=np.array(x);y=np.array(y);a=np.r_[x,y];n=len(x);m=len(y);u=sum(ranks(a)[:n])-n*(n+1)/2
 ties=np.array(list(collections.Counter(a).values()));v=n*m/12*((n+m+1)-sum(ties**3-ties)/((n+m)*(n+m-1)))
 z=(u-n*m/2-np.sign(u-n*m/2)*.5)/math.sqrt(v)
 return u,math.erfc(abs(z)/math.sqrt(2))
def sr(d):
 a=np.array(d);a=a[a!=0];v=np.array([float(format(abs(x),'.12g')) for x in a]);r=ranks(v);w=sum(r[a>0]);n=len(a);t=np.array(list(collections.Counter(v).values()));var=n*(n+1)*(2*n+1)/24-sum(t**3-t)/48
 z=(w-n*(n+1)/4-np.sign(w-n*(n+1)/4)*.5)/math.sqrt(var)
 return w,math.erfc(abs(z)/math.sqrt(2))
wide=read('tables/patient_data_paired16.csv');long=read('tables/patient_data_long32.csv');met=['CD8_percent','CYT','SARDH','PIPOX','GNMT','DMGDH']
assert len(wide)==len({r['patient_id'] for r in wide})==16 and len(long)==len({r['sample_id'] for r in long})==32
assert collections.Counter(r['response'] for r in wide)=={'R':9,'NR':7}
byid={(r['patient_id'],r['timepoint']):r for r in long}
for r in wide:
 for m in met:
  for t in ['PRE','EDT']:eq(r[t+'_'+m],byid[(r['patient_id'],t)][m])
  eq(r['Delta_'+m],float(r['EDT_'+m])-float(r['PRE_'+m]))
for r in read('tables/paired_statistics.csv'):
 vals=[float(x['Delta_'+r['metric']]) for x in wide if r['group']=='All' or x['response']==r['group']]
 v,pv=sr(vals);eq(v,r['V']);eq(pv,r['P']);eq(np.median(vals),r['delta_median'])
for r in read('tables/delta_R_vs_NR_statistics.csv'):
 x=[float(a['Delta_'+r['metric']]) for a in wide if a['response']=='R'];y=[float(a['Delta_'+r['metric']]) for a in wide if a['response']=='NR']
 u,pv=mw(x,y);eq(u,r['U_R']);eq(pv,r['P']);eq(np.median(x),r['R_delta_median']);eq(np.median(y),r['NR_delta_median'])
for r in read('tables/timepoint_R_vs_NR_statistics.csv'):
 x=[float(a[r['metric']]) for a in long if a['response']=='R' and a['timepoint']==r['timepoint']];y=[float(a[r['metric']]) for a in long if a['response']=='NR' and a['timepoint']==r['timepoint']]
 u,pv=mw(x,y);eq(u,r['U_R']);eq(pv,r['P'])
for name in ['paired_statistics','delta_R_vs_NR_statistics','timepoint_R_vs_NR_statistics']:
 rows=read('tables/'+name+'.csv')
 for family in set(r['BH_family'] for r in rows):
  subset=[r for r in rows if r['BH_family']==family]
  for r,q in zip(subset,bh([float(r['P']) for r in subset])):eq(q,r['BH_q'])
for m in met:
 z=read('tables/PRISM_'+m+'_paired_ID.csv');assert len(z)==16
 for a,b in zip(z,wide):
  assert a['patient_id']==b['patient_id'] and a['sample_PRE']==b['sample_PRE'] and a['sample_EDT']==b['sample_EDT']
  for t in ['PRE','EDT','Delta']:eq(a[t+'_'+m],b[t+'_'+m])
manifest=json.loads((p/'inputs/source_manifest.json').read_text())
for f,h in manifest.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
res={'status':'PASS','numerical_checks':checks,'all36_tests_independent_rank_normal_P_and_BH':True,'paired_ID_and_values':'PASS','source_hashes_unchanged':True}
(p/'qa/independent_verification.json').write_text(json.dumps(res,indent=2));print(res)
