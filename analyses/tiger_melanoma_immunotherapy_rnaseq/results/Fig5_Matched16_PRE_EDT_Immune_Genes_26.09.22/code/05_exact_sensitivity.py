from pathlib import Path
exec((Path(__file__).resolve().parent/'04_verify.py').read_text().split('wide=read(')[0])
from itertools import combinations
wide=read('tables/patient_data_paired16.csv');long=read('tables/patient_data_long32.csv')
rows=[]
for rr in read('tables/paired_statistics.csv'):
 m=rr['metric'];g=rr['group'];d=np.array([float(r['Delta_'+m]) for r in wide if g=='All' or r['response']==g]);d=d[d!=0]
 ranks2=np.rint(2*ranks(np.array([float(format(abs(x),'.12g')) for x in d]))).astype(int)
 # Exact conditional sign-flip distribution of signed-rank sum.
 sums=np.array([0],dtype=int)
 for r in ranks2:sums=np.r_[sums,sums+r]
 observed=sum(ranks2[d>0]);center=sum(ranks2)/2
 pv=float(np.mean(abs(sums-center)>=abs(observed-center)-1e-10))
 rows.append(dict(analysis='paired',metric=m,group=g,timepoint='',BH_family=rr['BH_family'],P_asymptotic=rr['P'],q_asymptotic=rr['BH_q'],P_exact=pv))
for name in ['delta_R_vs_NR_statistics','timepoint_R_vs_NR_statistics']:
 for rr in read('tables/'+name+'.csv'):
  if name.startswith('delta'):
   data=wide;key='Delta_'+rr['metric'];tp=''
  else:
   tp=rr['timepoint'];data=[r for r in long if r['timepoint']==tp];key=rr['metric']
  vals=np.array([float(r[key]) for r in data]);rank2=np.rint(2*ranks(vals)).astype(int);ix=[i for i,r in enumerate(data) if r['response']=='R'];n=len(ix)
  obs=sum(rank2[ix]);center=n*(len(vals)+1)
  sums=np.array([sum(rank2[list(c)]) for c in combinations(range(len(vals)),n)])
  pv=float(np.mean(abs(sums-center)>=abs(obs-center)-1e-10))
  rows.append(dict(analysis='delta_R_vs_NR' if tp=='' else 'timepoint_R_vs_NR',metric=rr['metric'],group='R_vs_NR',timepoint=tp,BH_family=rr['BH_family'],P_asymptotic=rr['P'],q_asymptotic=rr['BH_q'],P_exact=pv))
for fam in set(r['BH_family'] for r in rows):
 sub=[r for r in rows if r['BH_family']==fam]
 for r,q in zip(sub,bh([r['P_exact'] for r in sub])):r['q_exact']=q
with (p/'tables/exact_small_sample_sensitivity.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for r in rows:
 if (float(r['q_asymptotic'])<.05)!=(r['q_exact']<.05):print('Threshold sensitivity:',r)
print('Exact delta tests',[(r['metric'],r['q_exact']) for r in rows if r['analysis']=='delta_R_vs_NR'])
