#!/usr/bin/env python3
"""Independent count/statistical validation, then publication figure exports."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/mpl_immune_composition')
from pathlib import Path
import json,hashlib,csv,collections,itertools,sys
import numpy as np
import pandas as pd
base=Path(__file__).resolve().parents[1]
plan=json.loads((base/'ANALYSIS_PLAN.json').read_text())
checks=[]
def check(condition,name):
    if not condition: raise AssertionError(name)
    checks.append(name)
for k,item in plan['inputs'].items():
    check(hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256'],f'input unchanged: {k}')
md=pd.read_csv(plan['inputs']['patients']['path']);post=md.loc[md.treatment=='Post'].sort_values('Patient').reset_index(drop=True)
check(len(post)==12 and post.Patient.nunique()==12 and (post.group=='High').sum()==7,'original 12 patients / High7 Low5')
comp=pd.read_csv(plan['inputs']['composition']['path'])
# Recount the original cell annotation, including zero-aware dictionary lookups.
counts=collections.Counter();excluded=collections.Counter()
with open(plan['inputs']['cells']['path']) as f:
    for r in csv.DictReader(f):
        if r['Patient'] not in set(post.Patient): continue
        if r['final_lineage']=='Excluded residual doublet':
            excluded[r['Patient']]+=1; continue
        counts[(r['Patient'],r['final_lineage'])]+=1
sub=comp[comp.Patient.isin(post.Patient)]
for r in sub.itertuples():
    check(counts[(r.Patient,r.final_lineage)]==r.cells,f'source recount {r.Patient}/{r.final_lineage}')
for r in post.itertuples():
    check(sum(v for (pid,l),v in counts.items() if pid==r.Patient)==r.all_cells,f'all-cell denominator {r.Patient}')
check(sum(counts.values())==78192,'78192 cells')
lineages=plan['targets'];check(len(lineages)==13,'13 immune categories')
check(set(sub.final_lineage)-{'Epithelial','CAF'}==set(lineages),'all immune categories retained')
check(sum(v for (pid,l),v in counts.items() if l in lineages)==68093,'68093 immune cells')
values=pd.read_csv(base/'tables/patient_values.csv');stats=pd.read_csv(base/'tables/statistics.csv').set_index('cell_type').loc[lineages]
observed=(post.group=='High').astype(int).to_numpy()
hist=(post.histology=='Squamous').astype(int).to_numpy()
ixs=[np.flatnonzero(hist==h) for h in [0,1]]
choices=[list(itertools.combinations(ix,int(observed[ix].sum()))) for ix in ixs]
assign=[]
for a,z in itertools.product(*choices):
    group=np.zeros(12);group[list(a)+list(z)]=1;assign.append(group)
check(len(assign)==300,'300 exhaustive stratified assignments')
def ols(y,g):
    X=np.column_stack([np.ones(12),hist,g]);beta=np.linalg.lstsq(X,y,rcond=None)[0]
    residual=y-X@beta;sigma=(residual@residual)/(12-X.shape[1]);cov=sigma*np.linalg.inv(X.T@X)
    return beta[2],beta[2]/np.sqrt(cov[2,2])
P=[]
for l in lineages:
    num=np.array([counts[(pid,l)] for pid in post.Patient]);den=post.all_cells.to_numpy()
    check(np.all((num>0)&(num<den)),f'no pseudocount required {l}')
    pct=100*num/den;z=(pct-pct.mean())/pct.std(ddof=1);y=np.log(num/(den-num))
    r=values[values.cell_type==l].set_index('Patient').loc[post.Patient]
    check(np.allclose(r.percent,pct,rtol=0,atol=1e-11) and np.allclose(r.row_z,z,rtol=0,atol=1e-11),f'percent and row z {l}')
    check(list(r.group)==list(post.group) and list(r.histology)==list(post.histology),f'frozen metadata {l}')
    beta,t=ols(y,observed);ts=np.array([ols(y,g)[1] for g in assign]);ne=int((np.abs(ts)>=abs(t)-1e-12).sum());p=ne/300;P.append(p)
    st=stats.loc[l]
    check(np.isclose(st.log_ratio_beta,beta,rtol=0,atol=1e-11) and np.isclose(st.t_statistic,t,rtol=0,atol=1e-11) and ne==st.extreme_assignments and np.isclose(p,st.raw_P),f'independent OLS and permutation P {l}')
    check(np.isclose(pct[observed==0].mean(),st.mean_Low_percent) and np.isclose(pct[observed==1].mean(),st.mean_High_percent) and np.isclose(pct[observed==1].mean()-pct[observed==0].mean(),st.difference_pp),f'means and difference {l}')
P=np.array(P);order=np.argsort(P);q_sorted=np.minimum.accumulate((P[order]*len(P)/np.arange(1,len(P)+1))[::-1])[::-1];q=np.empty(len(P));q[order]=np.minimum(q_sorted,1)
check(np.allclose(q,stats.BH13_q,atol=1e-12),'independent BH13')
old=pd.read_csv(plan['inputs']['old_stats']['path']);old=old[old.scope=='all_singlets'].set_index('cell_type')
for l in ['CD8 T cell','Conventional DC']:
    check(np.isclose(stats.loc[l,'raw_P'],old.loc[l,'raw_P']) and np.isclose(stats.loc[l,'difference_pp'],old.loc[l,'raw_mean_difference_pp']),f'prior comparison reproduced {l}')
# Deterministic patient order preserves each original group; no clustering.
patients=post.assign(group_order=(post.group=='High').astype(int)).sort_values(['group_order','Patient']).Patient.tolist()
metadata=post.set_index('Patient').loc[patients]
mat=values.pivot(index='cell_type',columns='Patient',values='row_z').loc[lineages,patients]
raw=values.pivot(index='cell_type',columns='Patient',values='percent').loc[lineages,patients]
mat.to_csv(base/'tables/heatmap_row_z_ordered.csv');raw.to_csv(base/'tables/PRISM_percent_ordered.csv');metadata.to_csv(base/'tables/heatmap_column_annotations.csv')
(base/'qa/verification.json').write_text(json.dumps({'checks_passed':len(checks),'checks':checks,'independent_tests':'Source annotation recount, numpy OLS/permutation/BH13, original CD8/cDC P reproduction','versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}},indent=2)+'\n')
print(f'PASS {len(checks)} checks')
print(stats[['mean_Low_percent','mean_High_percent','difference_pp','raw_P','BH13_q']].to_string())
