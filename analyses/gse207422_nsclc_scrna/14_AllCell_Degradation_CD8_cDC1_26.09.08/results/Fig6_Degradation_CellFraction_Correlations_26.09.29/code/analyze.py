from pathlib import Path
import numpy as np
import pandas as pd
import itertools, hashlib, json, sys
O=Path(__file__).resolve().parents[1]; B=O.parent
sources=[B/'tables/01_FROZEN_allcell_patient_scores_groups.csv',B/'tables/00_patient_cell_composition.csv',B/'tables/01_composition_descriptive_correlations.csv',B/'Fig6_ImmuneComposition_26.09.21/tables/Patients_Cell_percent.csv']
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
s,c,old,existing=[pd.read_csv(p) for p in sources]
assert len(s)==15 and s.Patient.is_unique and s.Sample.is_unique
assert not c.duplicated(['Patient','final_lineage']).any()
m=c.merge(s[['Patient','Sample','all_cells']],on=['Patient','Sample'],validate='many_to_one')
assert np.allclose(m.cells/m.all_cells,m.fraction,atol=1e-14)
assert np.allclose(m.groupby('Patient').cells.sum().sort_index(),s.set_index('Patient').all_cells.sort_index())
assert s.all_cells.sum()==91844
names=['CD8 T cell','Conventional DC','Epithelial']; keys=['CD8','cDC','Epithelial']
d=s[['Patient','Sample','treatment','histology','group','degradation_mean_z','all_cells']].copy()
for nm,k in zip(names,keys):
 a=c.loc[c.final_lineage==nm].set_index('Patient')
 d[k+'_cells']=d.Patient.map(a.cells); d[k+'_percent']=d.Patient.map(a.fraction)*100
post=d[d.treatment=='Post'].copy(); assert len(post)==12 and post.all_cells.sum()==78192
v=post.merge(existing,on='Patient',suffixes=('_new','_old'),validate='one_to_one')
assert np.allclose(v.degradation_mean_z_new,v.degradation_mean_z_old,atol=1e-14)
assert (v.group_new==v.group_old).all()
for nm,k in zip(names,keys): assert np.allclose(v[k+'_percent'],v[nm+' (%)'],atol=1e-12)
assert d.notna().all().all()
d.to_csv(O/'Patient_data.csv',index=False,encoding='utf-8-sig')

def bh(p):
 p=np.asarray(p); order=np.argsort(p); q=np.empty(len(p)); q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]); return q

def ranks(a): return pd.DataFrame(a).rank(method='average').to_numpy()
def centered(a): return a-a.mean(axis=0)
def coefficients(x,y): return (x[:,None]*y).sum(axis=0)/np.sqrt((x*x).sum()*(y*y).sum(axis=0))
def marginal(x,y,seed,N=1000000):
 rng=np.random.default_rng(seed); obs=coefficients(x,y); den=np.sqrt((x*x).sum()*(y*y).sum(axis=0)); counts=np.zeros(y.shape[1],dtype=np.int64)
 for start in range(0,N,10000):
  n=min(10000,N-start); perm=rng.permuted(np.tile(x,(n,1)),axis=1)
  rr=np.einsum('bi,ij->bj',perm,y)/den
  counts+=(np.abs(rr)>=np.abs(obs)-1e-12).sum(axis=0)
 return obs,counts,(counts+1)/(N+1)
def adjusted(dd):
 cols=['histology'] if len(dd)==12 else ['histology','treatment']
 X=np.column_stack([np.ones(len(dd)),pd.get_dummies(dd[cols],drop_first=True,dtype=float).to_numpy()])
 Z=ranks(dd[['degradation_mean_z']+[k+'_percent' for k in keys]].to_numpy())
 res=Z-X@np.linalg.lstsq(X,Z,rcond=None)[0]
 return res[:,0],res[:,1:],cols
def exact_strat(dd,x,y,cols):
 obs=coefficients(x,y); den=np.sqrt((x*x).sum()*(y*y).sum(axis=0)); totals=np.zeros((1,3))
 strata=[]
 for label,idx in dd.reset_index(drop=True).groupby(cols,sort=True).indices.items():
  ix=np.asarray(idx); perms=np.array(list(itertools.permutations(ix)),dtype=int)
  contrib=np.einsum('bi,ij->bj',x[perms],y[ix])
  totals=(totals[:,None,:]+contrib[None,:,:]).reshape(-1,3); strata.append([str(label),len(ix)])
 rr=totals/den; count=(np.abs(rr)>=np.abs(obs)-1e-12).sum(axis=0)
 return obs,count,count/len(rr),len(rr),strata
rows=[]; loo=[];checks=[]
for scope,dd,seed in [('Post12',post,26092901),('All15',d,26092902)]:
 Z=centered(ranks(dd[['degradation_mean_z']+[k+'_percent' for k in keys]].to_numpy())); x,y=Z[:,0],Z[:,1:]
 rho,count,p=marginal(x,y,seed); q=bh(p)
 _,cnt2,p2=marginal(x,y,seed+100)
 se=np.sqrt(p*(1-p)/1000000); se2=np.sqrt(p2*(1-p2)/1000000)
 assert np.all(np.abs(p-p2)<6*np.sqrt(se**2+se2**2)+2e-6)
 for j,k in enumerate(keys):
  rows.append(dict(scope=scope,analysis='Spearman',cell_type=names[j],key=k,n=len(dd),rho=rho[j],p=p[j],BH_q=q[j],permutations=1000000,extreme_count=int(count[j]),seed=seed,MC_SE=se[j],adjustment='None'))
  checks.append(dict(scope=scope,key=k,p_primary=p[j],p_second_seed=p2[j],second_seed=seed+100))
 if scope=='All15':
  for nm,r in zip(names,rho): assert abs(r-old.set_index('final_lineage').loc[nm,'rho_spearman'])<1e-12
 ax,ay,cols=adjusted(dd); ar,ac,ap,N,strata=exact_strat(dd,ax,ay,cols)
 for j,k in enumerate(keys): rows.append(dict(scope=scope,analysis='Partial_Spearman',cell_type=names[j],key=k,n=len(dd),rho=ar[j],p=ap[j],BH_q=bh(ap)[j],permutations=N,extreme_count=int(ac[j]),seed='',MC_SE=0,adjustment='+'.join(cols)))
 print(scope,'strata',strata,'exact assignments',N,flush=True)
 for patient in dd.Patient:
  sub=dd[dd.Patient!=patient]; z=centered(ranks(sub[['degradation_mean_z']+[k+'_percent' for k in keys]].to_numpy())); r=coefficients(z[:,0],z[:,1:])
  # retain scope-specific adjustment set when deleting one patient
  X=np.column_stack([np.ones(len(sub)),pd.get_dummies(sub[cols],drop_first=True,dtype=float).to_numpy()]); res=z-X@np.linalg.lstsq(X,z,rcond=None)[0]; a=coefficients(res[:,0],res[:,1:])
  for j,k in enumerate(keys): loo.append(dict(scope=scope,omitted_patient=patient,key=k,rho=r[j],partial_rho=a[j]))
stats=pd.DataFrame(rows);stats.to_csv(O/'Correlation_statistics.csv',index=False)
pd.DataFrame(loo).to_csv(O/'Leave_one_out.csv',index=False)
pd.DataFrame(checks).to_csv(O/'Monte_Carlo_check.csv',index=False)
assert hashes=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
(O/'Source_manifest.json').write_text(json.dumps({'sha256':hashes,'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},indent=2))
print(stats.to_string(index=False))
