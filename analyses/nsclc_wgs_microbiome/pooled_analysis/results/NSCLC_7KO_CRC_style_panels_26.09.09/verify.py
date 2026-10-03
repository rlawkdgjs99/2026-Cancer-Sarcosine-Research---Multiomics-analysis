#!/usr/bin/env python3
"""Independent NumPy/stdlib checks against the R analysis; read-only outputs.
Run using a Python environment with NumPy (verified with Python3.11/NumPy2.1.3).
Checks every estimate/p/q; bootstrap endpoints checked separately in R.
"""
from pathlib import Path
import csv, math, hashlib, json
import numpy as np

OUT=Path(__file__).resolve().parent
CSV=OUT/'csv'
def rows(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
sources=rows(CSV/'Input_source_SHA256.csv')
for x in sources:
    assert hashlib.sha256(Path(x['Source']).read_bytes()).hexdigest()==x['SHA256']
def matrix(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        r=csv.reader(f);h=next(r);d=list(r)
    return h,np.array([x[:4] for x in d]),np.array([x[4:] for x in d],dtype=float)
meta=rows(sources[0]['Source'])
sh,sm,S=matrix(sources[1]['Source']);kh,km,K=matrix(sources[2]['Source']);sch,scm,SC=matrix(sources[3]['Source'])
assert len(meta)==824 and len(set(x['Run_ID'] for x in meta))==824
assert np.array_equal(sm,km) and np.array_equal(sm,scm)
assert (S>=0).all() and np.isfinite(S).all() and np.max(np.abs(S.sum(1)-100))<.001
assert np.max(np.abs(K[:,:5].sum(1)-SC[:,0]))<1e-12
assert np.max(np.abs(K[:,5:].sum(1)-SC[:,1]))<1e-12
assert np.max(np.abs(np.log2((SC[:,1]+1e-8)/(SC[:,0]+1e-8))-SC[:,2]))<1e-12
co=['PRJNA751792','PRJNA1023797','PRJEB22863'];kos=kh[4:];taxa=sh[4:]
ixs={'Pooled':np.arange(824),**{c:np.where(sm[:,0]==c)[0] for c in co}}
assert [len(ixs[c]) for c in co]==[338,421,65]
assert [np.sum(sm[ixs[c],2]=='R') for c in co]==[174,225,33]
ki={x:i for i,x in enumerate(kos)};si={x:i for i,x in enumerate(taxa)}
def ranks(x):
    _,inv,n=np.unique(x,return_inverse=True,return_counts=True)
    return (np.cumsum(n)-(n-1)/2)[inv]
def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);good=np.where(np.isfinite(p))[0]
    order=good[np.argsort(p[good])];m=len(order)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*m/np.arange(1,m+1))[::-1])[::-1]);return q
def wilcox(x,y):
    n,m=len(x),len(y);z=np.r_[x,y];rk=ranks(z);u=rk[:n].sum()-n*(n+1)/2
    _,cnt=np.unique(z,return_counts=True)
    variance=n*m/12*((n+m+1)-np.sum(cnt.astype(float)**3-cnt)/((n+m)*(n+m-1)))
    effect=2*u/(n*m)-1
    if variance<=0:return math.nan,math.nan
    dev=u-n*m/2;zz=(dev-.5*np.sign(dev))/math.sqrt(variance)
    return effect,math.erfc(abs(zz)/math.sqrt(2))
# Regularized incomplete beta continued fraction for two-sided t-tail p.
def cf(a,b,x):
    tiny=1e-300;qab=a+b;qap=a+1;qam=a-1;c=1.;d=1-qab*x/qap
    if abs(d)<tiny:d=tiny
    d=1/d;h=d
    for m in range(1,1001):
        aa=m*(b-m)*x/((qam+2*m)*(a+2*m))
        d=1+aa*d;c=1+aa/c
        if abs(d)<tiny:d=tiny
        if abs(c)<tiny:c=tiny
        d=1/d;h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+2*m)*(qap+2*m))
        d=1+aa*d;c=1+aa/c
        if abs(d)<tiny:d=tiny
        if abs(c)<tiny:c=tiny
        d=1/d;delta=d*c;h*=delta
        if abs(delta-1)<3e-14:return h
    raise ArithmeticError('Incomplete beta did not converge')
def ibeta(x,a,b):
    if x<=0:return 0.
    if x>=1:return 1.
    bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
    if x<(a+1)/(a+b+2):return bt*cf(a,b,x)/a
    return 1-bt*cf(b,a,1-x)/b
def spearman(x,y):
    x=ranks(x);y=ranks(y);x-=x.mean();y-=y.mean();den=np.linalg.norm(x)*np.linalg.norm(y)
    if den==0:return math.nan,math.nan
    r=float(x@y/den);return r,ibeta(max(0,1-r*r),(len(x)-2)/2,.5)
maxerr={};counts={}
def check(v,s,kind,tol=2e-10):
    s=float(s) if str(s) else math.nan
    if math.isnan(v):assert math.isnan(s),(kind,v,s);return
    assert math.isfinite(s),(kind,v,s)
    err=abs(v-s);maxerr[kind]=max(maxerr.get(kind,0),err)
    assert err<tol,(kind,v,s,err)
    if 0<abs(v)<1e-7:assert abs((s-v)/v)<2e-7,(kind,'relative',v,s)

forests=rows(CSV/'KO_7genes_3cohorts_forest.csv')
for k in kos:
    group=[r for r in forests if r['KO']==k];ps=[]
    for r in group:
        ix=ixs[r['Cohort']];a=K[ix,ki[k]];nr=sm[ix,2]=='NR';est,p=wilcox(a[nr],a[~nr]);ps.append(p)
        check(est,r['effect_NR_minus_R'],'forest_effect');check(p,r['p_value'],'forest_p')
        if np.isfinite(est):
            # Second independent point-estimate method: explicit pair comparisons.
            pair=(a[nr,None]>a[~nr]).sum()-(a[nr,None]<a[~nr]).sum()
            check(pair/(nr.sum()*(~nr).sum()),r['effect_NR_minus_R'],'forest_paircount')
            assert -1<=float(r['CI95_low'])<=float(r['CI95_high'])<=1
            assert float(r['CI95_low'])<=est<=float(r['CI95_high'])
            assert r['fill_status']=='Filled'
        else:
            assert r['estimable']=='FALSE' and r['CI95_low']==r['CI95_high']==r['BH_q']==''
    for r,q in zip(group,bh(ps)):check(q,r['BH_q'],'forest_q')
counts['forest_rows']=len(forests);counts['all_zero']=sum(r['status']=='All zero' for r in forests)
assert counts['all_zero']==4
pooled=rows(CSV/'KO_7genes_pooled_lollipop.csv');ps=[]
for r in pooled:
    j=ki[r['KO']];nr=sm[:,2]=='NR';a=K[nr,j];b=K[~nr,j];e,p=wilcox(a,b);ps.append(p)
    check(e,r['effect_NR_minus_R'],'pooled_KO_effect');check(p,r['p_value'],'pooled_KO_p')
    check(np.log2((a.mean()+1e-8)/(b.mean()+1e-8)),r['log2FC_NR_vs_R'],'pooled_KO_log2FC')
    check(np.mean(K[:,j]>0),r['prevalence'],'pooled_prevalence')
for r,q in zip(pooled,bh(ps)):check(q,r['BH_q'],'pooled_KO_q')

da=rows(CSV/'Species_DA_all_eligible.csv')
assoc=rows(CSV/'Sarcosine_associated_species_DA_pooled_and_3cohorts.csv')
for source,qname in [(da,'BH_q_global'),(assoc,'BH_q_selected')]:
    for scope,ix in ixs.items():
        group=[r for r in source if r['Analysis']==scope];ps=[]
        if source is da:
            eligible=np.where(np.mean(S[ix]>0,axis=0)>=.1)[0]
            assert set(taxa[j] for j in eligible)==set(r['Species_full'] for r in group)
        for r in group:
            x=S[ix,si[r['Species_full']]];nr=sm[ix,2]=='NR';a=x[nr];b=x[~nr];_,p=wilcox(a,b);ps.append(p)
            check(p,r['p_value'],'species_DA_p')
            check(np.log2((a.mean()+1e-6)/(b.mean()+1e-6)),r['log2FC_NR_vs_R'],'species_DA_log2FC')
        for r,q in zip(group,bh(ps)):check(q,r[qname],'species_DA_q')
counts['species_DA_tests']=len(da);counts['selected_DA_rows']=len(assoc)
cor=rows(CSV/'Species_7KO_Spearman_all_tests.csv');score=rows(CSV/'Species_pathway_score_Spearman_all_tests.csv')
for source in [cor,score]:
    families={}
    for r in source:
        scope=r['Analysis'];ix=ixs[scope];j=si[r['Species_full']]
        x=K[ix,ki[r['KO']]] if source is cor else SC[ix,0 if r['Score']=='Degradation' else 1]
        rho,p=spearman(x,S[ix,j]);check(rho,r['rho'],'spearman_rho');check(p,r['p_value'],'spearman_p')
        key=(scope,'7KO') if source is cor else (scope,r['Score'])
        families.setdefault(key,[]).append((r,p))
    for g in families.values():
        qs=bh([p for _,p in g])
        for (r,p),q in zip(g,qs):check(q,r['BH_q'],'spearman_q')
counts['KO_species_rows']=len(cor);counts['score_species_rows']=len(score)
# Selection reconstructed from complete testing universe, no selected-only BH.
passing={}
for r in cor:
    if r['Analysis']=='Pooled' and r['BH_q'] and float(r['BH_q'])<.05 and abs(float(r['rho']))>.3:
        passing.setdefault(r['Species_full'],[]).append(r)
ranked=sorted(passing,key=lambda s:(-max(abs(float(r['rho'])) for r in passing[s]),s))
top=ranked[:20];globalq={r['Species_full']:float(r['BH_q_global']) for r in da if r['Analysis']=='Pooled'}
selected=[s for s in top if globalq[s]<.05][:6]
assert selected==[r['Species_full'] for r in rows(CSV/'Scatter_selected_species.csv')]
assert len(selected)==5 and len(top)==12
for r in assoc:
    assert r['Species_full'] in top
    assert r['Plot_bar']==str(bool(r['BH_q_selected'] and float(r['BH_q_selected'])<.05)).upper()
member=rows(CSV/'Venn_3cohorts_membership.csv');regions=rows(CSV/'Venn_3cohorts_region_counts.csv')
assert len(regions)==28
for cat in ['Production-associated','Degradation-associated','NR-enriched','R-enriched']:
    expected={}
    for c in co:
        if cat.endswith('-associated'):
            expected[c]={r['Species_full'] for r in score if r['Analysis']==c and r['Score']==cat.split('-')[0] and float(r['rho'])>.3 and float(r['BH_q'])<.05}
        else:
            expected[c]={r['Species_full'] for r in da if r['Analysis']==c and float(r['BH_q_global'])<.05 and (float(r['log2FC_NR_vs_R'])>1 if cat=='NR-enriched' else float(r['log2FC_NR_vs_R'])< -1)}
        actual={r['Species_full'] for r in member if r['Category']==cat and r[c]=='TRUE'}
        assert actual==expected[c]
    for r in [r for r in regions if r['Category']==cat]:
        n=sum(all((s in expected[c])==(r[c]=='TRUE') for c in co) for s in set.union(*expected.values()))
        assert n==int(r['N'])
raw=rows(CSV/'Scatter_7KO_and_selected_species_per_Run.csv')
assert [r['Run_ID'] for r in raw]==list(sm[:,1])
for i,r in enumerate(raw):
    for j,k in enumerate(kos):assert float(r[k])==K[i,j]
    for s in selected:assert float(r[s])==S[i,si[s]]
print(json.dumps({'status':'PASS','counts':counts,'selected_species':[s.split('|s__')[-1] for s in selected],'max_absolute_errors':maxerr},indent=2))
