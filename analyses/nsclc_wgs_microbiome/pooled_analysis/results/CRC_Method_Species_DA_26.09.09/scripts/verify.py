#!/usr/bin/env python3
"""Independent pooled matrix, rank/tie/continuity, BH and plot selection checks.
Uses Python/numpy, not R's wilcox.test/p.adjust implementations.
"""
from pathlib import Path
import csv, json, hashlib, math, unicodedata
import numpy as np

out = Path(__file__).resolve().parents[1]
base = out.parents[2]
project = base.parents[1]
used = next(p for p in base.iterdir() if unicodedata.normalize('NFC', p.name) == '사용데이터_모음')
def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))
taxa, blocks, groups, runs = [], [], [], []
for co in ('PRJNA751792','PRJNA1023797','PRJEB22863'):
    with (used/co/(co+'_species_relative_abundance_matrix.csv')).open(newline='') as f:
        rr=csv.reader(f); header=next(rr); rows=list(rr)
    names=header[1:]; ids=[r[0] for r in rows]
    matrix=np.asarray([r[1:] for r in rows],dtype=float)
    meta={r['Run ID']:r for r in read(used/co/(co+'_patient_metadata_WGS_NSCLC_R_NR.csv'))}
    assert len(meta)==len(ids) and set(ids)==set(meta)
    groups.extend(meta[i]['Analysis Group'] for i in ids);runs.extend(ids)
    blocks.append((names,matrix));taxa.extend(x for x in names if x not in taxa)
assert len(set(runs))==824 and groups.count('R')==432 and groups.count('NR')==392
index={t:i for i,t in enumerate(taxa)}
matrix=np.zeros((824,len(taxa)));offset=0
for names,block in blocks:
    matrix[offset:offset+len(block),[index[n] for n in names]]=block
    offset+=len(block)
assert np.isfinite(matrix).all() and (matrix>=0).all()
assert np.max(np.abs(matrix.sum(1)-100))<.001
prev=(matrix>0).mean(0); eligible=np.flatnonzero(prev>=.1)
stats=read(out/'NSCLC_species_DA_results.csv')
assert len(eligible)==len(stats)==541
saved={r['Species_full']:r for r in stats}
assert set(saved)=={taxa[i] for i in eligible}
ridx=np.array(groups)=='R';nidx=~ridx
n1,n2=int(nidx.sum()),int(ridx.sum());N=n1+n2
derived={}
for j in eligible:
    values=np.concatenate([matrix[nidx,j],matrix[ridx,j]])
    unique,inverse,counts=np.unique(values,return_inverse=True,return_counts=True)
    midranks=np.cumsum(counts)-(counts-1)/2
    ranks=midranks[inverse]
    U=ranks[:n1].sum()-n1*(n1+1)/2
    variance=n1*n2/12*(N+1-np.sum(counts**3-counts)/(N*(N-1)))
    z=max(0,abs(U-n1*n2/2)-.5)/math.sqrt(variance)
    p=math.erfc(z/math.sqrt(2))
    meanR=matrix[ridx,j].mean();meanNR=matrix[nidx,j].mean()
    lfc=math.log2((meanNR+1e-6)/(meanR+1e-6))
    derived[taxa[j]]={'U':U,'p':p,'lfc':lfc,'meanR':meanR,'meanNR':meanNR,'prev':prev[j]}
order=sorted(derived,key=lambda k:derived[k]['p']);running=1
for pos in range(len(order)-1,-1,-1):
    k=order[pos];running=min(running,derived[k]['p']*len(order)/(pos+1));derived[k]['q']=running
errors={k:0. for k in ('mean','log2FC','U','p_relative','q_relative','negative_log10_q')}
eligible_hits={'R higher':[],'NR higher':[]}
for species,v in derived.items():
    s=saved[species]
    errors['mean']=max(errors['mean'],abs(v['meanR']-float(s['Mean_R'])),abs(v['meanNR']-float(s['Mean_NR'])))
    errors['log2FC']=max(errors['log2FC'],abs(v['lfc']-float(s['log2FC_NR_vs_R'])))
    errors['U']=max(errors['U'],abs(v['U']-float(s['Wilcoxon_U_NR'])))
    errors['p_relative']=max(errors['p_relative'],abs(v['p']/float(s['p_value'])-1))
    errors['q_relative']=max(errors['q_relative'],abs(v['q']/float(s['BH_q'])-1))
    errors['negative_log10_q']=max(errors['negative_log10_q'],abs(-math.log10(v['q'])-float(s['neg_log10_BH_q'])))
    assert abs(v['prev']-float(s['prevalence']))<1e-14
    cls='Below display thresholds'
    if v['q']<.05 and abs(v['lfc'])>1:
        cls='R higher' if v['lfc']<0 else 'NR higher';eligible_hits[cls].append(species)
    assert s['Display_class']==cls
assert max(errors.values())<1e-10,errors
selected=[]
for direction,keys in eligible_hits.items():
    selected+=sorted(keys,key=lambda k:(derived[k]['q'],k))[:15]
bars=read(out/'NSCLC_species_DA_bar_data.csv')
assert {r['Species_full'] for r in bars}==set(selected)
assert [r['Species_full'] for r in bars]==sorted(selected,key=lambda k:(derived[k]['lfc'],k))
labels=sorted(sum(eligible_hits.values(),[]),key=lambda k:(derived[k]['q'],k))[:15]
assert {r['Species_full'] for r in stats if r['volcano_label']=='TRUE'}==set(labels)
for r in read(out/'qa/input_sha256.csv'):
    assert hashlib.sha256((project/r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
record={'independent_python_verification':'passed','runs':824,'species':541,
        'max_errors':errors,'bar_rows':len(bars),'R_hits':len(eligible_hits['R higher']),
        'NR_hits':len(eligible_hits['NR higher']),'all_input_hashes_unchanged':True,
        'method':'independent midrank/U, tie-corrected normal variance, continuity correction, erfc P, BH step-up'}
(out/'qa/independent_verification.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2))
