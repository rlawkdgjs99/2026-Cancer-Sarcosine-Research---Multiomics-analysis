"""Independent numerical verification from NSCLC exports and original CRC HGMT tables.
No estimates or selections are changed. SciPy rank statistics and NumPy DL implementation
are compared against the R deliverable tables. Deterministic; no new random draws.
"""
from pathlib import Path
import csv, hashlib, json, platform, sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats

OUT = Path(__file__).resolve().parent.parent
ROOT = OUT.parents[5]
BASE = OUT.parents[3]
TABLES = OUT / "tables"
manifest = pd.read_csv(OUT/"qa/input_sha256.csv")
checks = []
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()
def close(name, a, b, tol=1e-9):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    assert a.shape == b.shape, (name,a.shape,b.shape)
    error=float(np.max(np.abs(a-b))) if a.size else 0.
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and error<tol,(name,error)
    checks.append(dict(check=name,n=int(a.size),max_abs_error=error,tolerance=tol))
def bh(p):
    p=np.asarray(p); order=np.argsort(p); q=np.empty_like(p)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q
def rd(name): return pd.read_csv(TABLES/(name+".csv"))
def assoc(m,y,g,co,adjust=False):
    prev=(m>0).mean(); m=m.loc[:,prev>=.10]
    x=stats.rankdata(m.to_numpy(),axis=0); y=stats.rankdata(y)
    if adjust:
        # Intercept + binary phenotype: projection residuals equal within-group
        # centered ranks. This independently checks the R QR implementation.
        for group in np.unique(g):
            ix=g==group
            x[ix]-=x[ix].mean(axis=0); y[ix]-=y[ix].mean()
    else:
        x-=x.mean(axis=0); y-=y.mean()
    rho=np.sum(x*y[:,None],axis=0)/np.sqrt((x*x).sum(axis=0)*(y*y).sum())
    cov=int(adjust); n=len(y); df=n-cov-2
    p=2*stats.t.sf(np.abs(rho*np.sqrt(df/np.maximum(1-rho*rho,np.finfo(float).eps))),df)
    z=pd.DataFrame(dict(cohort=co,species=m.columns,rho=rho,p_value=p,p_adj_BH=bh(p),
        prevalence=prev[m.columns].to_numpy(),n=n,fisher_z_variance=1/(n-cov-3)))
    assert np.all(np.isfinite(rho)) and np.all(np.abs(rho)<1)
    return z
def compare_table(name,computed):
    target=rd(name); fields=["rho","p_value","p_adj_BH","prevalence","n","fisher_z_variance"]
    keys=["cohort","species"]
    a=computed.set_index(keys).sort_index(); b=target.set_index(keys).sort_index()
    assert a.index.equals(b.index),name
    if np.max(np.abs(a[fields].to_numpy()-b[fields].to_numpy()))>=1e-9:
        detail=a[fields].join(b[fields],lsuffix='_Python',rsuffix='_R')
        for f in fields:detail[f+'_difference']=a[f]-b[f]
        detail.to_csv(OUT/'qa'/(name+'_verification_differences.csv'))
        print(name, {f:float((a[f]-b[f]).abs().max()) for f in fields},flush=True)
    close(name,a[fields],b[fields])
def dl(per):
    rows=[]
    for sp,g in per.groupby("species"):
        if len(g)<2: continue
        z=np.arctanh(g.rho.to_numpy()); v=g.fisher_z_variance.to_numpy()
        w=1/v; mu=np.dot(w,z)/w.sum(); Q=np.dot(w,(z-mu)**2)
        tau=max(0,(Q-len(z)+1)/(w.sum()-np.dot(w,w)/w.sum()))
        rw=1/(v+tau); mu=np.dot(rw,z)/rw.sum(); se=1/np.sqrt(rw.sum())
        rows.append(dict(species=sp,k=len(z),pooled_z=mu,pooled_rho=np.tanh(mu),se_z=se,
           ci_lb=np.tanh(mu-stats.norm.ppf(.975)*se),ci_ub=np.tanh(mu+stats.norm.ppf(.975)*se),
           p_value=2*stats.norm.sf(abs(mu/se)),tau2=tau))
    d=pd.DataFrame(rows); d["q_value_BH"]=bh(d.p_value)
    return d
for r in manifest.itertuples():
    assert sha(ROOT/r.path)==r.sha256,("Input changed before verification",r.path)
scores=rd("Fig3e_per_Run_scores")
deg=["K00301","K00302","K00303","K00305","K00306"]; prod=["K00315","K08688"]
families={key:[] for key in ["i_deg","i_prod","j_NSCLC","j_CRC","k_CRC"]}
raw_ko_files=[]; absences=[]
for co,s in scores.groupby("Cohort",sort=False):
    path=BASE/"사용데이터_모음"/co/(co+"_KEGG_KO_relative_abundance_matrix.csv")
    raw_ko_files.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
    ko=pd.read_csv(path,index_col="Run_ID",float_precision='round_trip')
    assert ko.index.is_unique and set(ko.index)==set(s.Run_ID)
    k=ko.reindex(index=s.Run_ID,columns=deg+prod,fill_value=0)
    absences.append(dict(cohort=co,absent_columns=sorted(set(deg+prod)-set(ko.columns))))
    vals=np.c_[k[deg].sum(axis=1),k[prod].sum(axis=1)]
    close(co+"_raw_KO_to_scores",vals,s[["Degradation","Production"]],tol=1e-12)
    close(co+"_raw_ratio",np.log2((vals[:,1]+1e-8)/(vals[:,0]+1e-8)),s.Production_Degradation)
    m=pd.read_csv(BASE/"사용데이터_모음"/co/(co+"_species_relative_abundance_matrix.csv"),index_col="Run_ID",float_precision='round_trip').loc[s.Run_ID]
    m.columns=[col.rsplit("|s__",1)[1] for col in m.columns]
    assert m.columns.is_unique and np.isfinite(m.to_numpy()).all() and (m.to_numpy()>=0).all()
    close(co+"_species_sum",m.sum(axis=1),np.full(len(m),100),tol=.001)
    for label,y,adjust in [("i_deg",vals[:,0],False),("i_prod",vals[:,1],False),("j_NSCLC",vals[:,0],True)]:
        families[label].append(assoc(m,y,s.Group.to_numpy(),co,adjust))
est=rd("Fig3e_statistics"); p=[]; effect=[]
for col in est.ID:
    r=scores.loc[scores.Group=="R",col];nr=scores.loc[scores.Group=="NR",col]
    u=stats.mannwhitneyu(r,nr,alternative="two-sided",method="asymptotic",use_continuity=True)
    p.append(u.pvalue);effect.append(2*u.statistic/(len(r)*len(nr))-1)
close("e_Wilcoxon_P",p,est.wilcox_p);close("e_BH",bh(p),est.q);close("e_effect",effect,est.effect)
print("Verified NSCLC score reconstruction, e tests and within-cohort associations.",flush=True)

crc_cache=ROOT/manifest.loc[manifest.path.str.endswith("sarcosine_KO_per_sample_pooled.csv"),"path"].iloc[0]
cache=pd.read_csv(crc_cache,float_precision='round_trip').set_index("Run.ID")
roster=[]
for co in ["PRJEB10878","PRJEB27928","PRJEB6070","PRJNA429097"]:
    meta_rel=manifest.loc[manifest.path.str.contains(co+"_")&manifest.path.str.contains("/selected_project_"),"path"].iloc[0]
    folder=(ROOT/meta_rel).parent
    # HGMT rows include a trailing extra field; trim to declared header, as the source reader does.
    with (ROOT/meta_rel).open() as f:
        next(f); rows=csv.reader(f,delimiter="\t"); header=next(rows)
        mta=pd.DataFrame([r[:len(header)]+[""]*max(0,len(header)-len(r)) for r in rows if len(r)>1],columns=header)
    mta=mta[(mta["Assay type"]=="WGS") & mta["Phenotype name"].isin(["Health","Colorectal Neoplasms"])].copy()
    assert mta["Run ID"].is_unique
    mta["group"]=np.where(mta["Phenotype name"]=="Health","Healthy","Cancer")
    ena_path=ROOT/manifest.loc[manifest.path.str.endswith(co+"_ena_runs.tsv"),"path"].iloc[0]
    ena=pd.read_csv(ena_path,sep="\t",dtype=str)
    mta=mta.merge(ena,left_on="Run ID",right_on="run_accession",validate="one_to_one")
    assert ((mta["Sample name"]==mta.sample_accession)|(mta["Sample name"]==mta.sample_alias)).all()
    bpath=next(folder.glob("Bacteria_*.txt"))
    b=pd.read_csv(bpath,sep="\t",float_precision='round_trip'); b.columns=["taxa","run","abundance"]
    b=b[b.run.isin(mta["Run ID"]) & b.taxa.str.contains("|s__",regex=False)&~b.taxa.str.contains("|t__",regex=False)].copy()
    assert np.isfinite(b.abundance).all() and (b.abundance>=0).all()
    b["species"]=b.taxa.str.rsplit("|s__",n=1).str[-1]
    mat=b.pivot_table(index="run",columns="species",values="abundance",aggfunc="sum",fill_value=0)
    ids=sorted(set(mat.index)&set(cache.loc[cache.Cohort==co].index))
    runmeta=mta.set_index("Run ID").loc[ids]
    score=cache.loc[ids,deg].sum(axis=1)
    assert (cache.loc[ids,"Group"].to_numpy()==runmeta.group.to_numpy()).all()
    families["k_CRC"].append(assoc(mat.loc[ids],score.to_numpy(),runmeta.group.to_numpy(),co))
    paired=runmeta.loc[runmeta.library_layout=="PAIRED"]
    raw=mat.loc[paired.index].copy()
    raw.index=paired.sample_accession
    # R colMeans accumulates in extended precision. Match that numeric
    # numeric operation to avoid breaking exact ties by machine-epsilon errors.
    agg=pd.DataFrame({bs:gr.to_numpy().astype(np.longdouble).mean(axis=0).astype(float)
        for bs,gr in raw.groupby(level=0,sort=True)}).T
    agg.columns=raw.columns
    agscore=pd.Series(score.loc[paired.index].to_numpy(),index=paired.sample_accession).groupby(level=0).mean()
    groups=paired.groupby("sample_accession").group.first().reindex(agg.index)
    assert (paired.groupby("sample_accession").group.nunique()==1).all()
    families["j_CRC"].append(assoc(agg,agscore.loc[agg.index].to_numpy(),groups.to_numpy(),co,True))
    roster.extend(dict(cohort=co,BioSample=bs,group=groups.loc[bs],n_runs=int((paired.sample_accession==bs).sum()),score=agscore.loc[bs]) for bs in agg.index)
    print(f"Verified raw CRC preprocessing: {co}, {len(ids)} Runs, {len(agg)} paired-library BioSamples.",flush=True)
families={k:pd.concat(v,ignore_index=True) for k,v in families.items()}
table_names={"i_deg":"Fig3i_Degradation_percohort","i_prod":"Fig3i_Production_percohort",
"j_NSCLC":"Fig3j_NSCLC_percohort","j_CRC":"Fig3j_CRC_percohort","k_CRC":"Fig3k_CRC_correlations"}
for key,name in table_names.items():compare_table(name,families[key])
rr=pd.DataFrame(roster).set_index(["cohort","BioSample"]).sort_index()
tr=rd("Fig3j_CRC_BioSample_roster").set_index(["cohort","BioSample"]).sort_index()
assert rr.index.equals(tr.index) and (rr.group==tr.group).all()
close("CRC_BioSample_roster",rr[["n_runs","score"]],tr[["n_runs","score"]])
meta_names={"i_deg":"Fig3i_Degradation_meta","i_prod":"Fig3i_Production_meta","j_NSCLC":"Fig3j_NSCLC_meta","j_CRC":"Fig3j_CRC_meta"}
metas={}
for key,name in meta_names.items():
    a=dl(families[key]).set_index("species").sort_index(); b=rd(name).set_index("species").sort_index()
    assert a.index.equals(b.index)
    close(name+"_all_DL_fits",a,b[a.columns])
    metas[key]=a
j=rd("Fig3j_plot_data"); js=rd("Fig3j_statistics").iloc[0]
for dis,key in [("CRC","j_CRC"),("NSCLC","j_NSCLC")]:
    close("j_plotted_"+dis,metas[key].loc[j.species,"pooled_rho"],j[dis+"_rho"])
close("j_concordance_rho",[stats.spearmanr(j.CRC_rho,j.NSCLC_rho).statistic],[js.spearman_rho])
close("j_sign_concordance",[(np.sign(j.CRC_rho)==np.sign(j.NSCLC_rho)).mean()],[js.sign_concordance])
null=rd("Fig3j_permutation_null").rho
assert len(null)==1000 and np.isfinite(null).all()
close("j_saved_permutation_P",[(1+(null.abs()>=abs(js.spearman_rho)).sum())/1001],[js.conditional_permutation_p])
def core(per,strict,n):
    selected=per[(per.rho>.3)&(per["p_adj_BH" if strict else "p_value"]<.05)]
    count=selected.groupby("species").cohort.nunique()
    return set(count[count==n].index)
cs=core(families["k_CRC"],True,4); ns=core(families["i_deg"],False,3); nsbh=core(families["i_deg"],True,3)
km=rd("Fig3k_plot_data")
assert set(km.loc[km.in_CRC,"species"])==cs and set(km.loc[km.in_NSCLC,"species"])==ns
assert set(km.loc[km.shared,"species"])==cs&ns
assert not nsbh
universe=set(families["k_CRC"].species)&(set(metas["i_deg"].index)|set(metas["i_prod"].index))
ks=rd("Fig3k_statistics").iloc[0]
assert len(universe)==ks.universe==487
close("k_descriptive_hypergeom",[stats.hypergeom.sf(len(cs&ns)-1,len(universe),len(cs),len(ns))],[ks.hypergeom_p])
# Every statistical value plotted in i must come from the corresponding full-family analysis.
idata=rd("Fig3i_plot_data")
for score,key in [("Degradation","i_deg"),("Production","i_prod")]:
    close("i_plotted_"+score,metas[key].loc[idata.species,"pooled_rho"],idata[score+"_rho"])
ip=[stats.mannwhitneyu(idata.loc[idata.enriched=="R",s+"_rho"],idata.loc[idata.enriched=="NR",s+"_rho"],method="asymptotic").pvalue for s in ["Degradation","Production"]]
itest=rd("Fig3i_descriptive_species_comparisons")
close("i_species_comparison_P",ip,itest.p);close("i_species_comparison_BH",bh(ip),itest.q_two_comparisons)
for r in manifest.itertuples():assert sha(ROOT/r.path)==r.sha256
protected=pd.read_csv(OUT/"qa/protected_files.csv")
for r in protected.itertuples():assert sha(ROOT/r.path)==r.sha256
for r in raw_ko_files:assert sha(ROOT/r["path"])==r["sha256"]
report={"status":"PASS","checks":checks,"CRC_BioSample_n":rr.groupby(level=0).size().to_dict(),
"source_files_unchanged":len(manifest),"protected_files_unchanged":len(protected),"raw_KO_files":raw_ko_files,
"absent_KO_columns":absences,"core_shared":sorted(cs&ns),"uniform_BH_NSCLC_core_n":len(nsbh),
"versions":{"Python":sys.version,"NumPy":np.__version__,"pandas":pd.__version__,"SciPy":scipy.__version__,"platform":platform.platform()},
"boundary":"Independently recomputed reported statistics from source inputs. Saved R permutation null was checked; no separate random permutation sequence was used. Overlapping NSCLC projects remain a limitation."}
(OUT/"qa/independent_numerical_validation.json").write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS: {len(checks)} independent numerical checks; all source and protected hashes unchanged.",flush=True)
