#!/usr/bin/env python3
"""Independent per-cohort raw-matrix, rank/tie/continuity, BH and selection audit.
Does not invoke R's statistical functions. Numerical tolerance: 1e-10.
"""
from pathlib import Path
import csv, json, hashlib, math, unicodedata
import numpy as np

out = Path(__file__).resolve().parents[1]
base = out.parents[3]
project = base.parents[1]
used = next(p for p in base.iterdir() if unicodedata.normalize("NFC", p.name) == "사용데이터_모음")
def read(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f))
records = {}
for co, countR, countNR in (("PRJNA751792",174,164),("PRJNA1023797",225,196),("PRJEB22863",33,32)):
    with (used/co/(co+"_species_relative_abundance_matrix.csv")).open(newline="") as f:
        rr=csv.reader(f);header=next(rr);rows=list(rr)
    taxa=header[1:];runs=[r[0] for r in rows]
    matrix=np.asarray([r[1:] for r in rows],dtype=float)
    md=read(used/co/(co+"_patient_metadata_WGS_NSCLC_R_NR.csv"))
    meta={r["Run ID"]:r for r in md}
    assert len(md)==len(meta)==len(runs)==len(set(runs))==countR+countNR
    assert len(set(r["Sample name"] for r in md))==len(md)
    assert set(runs)==set(meta) and len(taxa)==len(set(taxa))
    groups=np.array([meta[r]["Analysis Group"] for r in runs])
    ridx=groups=="R";nidx=groups=="NR"
    assert int(ridx.sum())==countR and int(nidx.sum())==countNR
    assert np.isfinite(matrix).all() and (matrix>=0).all()
    assert np.max(np.abs(matrix.sum(1)-100))<.001
    prev=(matrix>0).mean(0);eligible=np.flatnonzero(prev>=.1)
    saved={r["Species_full"]:r for r in read(out/co/(co+"_species_DA_results.csv"))}
    assert len(saved)==len(eligible) and set(saved)=={taxa[j] for j in eligible}
    n1,n2=countNR,countR;N=n1+n2;derived={}
    for j in eligible:
        values=np.concatenate([matrix[nidx,j],matrix[ridx,j]])
        unique,inverse,counts=np.unique(values,return_inverse=True,return_counts=True)
        assert len(unique)>1
        midranks=np.cumsum(counts)-(counts-1)/2
        U=midranks[inverse][:n1].sum()-n1*(n1+1)/2
        variance=n1*n2/12*(N+1-np.sum(counts**3-counts)/(N*(N-1)))
        z=max(0,abs(U-n1*n2/2)-.5)/math.sqrt(variance)
        p=math.erfc(z/math.sqrt(2))
        meanR=matrix[ridx,j].mean();meanNR=matrix[nidx,j].mean()
        lfc=math.log2((meanNR+1e-6)/(meanR+1e-6))
        derived[taxa[j]]={"U":U,"p":p,"lfc":lfc,"meanR":meanR,"meanNR":meanNR,
                          "prev":prev[j],"prevR":(matrix[ridx,j]>0).mean(),
                          "prevNR":(matrix[nidx,j]>0).mean()}
    order=sorted(derived,key=lambda k:derived[k]["p"]);running=1
    for pos in range(len(order)-1,-1,-1):
        k=order[pos];running=min(running,derived[k]["p"]*len(order)/(pos+1))
        derived[k]["q"]=running
    errors={k:0. for k in ("mean","log2FC","U","p_relative","q_relative","negative_log10_q")}
    hits={"R higher":[],"NR higher":[]}
    for species,v in derived.items():
        s=saved[species]
        assert s["Cohort"]==co and int(s["n_R"])==countR and int(s["n_NR"])==countNR
        errors["mean"]=max(errors["mean"],abs(v["meanR"]-float(s["Mean_R"])),abs(v["meanNR"]-float(s["Mean_NR"])))
        errors["log2FC"]=max(errors["log2FC"],abs(v["lfc"]-float(s["log2FC_NR_vs_R"])))
        errors["U"]=max(errors["U"],abs(v["U"]-float(s["Wilcoxon_U_NR"])))
        errors["p_relative"]=max(errors["p_relative"],abs(v["p"]/float(s["p_value"])-1))
        errors["q_relative"]=max(errors["q_relative"],abs(v["q"]/float(s["BH_q"])-1))
        errors["negative_log10_q"]=max(errors["negative_log10_q"],abs(-math.log10(v["q"])-float(s["neg_log10_BH_q"])))
        for pykey,col in (("prev","prevalence"),("prevR","prevalence_R"),("prevNR","prevalence_NR")):
            assert abs(v[pykey]-float(s[col]))<1e-14
        cls="Below display thresholds"
        if v["q"]<.05 and abs(v["lfc"])>1:
            cls="R higher" if v["lfc"]<0 else "NR higher";hits[cls].append(species)
        assert s["Display_class"]==cls
    assert max(errors.values())<1e-10,(co,errors)
    selected=[]
    for direction,keys in hits.items():
        selected+=sorted(keys,key=lambda k:(derived[k]["q"],k))[:15]
    bars=read(out/co/(co+"_species_DA_bar_data.csv"))
    assert {r["Species_full"] for r in bars}==set(selected)
    assert [r["Species_full"] for r in bars]==sorted(selected,key=lambda k:(derived[k]["lfc"],k))
    labels=sorted(sum(hits.values(),[]),key=lambda k:(derived[k]["q"],k))[:15]
    assert {k for k,s in saved.items() if s["volcano_label"]=="TRUE"}==set(labels)
    assert {k for k,s in saved.items() if s["bar_selected"]=="TRUE"}==set(selected)
    records[co]={"independent_verification":"passed","n_R":countR,"n_NR":countNR,
                 "tested_species":len(eligible),"q_below_05":sum(v["q"]<.05 for v in derived.values()),
                 "R_hits":len(hits["R higher"]),"NR_hits":len(hits["NR higher"]),
                 "bars":len(bars),"labels":len(labels),"max_errors":errors}
protected=read(out/"qa/input_sha256.csv")
for r in protected:
    assert hashlib.sha256((project/r["path"]).read_bytes()).hexdigest()==r["sha256"],r["path"]
record={"cohorts":records,"protected_files":len(protected),"all_protected_hashes_unchanged":True,
        "method":"Independent midrank/U, tie-corrected normal variance, continuity correction, erfc P and BH step-up",
        "numerical_tolerance":1e-10}
(out/"qa/independent_verification.json").write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps(record,indent=2))

