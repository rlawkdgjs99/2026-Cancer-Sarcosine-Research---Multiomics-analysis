"""Independent Python numerical and raster QA for the separate CRC-method branch.

No source data, existing figure or manuscript is modified. Only this branch's
three new PNGs are standardized, without changing any decoded RGB value.
"""
from pathlib import Path
import hashlib
import json
import struct
import subprocess
import tempfile
import unicodedata
import math
import numpy as np
import pandas as pd
from PIL import Image, ImageCms

POOL = Path(__file__).resolve().parent.parent
OUT = POOL / 'results/CRC_Method_7KO_NSCLC_26.09.07'
TABLE = OUT / 'tables'
QA = OUT / 'qa'
TARGET = ['K00301','K00302','K00303','K00305','K00306','K00315','K08688']
COHORTS = ['PRJNA751792','PRJNA1023797','PRJEB22863']
checks = []

def check(condition, label):
    assert bool(condition), label
    checks.append(label)

def same(a, b, label, atol=1e-12):
    check(np.allclose(a,b,rtol=1e-12,atol=atol,equal_nan=False),label)

def read(p, **kwargs):
    return pd.read_csv(p, encoding='utf-8-sig',float_precision='round_trip',**kwargs)

def sha(p):
    with open(p,'rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def child(parent, name):
    return next(p for p in parent.iterdir() if unicodedata.normalize('NFC',p.name)==name)

def bh(p):
    p=np.array(p,dtype=float); order=np.argsort(p); result=np.empty(len(p))
    result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return result

def effects(values, group):
    ranks=pd.DataFrame(values).rank(method='average',axis=0).to_numpy()
    n=group.sum(); m=(~group).sum()
    return 2*(ranks[group].sum(axis=0)-n*(n+1)/2)/(n*m)-1

def independent_wilcoxon_p(a,b):
    # Independent asymptotic Mann-Whitney formula, with tie variance and
    # continuity correction matching the specified R Wilcoxon test.
    combined=np.concatenate([a,b]); n=len(a); m=len(b); total=n+m
    ranks=pd.Series(combined).rank(method='average').to_numpy()
    u=ranks[:n].sum()-n*(n+1)/2
    counts=np.unique(combined,return_counts=True)[1].astype(float)
    variance=n*m/12*((total+1)-np.sum(counts**3-counts)/(total*(total-1)))
    z=(abs(u-n*m/2)-.5)/math.sqrt(variance)
    return min(1.,math.erfc(z/math.sqrt(2)))

manifest=read(QA/'input_sha256.csv')
check(all(sha(p)==s for p,s in zip(manifest.path,manifest.SHA256)), 'All 14 recorded input hashes match before QA')
used=child(POOL.parent,'사용데이터_모음')
frames=[]; metadata=[]
for cohort in COHORTS:
    x=read(used/cohort/f'{cohort}_KEGG_KO_relative_abundance_matrix.csv').set_index('Run_ID')
    md=read(used/cohort/f'{cohort}_patient_metadata_WGS_NSCLC_R_NR.csv').set_index('Run ID')
    check(set(x.index)==set(md.index),f'{cohort} original Run ID equality')
    check(x.index.is_unique and md.index.is_unique,f'{cohort} unique Run IDs')
    md=md.loc[x.index]
    check((md.Cohort==cohort).all(), f'{cohort} cohort label')
    frames.append(x);metadata.append(md)
ns=pd.concat(frames).fillna(0)
md=pd.concat(metadata)
del frames
group=(md['Analysis Group']=='R').to_numpy()
check(ns.shape[0]==824 and ns.index.is_unique,'824 unique Run records')
check(group.sum()==432 and (~group).sum()==392,'432 R and 392 NR')
check((md.groupby('Sample name').Cohort.nunique()>1).sum()==283,'283 cross-cohort repeated Sample names documented, no independence claim')
export=read(TABLE/'NSCLC_seven_KO_per_Run.csv').set_index('Run_ID')
check(list(export.index)==list(ns.index),'Per-Run export row order equals source')
check(list(export.Group)==list(md['Analysis Group']),'Export group labels equal original metadata')
same(export[TARGET],ns[TARGET],'All 5768 exported target cells reproduce source values',atol=1e-15)

scores=pd.DataFrame(index=ns.index)
scores['Degradation']=ns[TARGET[:5]].sum(axis=1)
scores['Production']=ns[TARGET[5:]].sum(axis=1)
scores['Production_Degradation']=np.log2((scores.Production+1e-8)/(scores.Degradation+1e-8))
saved=read(TABLE/'NSCLC_CRC_defined_pathway_scores.csv').set_index('Run_ID')
same(scores,saved[scores.columns],'All 2472 pathway cells reproduce CRC five+two+ratio definition')
check(np.isfinite(scores).all().all(),'All 824 ratios finite, no unreported exclusions')

for label,values,statsfile in [
    ('7-KO',ns[TARGET],'NSCLC_seven_KO_pooled_statistics.csv'),
    ('3-score',scores,'NSCLC_pathway_pooled_statistics.csv')]:
    stats=read(TABLE/statsfile).set_index('ID')
    check(list(stats.index)==list(values.columns),f'{label} exact feature family')
    pvalues=[]
    for feature in values:
        vector=values[feature].to_numpy();a=vector[group];b=vector[~group]
        pairwise=np.sign(a[:,None]-b[None,:]).mean()
        independent_p=independent_wilcoxon_p(a,b)
        same(pairwise,stats.loc[feature,'effect'],f'{feature} direct pairwise effect')
        same(independent_p,stats.loc[feature,'wilcox_p'],f'{feature} independent tie-corrected two-sided P')
        pvalues.append(independent_p)
        boot=read(TABLE/f'bootstrap_{feature}.csv')
        check(len(boot)==5000 and boot.effect.between(-1,1).all(),f'{feature} 5000 valid bootstrap effects')
        same(np.quantile(boot.effect,[.025,.975],method='linear'),
             stats.loc[feature,['ci_low','ci_high']].astype(float),f'{feature} independently reconstructed percentile CI')
    same(bh(pvalues),stats.q,f'{label} independent BH family correction')

crc_path=next(Path(p) for p in manifest.path if str(p).endswith('CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv'))
crc=read(crc_path).set_index('Run_ID')
cg=(crc.Analysis_Group=='Healthy').to_numpy()
check(len(crc)==1647 and cg.sum()==745 and (~cg).sum()==902,'CRC 1647 Run records, 745 Healthy and 902 Cancer')
crc=crc[[c for c in crc if c.startswith('K') and len(c)==6]]
npv=(ns>0).mean();cpv=(crc>0).mean()
common=sorted(set(npv[npv>=.1].index)&set(cpv[cpv>=.1].index))
display=sorted(set(common)|set(TARGET))
cross=read(TABLE/'CRC_NSCLC_pooled_KO_concordance.csv').set_index('KO')
check(list(cross.index)==display,'Concordance target/background membership independently reconstructed')
check(cross.targeted.sum()==7,'All seven targets have actual finite positions')
same(effects(ns[display].to_numpy(),group),cross.NSCLC_effect,'All 6982 NSCLC plotted effects independently reconstructed')
same(effects(crc[display].to_numpy(),cg),cross.CRC_effect,'All 6982 CRC plotted effects independently reconstructed')
same(npv.loc[display],cross.NSCLC_prevalence,'NSCLC per-feature prevalence')
same(cpv.loc[display],cross.CRC_prevalence,'CRC per-feature prevalence')
bg=cross.loc[common]
summary=read(TABLE/'concordance_summary.csv').iloc[0]
check(len(common)==summary.n_background and len(display)==summary.n_display,'6978 background and 6982 total points')
same(np.corrcoef(bg.CRC_effect.rank(),bg.NSCLC_effect.rank())[0,1],summary.spearman_rho,'Independent descriptive Spearman rho')
agreement=np.sign(bg.CRC_effect)==np.sign(bg.NSCLC_effect)
check(agreement.sum()==summary.same_direction,'4476 directionally concordant background KOs')
same(100*agreement.mean(),summary.direction_agreement_percent,'Independent directional agreement percentage')
check(all(sha(p)==s for p,s in zip(manifest.path,manifest.SHA256)),'All recorded inputs unchanged after numerical QA')
(QA/'independent_numerical_validation.json').write_text(json.dumps({
    'status':'PASS','checks':checks,'numpy_version':np.__version__,'pandas_version':pd.__version__,
    'limit':'Arithmetic/implementation validation does not establish independence of overlapping Sample names or replace cohort-adjusted inference. CI quantiles were checked independently from the saved R bootstrap realizations; R RNG was not recreated in Python.'},indent=2)+'\n')

# Profile embedding does not perform a color transform. Original decoded RGB
# values are compared before and after. These are this branch's new files only.
profile=Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc').read_bytes()
profile_name=ImageCms.getProfileDescription(ImageCms.ImageCmsProfile(__import__('io').BytesIO(profile))).strip()
pngs=[('NSCLC_CRC7_Individual_KOs.png','PPT_INSERT_NSCLC_7KO.png'),
      ('NSCLC_CRC7_Pathway_Balance.png','PPT_INSERT_NSCLC_Scores.png'),
      ('CRC_NSCLC_CRC7_Concordance.png','PPT_INSERT_CRC_NSCLC.png')]
raster=[]
for name,transport in pngs:
    p=OUT/'figures'/name
    with Image.open(p) as im:
        im.load()
        if 'A' in im.getbands():
            assert np.asarray(im.getchannel('A')).min()==255
        rgb=im.convert('RGB'); before=np.asarray(rgb).copy()
    with tempfile.NamedTemporaryFile(suffix='.png',delete=False) as tmp:
        tmp_path=Path(tmp.name)
    rgb.save(tmp_path,format='PNG',icc_profile=profile,dpi=(600,600))
    new=tmp_path.read_bytes(); tmp_path.unlink()
    p.write_bytes(new)
    dest=OUT/'figures'/transport;dest.write_bytes(new)
    for fp in (p,dest):
        fp.chmod(0o644)
        subprocess.run(['/usr/bin/xattr','-c',str(fp)],check=True)
        subprocess.run(['/usr/bin/chflags','nohidden',str(fp)],check=True)
        with Image.open(fp) as reread:
            assert reread.mode=='RGB' and np.array_equal(before,np.asarray(reread))
            assert reread.info['icc_profile']==profile
            assert np.allclose(reread.info['dpi'],[600,600],atol=.001)
        width,height,depth,colortype,compression,filter_method,interlace=struct.unpack('>IIBBBBB',new[16:29])
        assert depth==8 and colortype==2 and interlace==0
    assert sha(p)==sha(dest)
    raster.append({'file':str(p),'transport':str(dest),'sha256':sha(p),'width':width,'height':height,
                   'mode':'RGB','depth':8,'interlace':0,'dpi':600,'icc':profile_name,
                   'raw_rgb_sha256':hashlib.sha256(before.tobytes()).hexdigest(),'pixel_change':False})
(QA/'raster_validation.json').write_text(json.dumps({'status':'PASS','files':raster},indent=2)+'\n')
print(f'PASS: {len(checks)} numerical checks; 3 RGB/ICC PNGs and byte-identical transport copies')
