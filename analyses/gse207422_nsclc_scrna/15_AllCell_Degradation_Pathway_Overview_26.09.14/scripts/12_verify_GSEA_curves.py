"""Independent numeric checks and lossless RGB/sRGB packaging of R figure exports."""
from pathlib import Path
import csv,json,hashlib,sys,math,re
import numpy as np
from PIL import Image,ImageCms
from pypdf import PdfReader
root=Path(sys.argv[1]); tables=root/'results/tables'; qa=root/'results/qa'; figs=root/'results/figures'
def read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
annotations=read(tables/'Fig6_exact_GSEA_curve_annotations.csv')
vertices=read(tables/'Fig6_exact_GSEA_curve_vertices.csv'); hits=read(tables/'Fig6_exact_GSEA_gene_hits.csv'); ranks=read(qa/'GSEA_curve_frozen_ranks.csv')
checks=[]
for a in annotations:
    c,p=a['compartment'],a['pathway']
    rk=[r for r in ranks if r['compartment']==c]
    hh=[h for h in hits if (h['compartment'],h['pathway'])==(c,p)]
    vv=[v for v in vertices if (v['compartment'],v['pathway'])==(c,p)]
    stats=np.array([float(r['statistic']) for r in rk]); N=len(stats)
    assert N==int(a['ranked_genes']) and np.all(np.diff(stats)<=0)
    ix=np.array([int(h['rank'])-1 for h in hh]); mask=np.zeros(N,dtype=bool);mask[ix]=True
    assert [rk[i]['gene'] for i in ix]==[h['gene'] for h in hh]
    factor=2**30/np.abs(stats).sum();factor=math.floor(factor) if factor>=1 else factor
    assert factor==float(a['weight_scale'])
    wt=np.rint(np.abs(stats)*factor); step=np.full(N,-1/(N-len(ix)));step[mask]=wt[mask]/wt[mask].sum()
    walk=np.r_[0,np.cumsum(step)]
    error=max(abs(walk[int(v['rank'])]-float(v['running_ES'])) for v in vv)
    assert error<1e-10 and abs(walk[-1])<1e-10
    es=walk.max() if walk.max()>-walk.min() else walk.min()
    assert abs(es-float(a['ES']))<1e-10
    assert all(abs(float(v['rank_percent'])-100*int(v['rank'])/N)<1e-10 for v in vv)
    peak=int(np.argmax(walk) if es>0 else np.argmin(walk))
    assert peak==int(a['peak_rank'])
    le=[h for h in hh if (int(h['rank'])<=peak if es>0 else int(h['rank'])>peak)]
    assert len(le)==int(a['leading_edge_n'])
    checks.append(dict(compartment=c,pathway=p,max_vertex_error=error,ES_error=abs(es-float(a['ES'])),leading_edge_n=len(le)))
protected=read(qa/'GSEA_curve_input_hashes.csv')
assert all(sha(r['path'])==r['sha256'] for r in protected)
profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
images=[]
for name in ['W','T','D','P']:
    path=figs/(name+'.png')
    with Image.open(path) as im:
        im.load();rgb=im.convert('RGB'); pixelsha=hashlib.sha256(rgb.tobytes()).hexdigest()
    rgb.save(path,format='PNG',dpi=(600,600),icc_profile=profile,optimize=True)
    with Image.open(path) as check:
        check.load();assert check.mode=='RGB' and hashlib.sha256(check.tobytes()).hexdigest()==pixelsha
        assert abs(check.info['dpi'][0]-600)<.1 and check.info['icc_profile']==profile
        size=check.size
        prev=check.copy();prev.thumbnail((1700,1700));prev.save(Path('/private/tmp')/('Fig6_GSEA_'+name+'_preview.png'))
    header=path.read_bytes()[:29];assert header[:8]==b'\x89PNG\r\n\x1a\n' and header[24]==8 and header[25]==2
    pdf=PdfReader(figs/(name+'.pdf'));txt='\n'.join(x.extract_text() for x in pdf.pages)
    expected=12 if name=='P' else 4
    assert len(pdf.pages)==1 and txt.count('BH q =')==expected and len(re.findall(r'NES -?\d+\.\d+',txt))==expected
    assert 'CAMERA' not in txt and 'raw p' not in txt.lower()
    (qa/('GSEA_'+name+'_PDF_text.txt')).write_text(txt)
    images.append(dict(file=str(path),dimensions=size,mode='RGB',bit_depth=8,dpi=600,sRGB=True,pixels_preserved=True,sha256=sha(path),pdf_sha256=sha(figs/(name+'.pdf')),pdf_stat_annotation_pairs=expected))
report=dict(status='PASS',original_inputs_unchanged=len(protected),curve_checks=checks,images=images,scope='Numerical curve reconstruction and file-format checks; native PowerPoint drag not tested.')
(qa/'GSEA_curve_final_verification.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps({'status':'PASS','curves':len(checks),'max_vertex_error':max(r['max_vertex_error'] for r in checks),'images':[(i['file'].split('/')[-1],i['dimensions']) for i in images]},ensure_ascii=False))
