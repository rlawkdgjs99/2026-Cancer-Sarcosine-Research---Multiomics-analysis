#!/usr/bin/env python3
"""Finalize and independently check the frozen-data Fig8a display."""
from pathlib import Path
import hashlib,json,shutil,struct,subprocess,tempfile
import numpy as np
import pandas as pd
from PIL import Image
import pypdfium2 as pdfium
from pypdf import PdfReader
base=Path(__file__).resolve().parents[1]
fig=base/'results/figures'; qa_path=base/'results/qa/Fig8a_Cell_Context_QA.json'
qa=json.loads(qa_path.read_text()); tmp=Path('/private/tmp/fig8a_artifact_qa');tmp.mkdir(exist_ok=True)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for v in qa['inputs'].values(): assert sha(v['path'])==v['sha256']
md=pd.read_csv(qa['inputs']['patients']['path']); ref=pd.read_csv(qa['inputs']['reference']['path']).set_index('gene')
for g in ['SARDH','PIPOX']:
 x=np.log2(1+md[g+'_raw']/md.effective_library*1e6)
 np.testing.assert_allclose(x,md[g+'_log2CPM'],rtol=0,atol=1e-12)
 np.testing.assert_allclose((x-ref.loc[g,'mean_log'])/ref.loc[g,'sd_log'],md[g+'_z'],rtol=0,atol=1e-12)
np.testing.assert_allclose((md.SARDH_z+md.PIPOX_z)/2,md.degradation_mean_z,rtol=0,atol=1e-12)
assert list(np.where(md.degradation_mean_z>md.cutoff,'High','Low'))==list(md.group)
assert abs(np.median(md.degradation_mean_z)-qa['original_median'])<1e-14
p=pd.read_csv(base/'results/tables/Fig8a_Patient_Groups.csv')
expected=md[md.treatment=='Post'].sort_values(['degradation_mean_z','Patient']).reset_index(drop=True)
pd.testing.assert_frame_equal(p,expected[list(p.columns)],check_dtype=False,rtol=0,atol=1e-12)
cells=pd.read_csv(qa['inputs']['cells']['path'],usecols=['cell_id','Patient','final_lineage','umap_1','umap_2'])
c=cells[(cells.final_lineage!='Excluded residual doublet')&cells.Patient.isin(p.Patient)]
assert len(c)==78192 and c.Patient.nunique()==12 and not c.cell_id.duplicated().any()
assert np.isfinite(c[['umap_1','umap_2']]).all().all()
assert c.groupby('Patient').size().to_dict()==p.set_index('Patient').all_cells.to_dict()
comp=pd.read_csv(base/'results/tables/Fig8a_Cell_Composition.csv')
assert c.final_lineage.value_counts().to_dict()==comp.set_index('final_lineage').N.to_dict()
master=fig/'Fig8a_Cell_Context.png'; png=fig/'A.png'
with Image.open(master) as im:
 im.load(); before=np.array(im.convert('RGB')); im=im.convert('RGB')
 im.save(master,format='PNG',dpi=(600,600),icc_profile=Path('/System/Library/ColorSync/Profiles/sRGB Profile.icc').read_bytes(),compress_level=6)
shutil.copyfile(master,png)
for f in [master,png,fig/'A.pdf']:
 subprocess.run(['/usr/bin/xattr','-c',str(f)],check=True)
 subprocess.run(['/usr/bin/chflags','nohidden',str(f)],check=True)
 f.chmod(0o644)
 assert not subprocess.check_output(['/usr/bin/xattr',str(f)]).strip()
with Image.open(png) as im:
 assert im.mode=='RGB' and im.size==(3840,3780) and im.info.get('icc_profile')
 assert all(abs(x-600)<.02 for x in im.info['dpi'])
 after=np.array(im);assert np.array_equal(before,after)
 assert (after[0,0]==[255,255,255]).all()
 assert png.read_bytes()[24]==8 and png.read_bytes()[25]==2 and png.read_bytes()[28]==0
assert sha(master)==sha(png)
raw=tmp/'libpng.rgb'
subprocess.run(['/usr/local/bin/Rscript','-e',"a<-commandArgs(TRUE);x<-png::readPNG(a[1]);stopifnot(dim(x)[3]==3);writeBin(as.raw(round(as.vector(aperm(x,c(3,2,1)))*255)),a[2])",str(png),str(raw)],check=True,cwd=str(tmp))
assert raw.read_bytes()==after.tobytes()
r=PdfReader(fig/'A.pdf');assert len(r.pages)==1
text=r.pages[0].extract_text();(tmp/'pdf_text.txt').write_text(text)
for word in ['GSE207422','78,192','SARDH','PIPOX','P02','P10','0.0021','Conventional DC']:
 assert word in text,word
assert 'CAMERA' not in text
pdf=pdfium.PdfDocument(str(fig/'A.pdf')); page=pdf[0];page.render(scale=2.5).to_pil().save(tmp/'pdf_render.png')
qa['independent_python_validation']='Passed: original RNA transform/z/reference, score/median labels, 12-row patient table, 78,192 unique cell count, patient and 15-lineage counts; no source hash changes.'
qa['transport_QA']={'dimensions':[3840,3780],'dpi':600,'mode':'RGB8','interlaced':False,'ICC':'sRGB','master_transport_bytes_identical':True,'Pillow_R_libpng_pixels_identical':True,'conversion_pixel_change':False,'extended_attributes_empty':True,'hidden':False,'native_Finder_drag_tested':False,'scratch_PPTX':'pending'}
qa['pdf_QA']={'pages':1,'required_text_found':True,'render':str(tmp/'pdf_render.png'),'vector_master':True}
qa['outputs']={x.name:{'path':str(x),'sha256':sha(x),'bytes':x.stat().st_size} for x in [master,png,fig/'A.pdf',base/'results/tables/Fig8a_Patient_Groups.csv',base/'results/tables/Fig8a_Cell_Composition.csv']}
qa['scripts']={x.name:sha(x) for x in [base/'scripts/07_plot_Fig8a_Cell_Context_26.09.17.R',Path(__file__).resolve()]}
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2)+'\n')
print('Independent numerical and PNG/PDF structural checks passed. Final PNG:',png)
