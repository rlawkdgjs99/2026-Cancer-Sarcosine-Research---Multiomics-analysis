from pathlib import Path
import csv,gzip,statistics,hashlib,json,struct,shutil,zlib
root=Path.cwd();base=next(root.glob('2024_Drug_Res_Updates_NSCLC/**/03_lee_fig3_style_FINAL'))
out=base/'results/figures_publication/Fig6_UMAP_Redraw_26.09.29'
def rows(p):
 with (gzip.open(p,'rt',encoding='utf-8-sig') if str(p).endswith('.gz') else p.open(encoding='utf-8-sig')) as f:return list(csv.DictReader(f))
m=rows(out/'tables/Patient_Production_Ranking.csv')
p=[]
for g in ['GNMT','DMGDH']:
 v=[float(r[g+'_log2CPM_ref15']) for r in m];a=statistics.mean(v);sd=statistics.stdev(v);p.append([(x-a)/sd for x in v])
for i,r in enumerate(m):assert abs((p[0][i]+p[1][i])/2-float(r['Production_score_ref15']))<1e-12
rank=sorted(m,key=lambda r:float(r['Production_score_ref15']))
assert {r['Patient'] for r in rank[:3]}=={'P15','P09','P08'}
assert {r['Patient'] for r in rank[-3:]}=={'P04','P07','P12'}
orig=rows(base/'results/figures_publication/UMAP_Contrast_26.09.20/tables/plotted_cells.csv.gz');orig={r['cell_id']:r for r in orig}
d=rows(out/'tables/Fig6d_Plotted_Cells.csv.gz');e=rows(out/'tables/Fig6e_Plotted_Cells.csv.gz')
assert len(e)==len(orig)==90512 and len(d)==41685
assert {r['cell_id'] for r in e}==set(orig)
assert {r['cell_id'] for r in d}=={k for k,v in orig.items() if v['Patient'] in {'P15','P09','P08','P04','P07','P12'}}
for rr in [d,e]:
 for r in rr:
  o=orig[r['cell_id']]
  for k in ['Patient','final_lineage']:assert r[k]==o[k]
  for k in ['umap_1','umap_2','production_module_shifted','degradation_module_shifted']:
   if k in r: assert abs(float(r[k])-float(o[k]))<1e-12
# Verify PNG signature, every chunk CRC, 8-bit truecolor, dimensions and terminal IEND.
pngs=[]
for f in (out/'figures').glob('*.png'):
 buf=f.read_bytes();assert buf[:8]==b'\x89PNG\r\n\x1a\n';pos=8;last=None
 while pos<len(buf):
  n=struct.unpack('>I',buf[pos:pos+4])[0];tag=buf[pos+4:pos+8];data=buf[pos+8:pos+8+n];crc=struct.unpack('>I',buf[pos+8+n:pos+12+n])[0]
  assert zlib.crc32(tag+data)&0xffffffff==crc
  if tag==b'IHDR':w,h,depth,ctype,_,_,_=struct.unpack('>IIBBBBB',data);assert depth==8 and ctype==2
  last=tag;pos+=n+12
 assert last==b'IEND' and pos==len(buf)
 pngs.append({'file':f.name,'width':w,'height':h,'bit_depth':depth,'colour_type':ctype,'sha256':hashlib.sha256(buf).hexdigest()})
short=root/'2024_Drug_Res_Updates_NSCLC/PPT_UMAP_260929';short.mkdir(exist_ok=True)
for nm,new in [('Fig6d_Prod_TopBottom3','D'),('Fig6e_White_Modules','E'),('Production_White','P'),('Degradation_White','G'),('Lineages_All15','L')]:
 for ext in ['png','pdf']:
  src=out/'figures'/f'{nm}.{ext}';dst=short/f'{new}.{ext}';shutil.copy2(src,dst);assert src.read_bytes()==dst.read_bytes()
checks={'independent_score_ranking_verified':True,'cell_IDs_coordinates_lineages_scores_match_original':True,'e_cells':len(e),'d_cells':len(d),'pngs':pngs,'short_path_copies_verified':True}
(out/'qa/Independent_Verification.json').write_text(json.dumps(checks,indent=2))
inputs=[base/'intermediate/01_cell_paper_style_scores.rds',base/'results/tables/plotdata_01_UMAP_display_caps.csv',base/'../27_Fig6_Current_Degradation_CSV_Export_26.09.23/staging/scRNA15_Metadata_Sarcosine.csv']
(out/'qa/Input_Hashes.json').write_text(json.dumps({str(p.resolve().relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2,ensure_ascii=False))
readme='''# Figure 6 UMAP redraw — 2026-09-29\n\nD.png / D.pdf: Production bottom3 vs top3 patients from all15, selected by whole-tumour pseudobulk mean of GNMT and DMGDH sample-SD z-scores of log2(1+TMM-CPM), reference15. High P12/P07/P04, 23,180 annotated cells; Low P08/P09/P15, 18,505 annotated cells. P08 is pretreatment, other five post-treatment. This is a descriptive selected-patient map, not a group-effect test, and does not replace the GSEA grouping. All selected eligible cells retained; no density equalization or downsampling. Same frozen UMAP axes on both maps.\n\nE.png / E.pdf: all90,512 annotated cells across15patients. Original shifted cell module scores unchanged. White at0, warm production and blue degradation. Original separate 99.5th-percentile caps, computed on92,053sourcecells, retained (P≈0.286,D≈1.101). Values above caps saturate only visually. Colour does not directly measure sarcosine concentration or flux. Lower scores may blend into the white background.\n\nP and G: separate production/degradation maps. L: full15 lineage map with refreshed colours. PNGs are opaque8-bit RGB at400dpi; PDFs also available. PNG chunks/CRCs checked; no native PowerPoint insertion test performed. Original PPT not modified.\n\nReproducible R script: 03_lee_fig3_style_FINAL/scripts/10_render_fig6_white_topbottom3_26.09.29.R. Tables, input hashes, software session and independent verification are under the canonical result folder.\n'''
(out/'README.md').write_text(readme)
(short/'README.md').write_text(readme+'\nCanonical results: '+str(out.relative_to(root))+'\n')
shutil.copy2('/private/tmp/verify_umap_redraw.py',out/'qa/verify.py')
print(json.dumps({'checks':'passed','short_path':str(short),'pngs':pngs},indent=2))
