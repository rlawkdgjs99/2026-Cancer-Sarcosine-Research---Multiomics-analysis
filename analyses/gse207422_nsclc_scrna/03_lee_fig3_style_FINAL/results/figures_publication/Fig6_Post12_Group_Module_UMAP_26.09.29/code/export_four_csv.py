from pathlib import Path
import pandas as pd, numpy as np, hashlib, json
root=Path(__file__).resolve().parents[5]
# root is analysis_sarcosine_FINAL_26.08.25
pub=root/'03_lee_fig3_style_FINAL/results/figures_publication'
nsclc=root.parent.parent
folder=next(x for x in nsclc.iterdir() if x.is_dir() and 'csv' in x.name)
out=next(x for x in folder.iterdir() if x.is_dir() and x.name.startswith('csv'))
legacy=pub/'Fig6_Post12_Group_UMAP_26.09.29'
mod=pub/'Fig6_Post12_Group_Module_UMAP_26.09.29'
gsea=root/'19_Degradation_Low_HR_Overview_26.09.20'
sources=[legacy/'tables/Plotted_Cells.csv.gz',legacy/'tables/Patient_Groups.csv',mod/'tables/Plotted_Cells.csv.gz',mod/'tables/Plot_Audit.csv',gsea/'results/tables/Low_HR_selected_plotdata.csv',gsea/'results/qa/plot_scales.csv',gsea/'results/qa/HR_patient_rosters.csv']
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
cell,meta,mods,audit,g,scales,roster=[pd.read_csv(p) for p in sources]
cols=['Patient','Sample','treatment','histology','group','degradation_mean_z']
md=meta[cols].rename(columns={'group':'Degradation_group','degradation_mean_z':'Patient_Degradation_score'})
assert md.Patient.is_unique and len(md)==12
lineages=['Epithelial','CAF','B cell','Plasma cell','CD4 T cell','CD8 T cell','Cycling T cell','NK cell','Mast cell','Neutrophil','Monocyte','Macrophage','Conventional DC','pDC']
colors=['#F02C55','#D65720','#C28B00','#878600','#78B300','#009940','#00AA85','#009DC9','#0071D9','#3E36CA','#8153D4','#B14DD5','#AD158F','#EE54B1']
cmap=dict(zip(lineages,colors))
exports={}
for group,n in [('High',42328),('Low',34662)]:
 a=cell.loc[cell.group==group,['cell_id','Patient','final_lineage','umap_1','umap_2']].merge(md,on='Patient',validate='many_to_one')
 a=a.rename(columns={'cell_id':'Cell_ID','final_lineage':'Cell_type','umap_1':'UMAP_1','umap_2':'UMAP_2'})
 a=a[['Patient','Sample','Cell_ID','Degradation_group','Patient_Degradation_score','treatment','histology','Cell_type','UMAP_1','UMAP_2']]
 a['Cell_colour_hex']=a.Cell_type.map(cmap)
 # Direct Prism XY range: common UMAP_1 X and one Y column per lineage.
 for lin in lineages:a['Y_'+lin]=a.UMAP_2.where(a.Cell_type==lin)
 assert len(a)==n and a.Cell_ID.is_unique and (a.Degradation_group==group).all()
 assert (a[['Y_'+lin for lin in lineages]].notna().sum(axis=1)==1).all()
 exports[group+'.csv']=a
m=mods.drop(columns='group').merge(md,on='Patient',validate='many_to_one').rename(columns={'cell_id':'Cell_ID','final_lineage':'Cell_type','umap_1':'UMAP_1','umap_2':'UMAP_2','production_module_shifted':'Production_shifted_module_score','degradation_module_shifted':'Degradation_shifted_module_score'})
m=m[['Patient','Sample','Cell_ID','Degradation_group','Patient_Degradation_score','treatment','histology','Cell_type','UMAP_1','UMAP_2','Production_shifted_module_score','Degradation_shifted_module_score']]
for metric in ['Production','Degradation']:
 cap=audit.loc[audit.metric==metric,'cap'].unique();assert len(cap)==1
 m[metric+'_display_cap']=cap[0];m[metric+'_display_value']=np.minimum(m[metric+'_shifted_module_score'],cap[0])
 m[metric+'_draw_order']=m.groupby('Degradation_group')[metric+'_shifted_module_score'].rank(method='first').astype(int)
 # Rank ties must follow original pmin(display_value), cell_id plotting order.
 for gr in ['Low','High']:
  ix=m[m.Degradation_group==gr].sort_values([metric+'_display_value','Cell_ID']).index
  m.loc[ix,metric+'_draw_order']=np.arange(1,len(ix)+1)
assert len(m)==76990 and m.Cell_ID.is_unique
exports['Modules_Low_High.csv']=m
assert len(g)==26 and g.groupby('compartment').size().to_dict()=={'CD8':6,'Epithelial':8,'Whole_tumour':8,'cDC':4}
g=g.merge(scales,on='compartment',validate='many_to_one')
g['Y_plot_position']=g.groupby('compartment').pathway.transform('size')+1-g.display_order
g['Contrast']='Degradation High vs Low';g['Enriched_group']='Low';g['treatment']='Post'
for gr in ['High','Low']:
 ids=roster[roster.group==gr].groupby('compartment').Patient.agg(lambda z:';'.join(sorted(z)))
 g[gr+'_patient_IDs']=g.compartment.map(ids)
assert (g.NES<0).all() and np.allclose(-np.log10(g.BH_q),g.neglog10_BH_q)
assert all(len(t.split(';'))==n for t,n in zip(g.leading_edge,g.leading_edge_n))
keep=['compartment','display_order','Y_plot_position','display_label','pathway','collection','NES','ES','BH_q','neglog10_BH_q','leading_edge_n','gene_set_size','n','High','Low','High_patient_IDs','Low_patient_IDs','Contrast','Enriched_group','treatment','bh_family_n','q_source_column','colour_max','NES_min','NES_max','point_area_count_max','leading_edge','measured_genes']
exports['GSEA_Degradation_Low.csv']=g[keep]
checks=[]
for name,df in exports.items():
 target=out/name
 if target.exists():raise FileExistsError(target)
 df.to_csv(target,index=False,encoding='utf-8-sig')
 re=pd.read_csv(target)
 assert re.shape==df.shape
 for col in df:
  if pd.api.types.is_numeric_dtype(df[col]):assert np.allclose(df[col],re[col],rtol=1e-12,atol=1e-12,equal_nan=True),col
  else:assert df[col].fillna('').astype(str).equals(re[col].fillna('').astype(str)),col
 checks.append({'file':name,'rows':len(df),'columns':len(df.columns),'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
assert hashes=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
(mod/'qa/Four_CSV_export_checks.json').write_text(json.dumps({'sources':hashes,'outputs':checks,'readback_verified':True},indent=2))
(folder/'CSV_사용안내.md').write_text("""# 그림당 CSV 하나 — 2026-09-29

csv파일들 폴더의 4개 CSV는 첨부 그림 4개에 대응합니다. 환자군은 기존 whole-tumour Degradation Low5/High7, Post12입니다.

- High.csv: High UMAP, 42,328세포/7명. Low.csv: Low UMAP, 34,662세포/5명.
- Modules_Low_High.csv: 76,990세포를 한 행씩 저장. Degradation_group으로 Low/High를 구분하고 Production/Degradation score 두 열로 네 패널을 재현합니다. 환자별 Degradation score와 세포별 shifted module score는 다른 변수입니다.
- GSEA_Degradation_Low.csv: 첨부 GSEA 4패널을 compartment로 구분한 26경로. Whole_tumour8, Epithelial8, CD8 6, cDC4행. NES는 원래 High-minus-Low contrast의 음수 값입니다.

## Prism에서 사용할 열

High/Low: X=UMAP_1, Y=Y_Epithelial부터 Y_pDC까지의 14개 열을 XY dataset으로 사용하면 세포 유형별로 나뉩니다. 빈 Y칸은 해당 세포가 그 유형이 아니라는 의미이며 0으로 채우지 않습니다. UMAP_2는 모든 세포의 원래 Y좌표이며, Cell_colour_hex는 원래 색입니다. Patient/Sample/Cell_ID/분류 및 score는 분석·추적용입니다.

Modules: 그룹별 X=UMAP_1, Y=UMAP_2. 색의 수치값은 각 지표의 display_value이고 원래 수치는 shifted_module_score입니다. display_value=min(score, display_cap)이며 Production cap=0.2863583, Degradation cap=1.1007067(정확한 값은 CSV). score=0은 흰색, 지표별 두 군에 같은 척도입니다. draw_order 순서로 낮은 색값부터 그렸습니다. 원래 팔레트는 Production #FFFFFF/#FFF3DD/#FFBF55/#E56A0A/#943100, Degradation #FFFFFF/#E4F4FF/#65B9EB/#1776C6/#073B8C이며 위치는 0/0.18/0.4/0.7/1입니다.

GSEA: X=NES, Y=Y_plot_position, 라벨=display_label, 점 면적=leading_edge_n, 색=neglog10_BH_q. display_order는 위에서 아래 순서입니다. q는 원래 전체 검정군에 대한 BH 값이며 이 26행만으로 재계산하지 않습니다. 각 패널의 환자수·ID·BH family 크기를 함께 보존했습니다. q값에서 원래 p값을 역산하지 않았습니다.

CSV는 데이터를 보존하며 Prism 그래프 서식까지 자동 설정하는 파일은 아닙니다. 세포별 연속 색상이나 GSEA의 크기·색 동시 매핑은 사용 중인 Prism 버전에 맞는 추가 설정이 필요합니다. 이 컴퓨터에서 Prism 가져오기/작도 동작을 검증한 것은 아닙니다. GSEA 자체의 재실행은 경로별 결과 CSV만으로 할 수 없으며 원래 발현/랭킹·유전자집합·코드가 필요합니다. UMAP의 세포행을 환자 독립 표본으로 취급하지 않습니다.

원래 그림과 자료를 변경하지 않았습니다. 원본 대조, ID/행수/좌표/점수 및 CSV 재읽기 검증을 완료했습니다.
""")
print(json.dumps(checks,indent=2));print(out)
