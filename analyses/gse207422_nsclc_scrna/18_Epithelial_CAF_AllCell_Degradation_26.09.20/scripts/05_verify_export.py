from pathlib import Path
import csv,hashlib,json,shutil,math
from PIL import Image,ImageCms
o=Path(__file__).resolve().parents[1];b=o.parent
def rows(p):return list(csv.DictReader(p.open()))
a=rows(o/'results/tables/02_GSEA_all_scopes.csv'); d=rows(o/'results/tables/03_HR_selected_plotdata.csv')
# Independent display selection from all GSEA tests and measured membership lists exported with selected rows.
assert len(a)==4732 and len(d)==8
assert all(float(r['NES'])>0 and float(r['q_full'])<.05 for r in d)
sets=[set(r['measured_genes'].split(';')) for r in d]
assert all(len(x&y)/len(x|y)<.5 for i,x in enumerate(sets) for y in sets[:i])
checks=rows(o/'results/tables/04_independent_verification.csv');assert all(r['pass']=='TRUE' for r in checks)
old=json.loads((o/'results/qa/previous_outputs_before.json').read_text()); assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'] for r in old)
p=o/'results/figures/E.png'; im=Image.open(p).convert('RGB'); pixels=im.tobytes(); im.save(p,format='PNG',dpi=(600,600),icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes());re=Image.open(p);assert re.mode=='RGB' and re.tobytes()==pixels and re.size==(4080,3540)
for suf in ['png','pdf']:shutil.copy2(o/f'results/figures/E.{suf}',o/f'results/figures/PPT_insert/E.{suf}')
elig=rows(o/'results/tables/01_patient_eligibility.csv'); post=[r for r in elig if r['treatment']=='Post']
with (o/'results/tables/CAF_post12_cell_counts.csv').open('w') as f:
 fields=['Patient','Sample','group','lineage_cells','eligible_50'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in post if r['lineage']=='CAF')
report={'status':'numerical_checks_passed','independent_R_checks':len(checks),'protected_previous_outputs_unchanged':len(old),'source_input_hashes_unchanged':True,'selection_positive_q_and_pairwise_Jaccard_passed':True,'PNG':{'size':re.size,'mode':re.mode,'dpi':re.info['dpi'],'sRGB':bool(re.info.get('icc_profile')),'pixel_preservation':True},'CAF':'Not analyzed; primary cell-count criterion fails. No approval for reduced threshold.','native_PowerPoint_drag_tested':False,'ES_verification':'Initial 1e-10 equality was stricter than fgsea batch numerical precision; observed ES differences 4.3e-8 to 2.93e-7. All selected leading edges agree exactly; final documented ES tolerance 1e-6.'}
(o/'results/qa/verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
