from pathlib import Path
import csv,json,hashlib,sys
cd8_only = "--cd8-only" in sys.argv
from datetime import datetime
from zoneinfo import ZoneInfo
base=Path(__file__).resolve().parent.parent
origin=base.parents[1];src=origin.parent
qa=base/'results/qa';tables=base/'results/tables'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
counts={}
for filename,column,label in [('06_verification_checks.csv','passed','existing_analysis_independent_checks'),('02_selection_checks.csv','pass','additional_selection_checks')]:
 rr=read(tables/filename);assert all(r[column]=='TRUE' for r in rr);counts[label]=len(rr)
a=json.loads((qa/('artifact_checks_CD8_repair.json' if cd8_only else 'artifact_checks.json')).read_text());assert all(r['pass_check'] for r in a);counts['numpy_and_artifact_checks']=len(a)
a=json.loads((qa/('ppt_transport_checks_CD8_repair.json' if cd8_only else 'ppt_transport_checks.json')).read_text());assert all(a.values());counts['PPT_checks']=len(a)
ports=json.loads((qa/('transport_manifest_CD8_repair.json' if cd8_only else 'transport_manifest.json')).read_text());lib={r['file']:r for r in read(qa/'R_libpng_pixel_checks.csv')}
for r in ports:
 p=Path(r['transport']);assert sha(p)==r['sha256'];assert r['pixel_sha256']==lib[p.name]['RGB_pixel_sha256']
counts['R_libpng_pixel_checks']=len(ports)
protected=read(tables/'00_protected_before.csv')
for r in protected:
 p=src/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes']),r['path']
counts['existing14_15_files_unchanged']=len(protected)
for filename,label in [('06_input_preservation.csv','source_inputs_verified'),('06_previous_artifact_preservation.csv','earlier09_13_files_verified')]:
 rr=read(tables/filename);assert all(r['unchanged']=='TRUE' for r in rr);counts[label]=len(rr)
assert sha(base/'ANALYSIS_DISPLAY_PLAN.md')==(base/'logs/display_plan_sha256.txt').read_text().strip()
result={'time_KST':datetime.now(ZoneInfo('Asia/Seoul')).strftime('%Y-%m-%d %H:%M KST'),'all_pass':True,'checks':counts,'inference_rerun_this_task':False,'author_PPT_DOCX_edited':False,'native_PowerPoint_tested':False,'PPT_saved_reimported_rendered_with':'@oai/artifact-tool'}
if cd8_only:
 result.update(task_scope='CD8-only in-place redraw of two author-selected figures, 2026-09-15', native_PowerPoint_tested=False, native_PowerPoint_attempted=True, native_insertion_and_Finder_drag_verified=False, native_blocker='CUA ScreenCaptureKit -3811 capture failure after scratch presentation creation; author confirms computer and external monitor operating normally', drag_error_cause='Unresolved; original PNG files decode normally. Original Preview/quarantine xattrs observed and cleared; causality not established.', current_transport_count=2, statistics_changed=False, original_q_family=5959, visual_PDF_and_reopened_PPTX_review=True)
 before=json.loads((qa/'before_CD8_repair.json').read_text())
 protected_tables=[r for r in before if '/results/tables/' in r['path'] and not r['path'].endswith('/artifact_manifest.csv')]
 for r in protected_tables: assert sha(Path(r['path']))==r['sha256'],r['path']
 result['unchanged_display_and_analysis_tables']=len(protected_tables)
 result['native_test_remaining']=True
 if (qa/'short_path_repair.json').exists():
  repair=json.loads((qa/'short_path_repair.json').read_text())
  assert all(repair['checks'].values())
  result.update(task_scope='CD8 transport filename shortening after author-reported continued drag failure', current_transport_names=['G.png','I.png'], current_name_checks=len(repair['checks']), author_drag_status=repair['author_drag_status'], native_test_remaining=False, native_UI_required=False, drag_issue_resolved=False, drag_error_cause='Full-path length is a supported hypothesis; root cause and short-name drag success not yet author-confirmed')
  if repair.get('author_drag_status') == 'confirmed_success':
   result.update(drag_issue_resolved=True, native_drag_confirmed_by='author', native_drag_confirmation_time_KST=repair.get('author_confirmation_time_KST'), drag_error_cause='Author confirms success after same-folder filename/full-path shortening; exact internal limit not established', native_blocker=None)
(qa/'final_verification_summary.json').write_text(json.dumps(result,indent=2)+'\n')
manifest=tables/'artifact_manifest.csv'
files=sorted(p for p in base.rglob('*') if p.is_file() and p!=manifest and p.name!='.DS_Store' and 'font_cache' not in p.parts)
with manifest.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);w.writeheader()
 for p in files:w.writerow({'path':str(p.relative_to(base)),'bytes':p.stat().st_size,'sha256':sha(p)})
rr=read(manifest)
for r in rr:
 p=base/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
print(json.dumps(result,indent=2));print('MANIFEST',len(rr),'files',sum(int(r['bytes']) for r in rr),'bytes',sha(manifest))
