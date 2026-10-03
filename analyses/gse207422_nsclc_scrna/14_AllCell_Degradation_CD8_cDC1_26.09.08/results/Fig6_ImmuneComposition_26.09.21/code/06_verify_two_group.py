from pathlib import Path
from PIL import Image,ImageCms
import csv,json,hashlib,subprocess
b=Path(__file__).resolve().parent.parent
v=list(csv.DictReader((b/'tables/patient_values.csv').open()));g=list(csv.DictReader((b/'tables/two_group_heatmap.csv').open(encoding='utf-8-sig')));s=list(csv.DictReader((b/'tables/statistics.csv').open()))
checks=0
for r in g:
 d=[x for x in v if x['cell_type']==r['cell_type'] and x['group']==r['group']];st=next(x for x in s if x['cell_type']==r['cell_type'])
 assert len(d)==int(r['n'])
 assert abs(sum(float(x['percent']) for x in d)/len(d)-float(r['mean_percent']))<1e-10
 assert abs(sum(float(x['row_z']) for x in d)/len(d)-float(r['mean_patient_z']))<1e-12
 assert float(r['BH13_q'])==float(st['BH13_q'])
 checks+=4
p=b/'figures/H2.png'
with Image.open(p) as im:
 im.load();im.convert('RGB').save(p,dpi=(600,600),icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes())
with Image.open(p) as im:
 assert im.mode=='RGB' and im.size==(5100,5160) and im.info['icc_profile'] and abs(im.info['dpi'][0]-600)<.1
subprocess.run(['/opt/homebrew/bin/pdftoppm','-scale-to','1500','-png','-singlefile',str(b/'figures/H2.pdf'),str(b/'qa/H2_PDF')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
qa=b/'qa/H2_QA.json';q=json.loads(qa.read_text())
for it in q['inputs']:
 assert hashlib.sha256(Path(it['path']).read_bytes()).hexdigest()==it['sha256']
q.update(independent_value_checks=checks,PNG='5100x5160 RGB8/sRGB 600dpi',PNG_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),visual_QA='PNG checked; PDF preview pending')
qa.write_text(json.dumps(q,ensure_ascii=False,indent=2))
print('Passed',checks,'independent numeric checks; source hashes retained; RGB8/sRGB600dpi verified.')
