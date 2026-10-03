#!/usr/bin/env python3
"""Read-only figure inspection plus derived QA renders; never edits analysis results."""
from pathlib import Path
import csv, hashlib, json, subprocess, xml.etree.ElementTree as ET
import pdfplumber
from PIL import Image
root=Path(__file__).resolve().parents[1]
fig=root/"results/figures"; dest=root/"results/qa"
checks=[]
def add(name,ok,detail=""):
    checks.append({"check":name,"passed":bool(ok),"detail":str(detail)})
    if not ok: raise AssertionError(f"{name}: {detail}")
def rows(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))
bases=sorted(p.stem for p in fig.glob("*.png"))
add("Six complete figure triplets",len(bases)==6 and all((fig/(b+e)).is_file() for b in bases for e in [".png",".pdf",".svg"]))
texts={}
for b in bases:
    png=fig/(b+".png"); pdf=fig/(b+".pdf"); svg=fig/(b+".svg")
    with Image.open(png) as im:
        im.load(); add(b+" PNG dimensions",im.width>=4000 and im.height>=2000,f"{im.width}x{im.height}")
        add(b+" PNG normal color",im.mode in ["RGB","RGBA"],im.mode)
        if im.mode=="RGBA": add(b+" PNG opaque alpha",im.getchannel("A").getextrema()==(255,255))
        dpi=im.info.get("dpi",(0,0));add(b+" PNG400dpi",all(399<v<401 for v in dpi),dpi)
        size=im.size
    with pdfplumber.open(pdf) as doc:
        add(b+" PDF one page",len(doc.pages)==1)
        page=doc.pages[0]; text=page.extract_text() or "";texts[b]=text
        # Cairo PDF rounds the media box to integer points; permit one point, not rescaling.
        add(b+" PDF canvas matches PNG",abs(size[0]/400*72-page.width)<=1 and abs(size[1]/400*72-page.height)<=1,(size,page.width,page.height))
        add(b+" PDF readable text",len(text)>200 and "\ufffd" not in text)
        bad=[c for c in page.chars if c.get("text","").strip() and (c["x0"]<-.5 or c["x1"]>page.width+.5 or c["top"]<-.5 or c["bottom"]>page.height+.5)]
        add(b+" PDF text within canvas",not bad,[(x["text"],x["x0"],x["top"]) for x in bad[:10]])
        add(b+" explicit group contrast","High" in text and "Low" in text)
    tree=ET.parse(svg);add(b+" SVG valid",tree.getroot().tag.endswith("svg"))
    add(b+" SVG no external asset",all(not(str(v).startswith(("http://","https://","file:"))) for e in tree.iter() for k,v in e.attrib.items() if k.endswith("href")))
    subprocess.run(["/opt/homebrew/bin/pdftoppm","-r","110","-singlefile","-png",str(pdf),str(dest/("pdf_"+b))],check=True,capture_output=True)
    with Image.open(dest/("pdf_"+b+".png")) as im:
        im.load();add(b+" PDF raster opened",im.width>1000 and im.height>600)
cd8=texts["Fig_02_CD8_AllCell_Groups_Exact_GSEA"]
dc=texts["Fig_03_cDC1_AllCell_Groups_Exact_GSEA"]
add("CD8 primary n labels","High 7 / Low 5" in cd8 and "16,229" in cd8)
add("cDC1 primary n labels","High 6 / Low 4" in dc and "111 cells" in dc)
add("Score figure both cell row labels",all(x in texts["Fig_05_CD8_cDC1_Primary_Patient_Pathway_Scores"] for x in ["CD8 T cells","cDC1-like cells"]))
tab=root/"results/tables"
allfocus=rows(tab/"05_focus_complete.csv")
primary=[r for r in allfocus if r["scope"]=="post_group" and r["role"]=="Fig6 exact"]
for ct,text in [("CD8",cd8),("cDC1",dc)]:
    for r in [r for r in primary if r["cell_type"]==ct]:
        q=float(r["q_global"]); f=f"{q:.1e}" if q<.001 else f"{q:.3f}"
        add(ct+" q text "+r["pathway"],f in text,f)
        add(ct+" NES text "+r["pathway"],f"{float(r['NES']):.2f}" in text)
    p=rows(tab/("plotdata_"+ct+"_primary_exact.csv"))
    left={(r["cell_type"],r["pathway"]):(float(r["NES"]),float(r["q_global"])) for r in p}
    right={(r["cell_type"],r["pathway"]):(float(r["NES"]),float(r["q_global"])) for r in primary if r["cell_type"]==ct}
    add(ct+" plot source equals results",left==right)
with (tab/"08_artifact_checks.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["check","passed","detail"]);w.writeheader();w.writerows(checks)
(dest/"artifact_qa_summary.txt").write_text(f"PASS\n{len(checks)} checks\nPDFs rendered for human/model visual review; programmatic bounds checks alone do not certify layout.\n")
(dest/"pdf_extracted_text.json").write_text(json.dumps(texts,ensure_ascii=False,indent=2))
print(f"Artifact checks PASS: {len(checks)}. Review six pdf_*.png renders before final delivery.")
