"""Lossless PNG transport packaging and PDF text-boundary checks; no plot/data edits."""
from pathlib import Path
import hashlib, json, io, struct, subprocess
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image, ImageCms
OUT=Path(__file__).resolve().parent.parent
QA=OUT/"qa"
FILES=[("e","Fig3e_Pathway_Scores_7KO"),("i","Fig3i_Species_Function_7KO"),("j","Fig3j_Cross_Cancer_Concordance_7KO"),("k","Fig3k_Shared_Species_7KO")]
dest=OUT/"figures/PPT_insert";dest.mkdir(exist_ok=True)
profile=Path("/System/Library/ColorSync/Profiles/sRGB Profile.icc").read_bytes()
results=[];pdf=[]
for panel,stem in FILES:
    source=OUT/"figures"/(stem+".png"); target=dest/("PPT_INSERT_Fig3"+panel+"_7KO.png")
    with Image.open(source) as image:
        image.load()
        assert image.mode in ["RGB","RGBA"]
        if image.mode=="RGBA":assert np.asarray(image.getchannel("A")).min()==255
        rgb=image.convert("RGB"); pixels=np.asarray(rgb).copy()
    b=io.BytesIO();rgb.save(b,format="PNG",icc_profile=profile,dpi=(600,600))
    target.write_bytes(b.getvalue());target.chmod(0o644)
    with Image.open(target) as reread:
        reread.load();assert reread.mode=="RGB"
        assert np.array_equal(np.asarray(reread),pixels)
        assert reread.info["icc_profile"]==profile
        assert np.allclose(reread.info["dpi"],[600,600],atol=.001)
    subprocess.run(["/usr/bin/xattr","-c",str(target)],check=True)
    subprocess.run(["/usr/bin/chflags","nohidden",str(target)],check=True)
    attrs=subprocess.check_output(["/usr/bin/xattr",str(target)],text=True).splitlines()
    # macOS can retain its protected provenance tag; it does not prevent insertion.
    assert not (set(attrs)-{'com.apple.provenance'}),attrs
    width,height,depth,colour,comp,filt,interlace=struct.unpack(">IIBBBBB",target.read_bytes()[16:29])
    assert (depth,colour,interlace)==(8,2,0)
    results.append(dict(panel=panel,file=str(source),transport=str(target),width=width,height=height,dpi=600,
        mode="RGB",icc="sRGB",depth=8,interlace=0,xattrs=attrs,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        raw_rgb_sha256=hashlib.sha256(pixels.tobytes()).hexdigest(),pixels_changed=False))
    bbox=QA/(stem+"_PDF_bbox.html")
    subprocess.run(["/opt/homebrew/bin/pdftotext","-bbox",str(OUT/"figures"/(stem+".pdf")),str(bbox)],check=True)
    xml=ET.parse(bbox);ns={"x":"http://www.w3.org/1999/xhtml"};pages=xml.findall(".//x:page",ns)
    assert len(pages)==1
    page=pages[0];w=float(page.attrib["width"]);h=float(page.attrib["height"])
    words=page.findall("x:word",ns);bad=[]
    for word in words:
        a={k:float(v) for k,v in word.attrib.items()}
        if a["xMin"]<0 or a["yMin"]<0 or a["xMax"]>w or a["yMax"]>h:bad.append(word.text)
    assert words and not bad,(stem,bad)
    # Render the actual vector PDF for visual glyph/connector inspection.
    subprocess.run(["/opt/homebrew/bin/pdftoppm","-png","-r","150","-singlefile",
        str(OUT/"figures"/(stem+".pdf")),str(QA/(stem+"_PDF_render"))],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    pdf.append(dict(panel=panel,pages=1,words=len(words),out_of_bounds=bad,width_pt=w,height_pt=h))
(QA/"raster_validation.json").write_text(json.dumps(dict(status="PASS",files=results,pdf=pdf),indent=2)+"\n")
print("PASS: four lossless RGB/sRGB 600-dpi PNGs; four vector PDFs with all text in bounds.")
