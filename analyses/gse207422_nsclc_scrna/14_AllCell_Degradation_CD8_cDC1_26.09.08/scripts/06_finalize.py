"""Final read-only source/legacy recheck, syntax checks and derived-artifact manifest.
Does not rerun or modify inferential results.
"""
from pathlib import Path
import ast
import csv
import hashlib
import subprocess

out = Path(__file__).resolve().parents[1]
src = out.parent
tables = out / "results/tables"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

checks = []
for filename, base, key in [
    ("00_inputs.csv", src, "path"),
    ("00_protected_before.csv", src, "path"),
    ("01_exposure_lock.csv", out, "file"),
]:
    entries = rows(tables / filename)
    for row in entries:
        p = base / row[key]
        actual = sha(p)
        checks.append(dict(family=filename, path=row[key], pass_check=actual == row["sha256"],
                           expected_sha256=row["sha256"], actual_sha256=actual))
        assert actual == row["sha256"], f"Protected file changed: {p}"

for p in sorted((out / "scripts").glob("*.py")):
    ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
r_files = sorted((out / "scripts").glob("*.R"))
expr = "for (p in commandArgs(TRUE)) parse(file=p); cat('R syntax PASS\\n')"
subprocess.run(["/usr/local/bin/Rscript", "-e", expr, *map(str, r_files)], check=True)

for filename, expected in [("06_verification_checks.csv", 118), ("08_artifact_checks.csv", 88)]:
    actual = len(rows(tables / filename))
    assert actual == expected, (filename, actual)
summary = rows(tables / "07_primary_cross_method_summary.csv")
assert len(summary) == 8
sig = {(r["cell_type"], r["label"]) for r in summary
       if float(r["q_global"]) < .05 and float(r["NES"]) > 0}
assert sig == {("CD8", "APC cross-presentation"), ("CD8", "TCR signaling"),
               ("CD8", "IFNG response"), ("cDC1", "IFNG response")}
assert all(float(r["CAMERA_q"]) > .05 and float(r["score_q"]) > .05 for r in summary)
figs = sorted((out / "results/figures").iterdir())
assert len(figs) == 18
assert all(sum(p.suffix == ext for p in figs) == 6 for ext in (".png", ".pdf", ".svg"))

with (tables / "09_final_preservation.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(checks[0]))
    writer.writeheader()
    writer.writerows(checks)
qa = out / "results/qa/final_integrity_summary.txt"
qa.write_text(
    f"PASS\nProtected hashes rechecked: {len(checks)}\n"
    f"R scripts parsed: {len(r_files)}\nPython scripts parsed: 2\n"
    "Primary exact-set interpretation and 6 PNG/PDF/SVG triplets verified.\n"
    "No inferential analysis rerun.\n", encoding="utf-8")
manifest = tables / "artifact_manifest.csv"
files = [p for p in sorted(out.rglob("*")) if p.is_file() and p != manifest
         and "font_cache" not in p.parts and p.name != "06_finalize.log"
         and "__pycache__" not in p.parts]
with manifest.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["path", "bytes", "sha256"])
    writer.writeheader()
    writer.writerows(dict(path=str(p.relative_to(out)), bytes=p.stat().st_size,
                          sha256=sha(p)) for p in files)
print(qa.read_text(), f"Manifest: {len(files)} files", sep="\n")

