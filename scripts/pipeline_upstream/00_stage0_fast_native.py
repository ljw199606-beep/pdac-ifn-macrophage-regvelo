#!/usr/bin/env python3
"""
Stage 0 (Native Windows NVMe Speed):
1. Build INPUT_MANIFEST.tsv with full SHA256 checksums across v0 baseline (read-only) and v1 data/scripts/results (read-only).
2. Extract text & tables from original v0 manuscript_pdac_regvelo.docx & all_tables_pdac_regvelo.docx for exact line-by-line baseline verification.
3. Save official NCBI GEO SOFT metadata for GSE272362, GSE154778, and GSE202051.
"""
import os
import json
import hashlib
import subprocess
from pathlib import Path
from datetime import datetime
import pandas as pd

OUT_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo")
V0_ROOT = Path(r"D:\09-04\09-04")

AUDIT_DIR = OUT_ROOT / "01_audit"
META_DIR = OUT_ROOT / "02_metadata"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(16 * 1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

print("[1/3] Computing native NVMe SHA256 checksums for INPUT_MANIFEST.tsv...", flush=True)
manifest_rows = []

v0_targets = [
    (V0_ROOT / "manuscript_pdac_regvelo.docx", "v0_baseline_manuscript_docx", "Original reviewed manuscript (v0 baseline)"),
    (V0_ROOT / "all_tables_pdac_regvelo.docx", "v0_baseline_tables_docx", "Original reviewed tables document (v0 baseline)"),
    (V0_ROOT / "manuscript_pdac_regvelo_rev2.pdf", "v0_baseline_pdf", "Original PDF snapshot (v0 baseline)"),
]
for p, ftype, usage in v0_targets:
    if p.exists():
        st = p.stat()
        manifest_rows.append({
            "file_path": str(p),
            "file_type": ftype,
            "version_date": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "size_bytes": st.st_size,
            "sha256": sha256_file(p),
            "analysis_purpose": usage,
        })

v0_code_dir = V0_ROOT / "原始代码" / "pdac_regvelo_project"
if v0_code_dir.exists():
    for p in sorted(v0_code_dir.rglob("*")):
        if p.is_file() and p.suffix in {".py", ".R", ".sh", ".md", ".json", ".csv"} and p.stat().st_size < 50_000_000:
            st = p.stat()
            manifest_rows.append({
                "file_path": str(p),
                "file_type": f"v0_code_{p.suffix.lstrip('.')}",
                "version_date": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                "size_bytes": st.st_size,
                "sha256": sha256_file(p),
                "analysis_purpose": "Original v0 pipeline script/output for baseline issue audit",
            })

for sub, purpose in [
    ("02_data", "Primary dataset / count matrix / prior network"),
    ("03_scripts", "v1 rebuild script"),
    ("04_results", "v1 rebuild intermediate/summary result"),
    ("06_composite_figures", "v1 rebuild figure"),
    ("07_manuscript", "v1 rebuild manuscript / supplementary table"),
]:
    base = V1_ROOT / sub
    if not base.exists():
        continue
    for p in sorted(base.rglob("*")):
        if not p.is_file():
            continue
        if "starsolo_out" in str(p) and p.name not in {"Summary.csv"}:
            continue
        if "raw" in p.parts and p.suffix in {".fastq", ".gz", ".sra", ".bam"}:
            continue
        st = p.stat()
        manifest_rows.append({
            "file_path": str(p),
            "file_type": p.suffix.lstrip(".") or "file",
            "version_date": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "size_bytes": st.st_size,
            "sha256": sha256_file(p),
            "analysis_purpose": purpose,
        })

df_manifest = pd.DataFrame(manifest_rows)
df_manifest.to_csv(AUDIT_DIR / "INPUT_MANIFEST.tsv", sep="\t", index=False)
print(f"  Saved INPUT_MANIFEST.tsv ({len(df_manifest)} files hashed).", flush=True)

print("[2/3] Extracting v0 baseline docx files...", flush=True)
import docx
for name in ["manuscript_pdac_regvelo.docx", "all_tables_pdac_regvelo.docx"]:
    p = V0_ROOT / name
    if p.exists():
        doc = docx.Document(str(p))
        lines = [f"=== PARAGRAPHS ({len(doc.paragraphs)}) ==="]
        for i, para in enumerate(doc.paragraphs):
            if para.text.strip():
                lines.append(f"[P{i}] {para.text.strip()}")
        lines.append(f"\n=== TABLES ({len(doc.tables)}) ===")
        for ti, tbl in enumerate(doc.tables):
            lines.append(f"\n--- Table {ti} ({len(tbl.rows)}x{len(tbl.columns)}) ---")
            for ri, row in enumerate(tbl.rows[:60]):
                vals = [c.text.strip().replace("\n", " ") for c in row.cells]
                lines.append(f"R{ri}: " + " | ".join(vals))
        out_txt = AUDIT_DIR / f"v0_extracted_{name.replace('.docx', '.txt')}"
        out_txt.write_text("\n".join(lines), encoding="utf-8")
        print(f"  Saved {out_txt.name} ({len(lines)} lines)", flush=True)

print("[3/3] Fetching GEO Series + GSM metadata via curl.exe...", flush=True)
for gse in ["GSE272362", "GSE154778", "GSE202051"]:
    for targ in ["self", "gsm"]:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ={targ}&form=text&view=brief"
        out = subprocess.check_output(["curl.exe", "-sL", url]).decode("utf-8", "replace")
        (META_DIR / f"geo_{gse}_{targ}_brief.soft").write_text(out, encoding="utf-8")
        if targ == "gsm":
            recs, cur = [], {}
            for ln in out.splitlines():
                if ln.startswith("^SAMPLE ="):
                    if cur:
                        recs.append(cur)
                    cur = {"gsm": ln.split("=", 1)[1].strip(), "chars": [], "supps": []}
                elif ln.startswith("!Sample_title ="):
                    cur["title"] = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_source_name_ch1 ="):
                    cur["source_name"] = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_characteristics_ch1 ="):
                    cur["chars"].append(ln.split("=", 1)[1].strip())
                elif ln.startswith("!Sample_supplementary_file"):
                    cur["supps"].append(ln.split("=", 1)[1].strip())
            if cur:
                recs.append(cur)
            pd.DataFrame([
                {
                    "gsm": r["gsm"],
                    "title": r.get("title", ""),
                    "source_name": r.get("source_name", ""),
                    "characteristics": " | ".join(r.get("chars", [])),
                    "n_supp_files": len(r.get("supps", [])),
                    "supp_files": " ; ".join(r.get("supps", [])),
                }
                for r in recs
            ]).to_csv(META_DIR / f"{gse.lower()}_gsm_metadata.tsv", sep="\t", index=False)
            print(f"  Saved {gse.lower()}_gsm_metadata.tsv ({len(recs)} GSMs)", flush=True)

print("[Stage 0 Native Inventory Complete]", flush=True)
