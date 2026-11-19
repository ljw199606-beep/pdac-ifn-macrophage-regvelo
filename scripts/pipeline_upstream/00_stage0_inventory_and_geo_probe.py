#!/usr/bin/env python3
"""
Stage 0: Baseline Input Inventory, SHA256 Checksums, Original v0 vs v1 Manuscript Audit,
and Live NCBI GEO Provenance / GSE272362 Spatial Matrix Availability Probe.
Strictly read-only on all upstream directories.
"""
import os
import sys
import json
import hashlib
import urllib.request
import re
from pathlib import Path
from datetime import datetime
import pandas as pd

OUT_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo")
V0_ROOT = Path("/mnt/d/09-04/09-04")

AUDIT_DIR = OUT_ROOT / "01_audit"
META_DIR = OUT_ROOT / "02_metadata"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

def sha256_file(path: Path, max_bytes: int = None) -> str:
    h = hashlib.sha256()
    read_so_far = 0
    with open(path, "rb") as f:
        while True:
            chunk = f.read(8 * 1024 * 1024)
            if not chunk:
                break
            if max_bytes is not None and read_so_far + len(chunk) > max_bytes:
                h.update(chunk[: max_bytes - read_so_far])
                return h.hexdigest() + f"_partial_{max_bytes}B"
            h.update(chunk)
            read_so_far += len(chunk)
    return h.hexdigest()

print("[1/5] Building INPUT_MANIFEST.tsv with SHA256 checksums...")
manifest_rows = []

# Strictly read-only original v0 baseline files (excluding any cNMF/redesign directories)
v0_targets = [
    (V0_ROOT / "manuscript_pdac_regvelo.docx", "v0_baseline_manuscript_docx", "Original reviewed manuscript (v0)"),
    (V0_ROOT / "all_tables_pdac_regvelo.docx", "v0_baseline_tables_docx", "Original reviewed tables document (v0)"),
    (V0_ROOT / "manuscript_pdac_regvelo_rev2.pdf", "v0_baseline_pdf", "Original PDF snapshot (v0)"),
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

# Original v0 code scripts (read-only)
v0_code_dir = V0_ROOT / "原始代码" / "pdac_regvelo_project"
if v0_code_dir.exists():
    for p in sorted(v0_code_dir.rglob("*")):
        if p.is_file() and p.suffix in {".py", ".R", ".sh", ".md", ".yaml", ".yml", ".json", ".csv"} and p.stat().st_size < 50_000_000:
            st = p.stat()
            manifest_rows.append({
                "file_path": str(p),
                "file_type": f"v0_code_{p.suffix.lstrip('.')}",
                "version_date": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                "size_bytes": st.st_size,
                "sha256": sha256_file(p),
                "analysis_purpose": "Original v0 pipeline script/result for baseline issue verification",
            })

# v1_rebuild data, scripts, results, manuscript (read-only)
for sub, purpose in [
    ("02_data", "Primary raw/processed dataset or prior reference"),
    ("03_scripts", "v1 rebuild analysis script"),
    ("04_results", "v1 rebuild intermediate/summary output"),
    ("06_composite_figures", "v1 rebuild composite figure"),
    ("07_manuscript", "v1 rebuild manuscript/supplementary table"),
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
        # Full SHA256 for files <= 500MB, full SHA256 for h5ad too (fast NVMe, only ~6GB total)
        digest = sha256_file(p) if st.st_size <= 3_000_000_000 else sha256_file(p, 64 * 1024 * 1024)
        manifest_rows.append({
            "file_path": str(p),
            "file_type": p.suffix.lstrip(".") or "file",
            "version_date": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "size_bytes": st.st_size,
            "sha256": digest,
            "analysis_purpose": purpose,
        })

df_manifest = pd.DataFrame(manifest_rows)
df_manifest.to_csv(AUDIT_DIR / "INPUT_MANIFEST.tsv", sep="\t", index=False)
print(f"  Saved INPUT_MANIFEST.tsv with {len(df_manifest)} entries.")

print("[2/5] Extracting original v0 manuscript & all_tables_pdac_regvelo.docx text for baseline verification...")
try:
    import docx
    for name in ["manuscript_pdac_regvelo.docx", "all_tables_pdac_regvelo.docx"]:
        p = V0_ROOT / name
        if p.exists():
            doc = docx.Document(str(p))
            lines = []
            lines.append(f"=== PARAGRAPHS ({len(doc.paragraphs)}) ===")
            for i, para in enumerate(doc.paragraphs):
                if para.text.strip():
                    lines.append(f"[P{i}] {para.text.strip()}")
            lines.append(f"\n=== TABLES ({len(doc.tables)}) ===")
            for ti, tbl in enumerate(doc.tables):
                lines.append(f"\n--- Table {ti} ({len(tbl.rows)}x{len(tbl.columns)}) ---")
                for ri, row in enumerate(tbl.rows[:40]):
                    vals = [c.text.strip().replace("\n", " ") for c in row.cells]
                    lines.append(f"R{ri}: " + " | ".join(vals))
            out_txt = AUDIT_DIR / f"v0_extracted_{name.replace('.docx', '.txt')}"
            out_txt.write_text("\n".join(lines), encoding="utf-8")
            print(f"  Extracted {out_txt.name} ({len(lines)} lines)")
except Exception as e:
    print(f"  Warning reading docx: {e}")

print("[3/5] Probing live NCBI GEO records (GSE272362, GSE154778, GSE202051)...")
geo_summary = {}
for gse in ["GSE272362", "GSE154778", "GSE202051"]:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            txt = resp.read().decode("utf-8", errors="replace")
        (META_DIR / f"geo_{gse}_self_brief.soft").write_text(txt, encoding="utf-8")
        supps = [ln.split("=", 1)[1].strip() for ln in txt.splitlines() if ln.startswith("!Series_supplementary_file")]
        samples = [ln.split("=", 1)[1].strip() for ln in txt.splitlines() if ln.startswith("!Series_sample_id")]
        title = [ln.split("=", 1)[1].strip() for ln in txt.splitlines() if ln.startswith("!Series_title")]
        summary_ln = [ln.split("=", 1)[1].strip() for ln in txt.splitlines() if ln.startswith("!Series_summary")]
        overall = [ln.split("=", 1)[1].strip() for ln in txt.splitlines() if ln.startswith("!Series_overall_design")]
        geo_summary[gse] = {
            "title": title[0] if title else "",
            "n_samples": len(samples),
            "samples": samples,
            "supplementary_files": supps,
            "overall_design": " ".join(overall),
            "summary": " ".join(summary_ln),
        }
        print(f"  {gse}: {len(samples)} GSMs, {len(supps)} Series supplementary files")
        for s in supps:
            print(f"    supp: {s}")
    except Exception as e:
        print(f"  Error fetching {gse}: {e}")

# Also fetch GSM brief metadata for GSE272362 so we know exact tissue site, patient ID, and per-GSM supplementary files!
try:
    url272 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272362&targ=gsm&form=text&view=brief"
    req = urllib.request.Request(url272, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=45) as resp:
        gsm_txt = resp.read().decode("utf-8", errors="replace")
    (META_DIR / "geo_GSE272362_gsm_brief.soft").write_text(gsm_txt, encoding="utf-8")
    gsm_records = []
    cur = {}
    for ln in gsm_txt.splitlines():
        if ln.startswith("^SAMPLE ="):
            if cur:
                gsm_records.append(cur)
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
        gsm_records.append(cur)
    df_gsm272 = pd.DataFrame([
        {
            "gsm": r["gsm"],
            "title": r.get("title", ""),
            "source_name": r.get("source_name", ""),
            "characteristics": " | ".join(r.get("chars", [])),
            "n_supp_files": len(r.get("supps", [])),
            "supp_files": " ; ".join(r.get("supps", [])),
        }
        for r in gsm_records
    ])
    df_gsm272.to_csv(META_DIR / "gse272362_gsm_metadata.tsv", sep="\t", index=False)
    print(f"  Parsed {len(df_gsm272)} GSM records for GSE272362:")
    print(df_gsm272[["gsm", "title", "source_name", "characteristics", "n_supp_files"]].to_string())
except Exception as e:
    print(f"  Error fetching GSE272362 GSMs: {e}")

# Also fetch GSM brief metadata for GSE202051 and GSE154778
for gse in ["GSE202051", "GSE154778"]:
    try:
        url_gsm = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief"
        req = urllib.request.Request(url_gsm, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=45) as resp:
            gsm_txt = resp.read().decode("utf-8", errors="replace")
        (META_DIR / f"geo_{gse}_gsm_brief.soft").write_text(gsm_txt, encoding="utf-8")
        gsm_records = []
        cur = {}
        for ln in gsm_txt.splitlines():
            if ln.startswith("^SAMPLE ="):
                if cur:
                    gsm_records.append(cur)
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
            gsm_records.append(cur)
        df_gsm = pd.DataFrame([
            {
                "gsm": r["gsm"],
                "title": r.get("title", ""),
                "source_name": r.get("source_name", ""),
                "characteristics": " | ".join(r.get("chars", [])),
                "n_supp_files": len(r.get("supps", [])),
            }
            for r in gsm_records
        ])
        df_gsm.to_csv(META_DIR / f"{gse.lower()}_gsm_metadata.tsv", sep="\t", index=False)
        print(f"  Saved {gse.lower()}_gsm_metadata.tsv ({len(df_gsm)} GSMs)")
    except Exception as e:
        print(f"  Error fetching {gse} GSMs: {e}")

(META_DIR / "geo_series_probe_summary.json").write_text(json.dumps(geo_summary, indent=2), encoding="utf-8")
print("[Done Stage 0 probe Part 1]")
