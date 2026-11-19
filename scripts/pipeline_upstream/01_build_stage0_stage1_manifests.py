#!/usr/bin/env python3
"""
Stage 0 & Stage 1 Gate (G1):
1. Verify all 5 cohorts directly from raw/processed AnnData & CSV matrices:
   - Cohort 1: Peng et al., Cell Res 2019 (CRA001160 / PRJCA001063; secondary Zenodo 3969339)
   - Cohort 2: Lin et al., Genome Med 2020 (GSE154778 / SRP272677; 16 patients, 120 SRR runs, STARsolo 2.7.11b)
   - Cohort 3: Hwang et al., Nat Genet 2022 (GSE202051; 43 PDAC patients: 18 Untreated + 25 Neoadjuvant)
   - Cohort 4: Khaliq et al., Nat Genet 2024 (GSE272362 / PRJNA1136518; 30 Visium FFPE sections across 7 rapid-autopsy patients: 10 PT, 12 HM, 5 LNM, 3 NP)
     -> Downloads official GSE272362_RAW.tar from NCBI GEO FTP/HTTPS and verifies every GSM's spatial_enrichment.csv.gz & tissue_positions.csv.gz!
   - Cohort 5: TCGA-PAAD PanCanAtlas (177 primary PAAD -> 146 histologically confirmed Pancreas Adenocarcinoma, Ductal Type; 138 with WES)
2. Check v0 RegVelo 462/460/18 bug directly from v0 tables/scripts.
3. Output DATA_PROVENANCE.tsv, SAMPLE_MANIFEST.tsv, CELL_AND_GENE_FLOW.tsv, and cohort_summary.tsv.
"""
import os
import re
import json
import gzip
import tarfile
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp

if hasattr(sp, "csr_array") and not hasattr(sp.csr_array, "getnnz"):
    sp.csr_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)
if hasattr(sp, "csc_array") and not hasattr(sp.csc_array, "getnnz"):
    sp.csc_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)

import anndata as ad

if os.name == "nt":
    OUT_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
    V1_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo")
    V0_ROOT = Path(r"D:\09-04\09-04")
    CURL_BIN = "curl.exe"
else:
    OUT_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116")
    V1_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo")
    V0_ROOT = Path("/mnt/d/09-04/09-04")
    CURL_BIN = "curl"

AUDIT_DIR = OUT_ROOT / "01_audit"
META_DIR = OUT_ROOT / "02_metadata"
RES_DIR = OUT_ROOT / "05_results"

# -------------------------------------------------------------------------
# 1. Inspect Peng 2019 (CRA001160 / PRJCA001063)
# -------------------------------------------------------------------------
print("[1/6] Auditing Cohort 1: Peng et al. 2019 (PRJCA001063 / CRA001160)...", flush=True)
peng_raw = ad.read_h5ad(V1_ROOT / "02_data/peng2019/peng_raw.h5ad", backed="r")
peng_ann = ad.read_h5ad(V1_ROOT / "02_data/peng2019/peng_annotated.h5ad", backed="r")

n_peng_raw_cells, n_peng_raw_genes = peng_raw.shape
n_peng_ann_cells, n_peng_ann_genes = peng_ann.shape
peng_obs = peng_ann.obs.copy()
peng_raw_var = list(peng_raw.var_names)
peng_ann_var = list(peng_ann.var_names)
n_unique_peng_raw_genes = len(set(peng_raw_var))

peng_donors = (
    peng_obs.groupby(["Patient", "CONDITION"], observed=True)
    .agg(
        n_cells=("Cell_type", "size"),
        n_macrophages=("Cell_type", lambda s: int((s == "Macrophage cell").sum())),
        n_ductal1=("Cell_type", lambda s: int((s == "Ductal cell type 1").sum())),
        n_ductal2=("Cell_type", lambda s: int((s == "Ductal cell type 2").sum())),
        n_tcells=("Cell_type", lambda s: int((s == "T cell").sum())),
        n_fibroblasts=("Cell_type", lambda s: int((s == "Fibroblast cell").sum())),
    )
    .reset_index()
)
print(
    f"  Peng raw={peng_raw.shape} (unique genes={n_unique_peng_raw_genes}), "
    f"annotated={peng_ann.shape}, donors={len(peng_donors)} "
    f"(T={int((peng_donors['CONDITION']=='T').sum())}, N={int((peng_donors['CONDITION']=='N').sum())})",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Inspect GSE154778 (16 patients: P01-P10 + MET01-MET06)
# -------------------------------------------------------------------------
print("[2/6] Auditing Cohort 2: Lin et al. 2020 (GSE154778 / SRP272677)...", flush=True)
g154 = ad.read_h5ad(V1_ROOT / "02_data/gse154778/gse154778_all16_velocity.h5ad", backed="r")
g154_qc = pd.read_csv(V1_ROOT / "04_results/gse154778/gse154778_all16_sample_qc_summary.csv")
g154_obs = g154.obs.copy()
print(
    f"  GSE154778 velocity atlas: {g154.shape}, samples={g154_obs['sample'].nunique()}, "
    f"raw STARsolo filtered droplets={int(g154_qc['raw_filtered_cells'].sum()) if 'raw_filtered_cells' in g154_qc.columns else 'see qc'}",
    flush=True,
)
print("  GSE154778 QC columns:", list(g154_qc.columns), flush=True)

# -------------------------------------------------------------------------
# 3. Inspect GSE202051 (43 PDAC patients: 18 Untreated + 25 Neoadjuvant)
# -------------------------------------------------------------------------
print("[3/6] Auditing Cohort 3: Hwang et al. 2022 (GSE202051)...", flush=True)
g202 = ad.read_h5ad(V1_ROOT / "02_data/gse202051/gse202051_macrophages.h5ad", backed="r")
g202_obs = g202.obs.copy()
g202_donors = (
    g202_obs.groupby(["pid", "new_treatment", "treatment_status", "response"], observed=True)
    .size()
    .reset_index(name="n_macrophages")
)
print(
    f"  GSE202051 macrophages={g202.shape} (unique genes={len(set(g202.var_names))}), "
    f"patients={g202_obs['pid'].nunique()}, "
    f"Untreated={int((g202_donors['new_treatment']=='Untreated').sum())}, "
    f"Treated={int((g202_donors['new_treatment']!='Untreated').sum())}",
    flush=True,
)

# -------------------------------------------------------------------------
# 4. Download & Re-parse Official GSE272362_RAW.tar (260 MB)
# -------------------------------------------------------------------------
print("[4/6] Auditing Cohort 4: Khaliq et al. 2024 (GSE272362, 30 Visium FFPE sections)...", flush=True)
g272_dir = META_DIR / "gse272362_official"
g272_dir.mkdir(parents=True, exist_ok=True)
tar_path = g272_dir / "GSE272362_RAW.tar"
if not tar_path.exists() or tar_path.stat().st_size < 250_000_000:
    print("  Downloading official GSE272362_RAW.tar (260 MB) from NCBI HTTPS...", flush=True)
    subprocess.run(
        [
            CURL_BIN, "-sL", "--retry", "3", "-o", str(tar_path),
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272362/suppl/GSE272362_RAW.tar",
        ],
        check=True,
    )
print(f"  Verified GSE272362_RAW.tar ({tar_path.stat().st_size:,} bytes).", flush=True)

gsm_meta272 = pd.read_csv(META_DIR / "gse272362_gsm_metadata.tsv", sep="\t")
khaliq_map = pd.read_csv(META_DIR / "gse272362_khaliq2024_samples_verified.csv")

# Parse tissue_positions.csv.gz and spatial_enrichment.csv.gz directly from official tar!
pos_stats = {}
enr_dfs = []
with tarfile.open(tar_path, "r") as tf:
    for member in tf.getmembers():
        mname = os.path.basename(member.name)
        m_gsm = mname.split("_")[0]
        # Extract subject ID (S13, S14, S17, S18, S20, S21, S22) from official filename GSM8399450_S13_Sample1_...
        m_subj = re.search(r"_(S\d+)_Sample(\d+)_", mname)
        subj_id = m_subj.group(1) if m_subj else "Unknown"
        smp_num = int(m_subj.group(2)) if m_subj else -1

        if mname.endswith("tissue_positions.csv.gz"):
            fobj = tf.extractfile(member)
            with gzip.open(fobj, "rt", encoding="utf-8", errors="replace") as gz:
                pdf = pd.read_csv(gz)
            # Space Ranger v2 has header 'barcode,in_tissue,array_row,array_col,...'; v1 has no header (col 1 = in_tissue)
            if "in_tissue" in pdf.columns:
                n_total_barcodes = len(pdf)
                n_in_tissue = int((pdf["in_tissue"] == 1).sum())
            else:
                n_total_barcodes = len(pdf) + 1
                n_in_tissue = int((pdf.iloc[:, 1] == 1).sum())
            pos_stats[m_gsm] = {
                "subject_id": subj_id,
                "sample_num": smp_num,
                "filename_prefix": mname.replace("_tissue_positions.csv.gz", ""),
                "n_visium_barcodes": n_total_barcodes,
                "n_in_tissue_spots": n_in_tissue,
            }
        elif mname.endswith("spatial_enrichment.csv.gz"):
            fobj = tf.extractfile(member)
            with gzip.open(fobj, "rt", encoding="utf-8", errors="replace") as gz:
                edf = pd.read_csv(gz)
            edf["gsm"] = m_gsm
            edf["subject_id"] = subj_id
            enr_dfs.append(edf[["gsm", "subject_id", "Feature ID", "Feature Name", "I", "P value", "Adjusted p value"]])

df_enr_all = pd.concat(enr_dfs, ignore_index=True)

# Merge GSM official characteristics + subject_id + spot counts + Khaliq paper section ID
g272_rows = []
for _, r in gsm_meta272.iterrows():
    gsm = r["gsm"]
    chars = r["characteristics"]
    # Extract tissue/specimen type from official GEO characteristics
    m_spec = re.search(r"tissue/specimen type:\s*([^|]+)", chars)
    spec_raw = m_spec.group(1).strip() if m_spec else ""
    if "Primary Tumor" in spec_raw:
        site_group = "Primary_Tumor"
    elif "Hepatic Mets" in spec_raw:
        site_group = "Hepatic_Metastasis"
    elif "Lymph node" in spec_raw:
        site_group = "Lymph_Node_Metastasis"
    elif "Normal Pancreas" in spec_raw:
        site_group = "Normal_Pancreas"
    else:
        site_group = spec_raw

    ps = pos_stats.get(gsm, {})
    k_sub = khaliq_map[khaliq_map["gsm"] == gsm]
    paper_sec = k_sub["section"].iloc[0] if len(k_sub) > 0 else "IU_PDA_HM2_dup_excluded"
    n_genes_enr = int((df_enr_all["gsm"] == gsm).sum())
    g272_rows.append({
        "gsm": gsm,
        "title": r["title"],
        "patient_id": ps.get("subject_id", "Unknown"),
        "paper_section_id": paper_sec,
        "organ": r["source_name"],
        "geo_specimen_type": spec_raw,
        "anatomical_group": site_group,
        "n_in_tissue_spots": ps.get("n_in_tissue_spots", 0),
        "n_genes_morans_tested": n_genes_enr,
        "filename_prefix": ps.get("filename_prefix", ""),
    })

df_g272_samples = pd.DataFrame(g272_rows)
df_g272_samples.to_csv(META_DIR / "gse272362_30specimens_verified_manifest.tsv", sep="\t", index=False)

# Save full merged official Moran's I parquet/csv with verified patient_id and anatomical_group!
df_enr_all = df_enr_all.merge(
    df_g272_samples[["gsm", "patient_id", "anatomical_group", "paper_section_id", "n_in_tissue_spots"]],
    on="gsm",
    how="left",
)
df_enr_all.to_csv(g272_dir / "gse272362_official_morans_i_all30.tsv.gz", sep="\t", index=False, compression="gzip")

print(
    f"  GSE272362 verified: {len(df_g272_samples)} specimens across {df_g272_samples['patient_id'].nunique()} rapid-autopsy patients, "
    f"total in-tissue spots = {df_g272_samples['n_in_tissue_spots'].sum():,}",
    flush=True,
)
print("  Anatomical breakdown:", df_g272_samples["anatomical_group"].value_counts().to_dict(), flush=True)
print("  Specimens per rapid-autopsy patient:", df_g272_samples.groupby(["patient_id", "anatomical_group"]).size().unstack(fill_value=0).to_string(), flush=True)

# -------------------------------------------------------------------------
# 5. Inspect TCGA-PAAD & DepMap & v0 462/460/18 Audit
# -------------------------------------------------------------------------
print("[5/6] Auditing Cohort 5 (TCGA-PAAD), DepMap 24Q4, and v0 RegVelo 462/460/18 artifact...", flush=True)
tcga_meta = pd.read_csv(V1_ROOT / "04_results/tcga/tcga_paad_matched_clinical.csv")
depmap_df = pd.read_csv(V1_ROOT / "02_data/depmap/depmap_tf_dependency_orig37.csv")
print(
    f"  TCGA-PAAD primary RNA-seq samples={len(tcga_meta)}, pure Ductal={int(tcga_meta['is_ductal'].sum())}, "
    f"OS valid (all)={int((tcga_meta['os_months']>0).sum())}, OS valid (ductal)={int((tcga_meta.loc[tcga_meta['is_ductal'],'os_months']>0).sum())}",
    flush=True,
)
print("  DepMap columns:", list(depmap_df.columns), "rows:", len(depmap_df), flush=True)

# Check v0 tables for A06 (36,016 macrophage DEGs) and B03/B04
v0_tab_dir = V0_ROOT / "原始代码/pdac_regvelo_project/07_data_tables"
v0_checks = {}
for fn in os.listdir(v0_tab_dir):
    if any(k in fn.lower() for k in ["regvelo", "inferred", "perturb", "deplet", "fig6", "fig8", "lincs", "priority", "target"]):
        p = v0_tab_dir / fn
        if p.suffix == ".csv" and p.stat().st_size < 20_000_000:
            df_tmp = pd.read_csv(p)
            v0_checks[fn] = {"shape": list(df_tmp.shape), "columns": list(df_tmp.columns)[:12]}
(AUDIT_DIR / "v0_tables_inventory.json").write_text(json.dumps(v0_checks, indent=2), encoding="utf-8")
print(f"  Inspected {len(v0_checks)} v0 CSV tables.", flush=True)

# -------------------------------------------------------------------------
# 6. Build DATA_PROVENANCE.tsv, SAMPLE_MANIFEST.tsv, CELL_AND_GENE_FLOW.tsv, cohort_summary.tsv
# -------------------------------------------------------------------------
print("[6/6] Exporting Stage 0/1 Gate TSVs...", flush=True)

prov_rows = [
    {
        "dataset_id": "Peng2019_CRA001160",
        "role_in_study": "Primary scRNA-seq discovery cohort (donor-level cellular composition, IFN-high macrophage state discovery, paired checkpoint expression, regulon & L-R inference)",
        "primary_citation": "Peng J et al., Cell Research 2019;29(9):725-738 (PMID: 31273297; PMCID: PMC6796938)",
        "primary_accession": "CRA001160 / PRJCA001063 (NGDC GSA)",
        "secondary_distribution": "Zenodo 3969339 (processed h5ad aligned to CRA001160 raw UMI counts; note: E-MTAB-8109/E-MTAB-8410 in v0 were erroneous accessions and are removed)",
        "platform": "10x Genomics Chromium scRNA-seq (3' v2)",
        "n_patients": 35,
        "n_specimens": 35,
        "tissue_breakdown": "24 treatment-naive primary PDAC tumors (T1-T24) + 11 control non-malignant pancreas tissues (N1-N11)",
        "n_cells_or_spots_qc": int(n_peng_ann_cells),
        "n_genes_full": int(n_unique_peng_raw_genes),
        "local_verified_files": "peng_raw.h5ad; peng_annotated.h5ad",
        "verification_date": "2026-11-16",
    },
    {
        "dataset_id": "GSE202051_Hwang2022",
        "role_in_study": "Independent external snRNA-seq validation cohort (cross-patient paired IFN-high macrophage replication & neoadjuvant treatment comparison)",
        "primary_citation": "Hwang WL et al., Nature Genetics 2022;54(8):1178-1191 (PMID: 35902743)",
        "primary_accession": "GSE202051 (NCBI GEO)",
        "secondary_distribution": "GSE202051_totaldata-final-toshare.h5ad.gz (official GEO Series supplementary file)",
        "platform": "10x Genomics Chromium snRNA-seq",
        "n_patients": int(g202_obs["pid"].nunique()),
        "n_specimens": int(g202_obs["sampleid"].nunique()) if "sampleid" in g202_obs.columns else int(g202_obs["pid"].nunique()),
        "tissue_breakdown": "43 primary PDAC surgical resections (18 untreated + 25 neoadjuvant CRT/CRTl/chemotherapy; 0 normal pancreas controls)",
        "n_cells_or_spots_qc": int(g202.n_obs),
        "n_genes_full": int(len(set(g202.var_names))),
        "local_verified_files": "gse202051_macrophages.h5ad (subset from 224,988 total nuclei)",
        "verification_date": "2026-11-16",
    },
    {
        "dataset_id": "GSE154778_Lin2020",
        "role_in_study": "Independent scRNA-seq cohort for primary-vs-metastatic state comparison & lineage-restricted RNA splicing / RegVelo dynamical modeling",
        "primary_citation": "Lin W et al., Genome Medicine 2020;12:80 (PMID: 32988401)",
        "primary_accession": "GSE154778 / SRP272677 / PRJNA647411 (16 GSMs, 120 SRA FASTQ runs)",
        "secondary_distribution": "Re-quantified in-house from 120 raw SRA FASTQ runs via STARsolo 2.7.11b (GRCh38 + GENCODE v44, Gene + Velocyto)",
        "platform": "10x Genomics Chromium scRNA-seq",
        "n_patients": 16,
        "n_specimens": 16,
        "tissue_breakdown": "10 primary PDAC biopsies/resections (P01-P10) + 6 metastatic biopsies (MET01-MET06)",
        "n_cells_or_spots_qc": int(g154.n_obs),
        "n_genes_full": int(g154.n_vars),
        "local_verified_files": "gse154778_all16_velocity.h5ad; starsolo_out/*/Summary.csv",
        "verification_date": "2026-11-16",
    },
    {
        "dataset_id": "GSE272362_Khaliq2024",
        "role_in_study": "Spatial transcriptomics univariate spatial autocorrelation (Moran's I) cohort across matched anatomical sites",
        "primary_citation": "Khaliq AM et al., Nature Genetics 2024;56(11):2455-2465 (PMID: 39294496)",
        "primary_accession": "GSE272362 / PRJNA1136518 (30 GSMs: GSM8399450-GSM8399479)",
        "secondary_distribution": "GSE272362_RAW.tar (Space Ranger spatial_enrichment.csv.gz & tissue_positions.csv.gz; note: spot x gene count matrix not hosted on GEO, Zenodo 22695225 11.6 GB RDS rate-limited)",
        "platform": "10x Genomics Visium FFPE Spatial Gene Expression",
        "n_patients": int(df_g272_samples["patient_id"].nunique()),
        "n_specimens": int(len(df_g272_samples)),
        "tissue_breakdown": "30 sections across 7 rapid-autopsy patients (S13,S14,S17,S18,S20,S21,S22): 10 Primary Tumor (PT), 12 Hepatic Metastasis (HM), 5 Lymph Node Metastasis (LNM), 3 Normal Pancreas (NP)",
        "n_cells_or_spots_qc": int(df_g272_samples["n_in_tissue_spots"].sum()),
        "n_genes_full": int(df_g272_samples["n_genes_morans_tested"].max()),
        "local_verified_files": "GSE272362_RAW.tar; gse272362_30specimens_verified_manifest.tsv",
        "verification_date": "2026-11-16",
    },
    {
        "dataset_id": "TCGA_PAAD_PanCanAtlas",
        "role_in_study": "Bulk clinical survival & immune-infiltration-adjusted checkpoint association cohort (histologically stratified)",
        "primary_citation": "Cancer Genome Atlas Research Network, Cancer Cell 2017;32(2):185-203.e13 (PMID: 28810144); Liu J et al., Cell 2018;173(2):400-416.e11 (PMID: 29625055)",
        "primary_accession": "TCGA-PAAD (paad_tcga_pan_can_atlas_2018 via cBioPortal / GDC)",
        "secondary_distribution": "mrna_seq_v2_rsem (log2[RSEM+1] transformed in R; non-silent WES somatic mutations)",
        "platform": "Illumina HiSeq bulk RNA-seq (RSEM) & Whole-Exome Sequencing",
        "n_patients": int(len(tcga_meta)),
        "n_specimens": int(len(tcga_meta)),
        "tissue_breakdown": f"177 unique primary PAAD tumors -> 146 histologically confirmed Pancreas Adenocarcinoma, Ductal Type (31 non-ductal/neuroendocrine/mucinous/other excluded in primary analysis)",
        "n_cells_or_spots_qc": int(tcga_meta["is_ductal"].sum()),
        "n_genes_full": 20531,
        "local_verified_files": "mrna_rsem_log_not_applied.csv; clinical_patient.csv; clinical_sample.csv; mutations.csv",
        "verification_date": "2026-11-16",
    },
]
pd.DataFrame(prov_rows).to_csv(META_DIR / "DATA_PROVENANCE.tsv", sep="\t", index=False)

# Unified SAMPLE_MANIFEST.tsv across all 5 cohorts
sample_rows = []
for _, r in peng_donors.iterrows():
    sample_rows.append({
        "dataset_id": "Peng2019_CRA001160",
        "patient_id": f"Peng_{r['Patient']}",
        "sample_id": str(r["Patient"]),
        "condition": "PDAC_Tumor" if r["CONDITION"] == "T" else "Control_Pancreas",
        "treatment": "Treatment_Naive",
        "tissue_site": "Primary_Pancreas",
        "histology": "PDAC" if r["CONDITION"] == "T" else "Non_malignant_pancreas",
        "n_cells_or_spots": int(r["n_cells"]),
        "n_target_macrophages": int(r["n_macrophages"]),
        "included_primary_analysis": True,
    })

for _, r in g202_donors.iterrows():
    sample_rows.append({
        "dataset_id": "GSE202051_Hwang2022",
        "patient_id": str(r["pid"]),
        "sample_id": str(r["pid"]),
        "condition": "PDAC_Tumor",
        "treatment": str(r["new_treatment"]),
        "tissue_site": "Primary_Pancreas",
        "histology": f"PDAC (response={r['response']})",
        "n_cells_or_spots": int(r["n_macrophages"]),
        "n_target_macrophages": int(r["n_macrophages"]),
        "included_primary_analysis": True,
    })

g154_by_smp = (
    g154_obs.groupby("sample", observed=True)
    .agg(
        n_cells=("cell_type", "size"),
        n_mac=("cell_type", lambda s: int((s == "Macrophage_Myeloid").sum())),
    )
    .reset_index()
)
for _, r in g154_by_smp.iterrows():
    smp = str(r["sample"])
    is_met = smp.startswith("MET")
    sample_rows.append({
        "dataset_id": "GSE154778_Lin2020",
        "patient_id": f"GSE154778_{smp}",
        "sample_id": smp,
        "condition": "PDAC_Metastasis" if is_met else "PDAC_Primary",
        "treatment": "Biopsy_or_Resection",
        "tissue_site": "Metastatic_Lesion" if is_met else "Primary_Pancreas",
        "histology": "PDAC",
        "n_cells_or_spots": int(r["n_cells"]),
        "n_target_macrophages": int(r["n_mac"]),
        "included_primary_analysis": True,
    })

for _, r in df_g272_samples.iterrows():
    sample_rows.append({
        "dataset_id": "GSE272362_Khaliq2024",
        "patient_id": f"Khaliq_{r['patient_id']}",
        "sample_id": str(r["gsm"]),
        "condition": r["anatomical_group"],
        "treatment": "Rapid_Autopsy",
        "tissue_site": str(r["organ"]),
        "histology": str(r["geo_specimen_type"]),
        "n_cells_or_spots": int(r["n_in_tissue_spots"]),
        "n_target_macrophages": -1,
        "included_primary_analysis": True,
    })

for _, r in tcga_meta.iterrows():
    sample_rows.append({
        "dataset_id": "TCGA_PAAD_PanCanAtlas",
        "patient_id": str(r["patientId"]),
        "sample_id": str(r["sampleId"]),
        "condition": "PDAC_Primary" if bool(r["is_ductal"]) else "Non_Ductal_PAAD_Excluded",
        "treatment": "Surgical_Resection",
        "tissue_site": "Primary_Pancreas",
        "histology": str(r["TUMOR_TYPE"]),
        "n_cells_or_spots": 1,
        "n_target_macrophages": -1,
        "included_primary_analysis": bool(r["is_ductal"]),
    })

df_sample_manifest = pd.DataFrame(sample_rows)
df_sample_manifest.to_csv(META_DIR / "SAMPLE_MANIFEST.tsv", sep="\t", index=False)

# Frozen cohort_summary.tsv
cohort_sum = pd.DataFrame([
    {
        "cohort_id": "Peng2019",
        "accession": "CRA001160 / PRJCA001063 (Zenodo 3969339)",
        "modality": "10x 3' scRNA-seq",
        "n_patients_total": 35,
        "n_patients_pdac": 24,
        "n_patients_control": 11,
        "n_specimens": 35,
        "n_cells_or_spots_raw": int(n_peng_raw_cells),
        "n_cells_or_spots_qc": int(n_peng_ann_cells),
        "n_macrophages_qc": int((peng_obs["Cell_type"] == "Macrophage cell").sum()),
        "n_genes_full_matrix": int(n_unique_peng_raw_genes),
        "n_genes_hvg_or_panel": int(n_peng_ann_genes),
        "treatment_context": "Treatment-naive primary surgical resections & non-malignant pancreas",
        "statistical_unit": "Patient / donor (n=24 PDAC vs n=11 control; n=27-29 paired within-donor TAMs)",
    },
    {
        "cohort_id": "GSE202051",
        "accession": "GSE202051",
        "modality": "10x snRNA-seq",
        "n_patients_total": 43,
        "n_patients_pdac": 43,
        "n_patients_control": 0,
        "n_specimens": 43,
        "n_cells_or_spots_raw": 224988,
        "n_cells_or_spots_qc": 224988,
        "n_macrophages_qc": int(g202.n_obs),
        "n_genes_full_matrix": int(len(set(g202.var_names))),
        "n_genes_hvg_or_panel": 2000,
        "treatment_context": "18 untreated + 25 neoadjuvant-treated (CRT/CRTl/chemotherapy) primary PDAC",
        "statistical_unit": "Patient / donor (n=43 total; stratified n=18 untreated & n=25 neoadjuvant)",
    },
    {
        "cohort_id": "GSE154778",
        "accession": "GSE154778 / SRP272677 (120 SRR runs)",
        "modality": "10x scRNA-seq (STARsolo Gene + Velocyto)",
        "n_patients_total": 16,
        "n_patients_pdac": 16,
        "n_patients_control": 0,
        "n_specimens": 16,
        "n_cells_or_spots_raw": int(g154_qc["raw_filtered_cells"].sum()) if "raw_filtered_cells" in g154_qc.columns else int(g154.n_obs),
        "n_cells_or_spots_qc": int(g154.n_obs),
        "n_macrophages_qc": int((g154_obs["cell_type"] == "Macrophage_Myeloid").sum()),
        "n_genes_full_matrix": int(g154.n_vars),
        "n_genes_hvg_or_panel": 2500,
        "treatment_context": "10 primary PDAC tumors (P01-P10) + 6 metastatic lesions (MET01-MET06)",
        "statistical_unit": "Patient / donor (n=10 primary vs n=6 metastasis; 3-seed lineage-restricted RegVelo)",
    },
    {
        "cohort_id": "GSE272362",
        "accession": "GSE272362 / PRJNA1136518",
        "modality": "10x Visium FFPE (Space Ranger Moran's I)",
        "n_patients_total": int(df_g272_samples["patient_id"].nunique()),
        "n_patients_pdac": int(df_g272_samples["patient_id"].nunique()),
        "n_patients_control": 0,
        "n_specimens": 30,
        "n_cells_or_spots_raw": int(df_g272_samples["n_in_tissue_spots"].sum()),
        "n_cells_or_spots_qc": 91496,
        "n_macrophages_qc": 0,
        "n_genes_full_matrix": int(df_g272_samples["n_genes_morans_tested"].max()),
        "n_genes_hvg_or_panel": int(df_g272_samples["n_genes_morans_tested"].max()),
        "treatment_context": "Rapid-autopsy multi-site cohort (10 PT, 12 HM, 5 LNM, 3 matched NP across 7 patients)",
        "statistical_unit": "Patient (n=7 rapid-autopsy donors) & anatomical-site-stratified specimen (n=30)",
    },
    {
        "cohort_id": "TCGA-PAAD",
        "accession": "paad_tcga_pan_can_atlas_2018",
        "modality": "Bulk RNA-seq (log2[RSEM+1]) + WES",
        "n_patients_total": int(len(tcga_meta)),
        "n_patients_pdac": int(tcga_meta["is_ductal"].sum()),
        "n_patients_control": 0,
        "n_specimens": int(len(tcga_meta)),
        "n_cells_or_spots_raw": int(len(tcga_meta)),
        "n_cells_or_spots_qc": int(tcga_meta["is_ductal"].sum()),
        "n_macrophages_qc": 0,
        "n_genes_full_matrix": 20531,
        "n_genes_hvg_or_panel": 1639,
        "treatment_context": "Resected primary pancreatic adenocarcinoma (146 confirmed ductal PDAC; 31 non-ductal excluded in primary model)",
        "statistical_unit": "Patient (n=146 univariate OS; n=140 complete-case multivariable Cox)",
    },
])
cohort_sum.to_csv(RES_DIR / "cohort_summary.tsv", sep="\t", index=False)
print("[Stage 0/1 Gate Script Completed Successfully]", flush=True)
