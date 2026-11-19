#!/usr/bin/env python3
"""
Stage 3 (Gate G3): Non-Circular IFN-High Macrophage State Rebuild, Quantile Q1 vs Q2 Separation,
Ambient RNA Soup Cancellation Proof, and Multi-Cohort Patient-Level Validation (Peng2019, GSE202051, GSE154778).
"""
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

if hasattr(sp, "csr_array") and not hasattr(sp.csr_array, "getnnz"):
    sp.csr_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)
if hasattr(sp, "csc_array") and not hasattr(sp.csc_array, "getnnz"):
    sp.csc_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)

import anndata as ad

if os.name == "nt":
    OUT_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
    V1_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo")
    V0_ROOT = Path(r"D:\09-04\09-04")
else:
    OUT_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116")
    V1_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo")
    V0_ROOT = Path("/mnt/d/09-04/09-04")

ST3_DIR = OUT_ROOT / "05_results/stage3_macrophage"
PB_DIR = ST3_DIR / "pseudobulk_inputs"
PB_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------------
# Signature definitions (locked in ANALYSIS_PLAN_LOCKED.md)
# -------------------------------------------------------------------------
IFN15_LEGACY = [
    "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11",
    "MX1", "MX2", "IFIT1", "IFIT3", "GBP1", "GBP5",
    "SOCS1", "SOCS3", "JAK1", "JAK2",
]
# Primary non-circular IFN-14 (excludes STAT1 and all checkpoints/HLA genes)
IFN14_PRIMARY = [g for g in IFN15_LEGACY if g != "STAT1"]
# Ultra-strict non-circular IFN-10 antiviral ISG core (excludes STAT1, IRF1, AND CXCL9/10/11)
IFN10_STRICT = ["MX1", "MX2", "IFIT1", "IFIT3", "GBP1", "GBP5", "SOCS1", "SOCS3", "JAK1", "JAK2"]

PRESPECIFIED_PANEL = [
    # Transcription factors (independent of IFN-10; STAT1 independent of IFN-14)
    ("TF_Regulator", "STAT1"), ("TF_Regulator", "IRF1"), ("TF_Regulator", "IRF7"), ("TF_Regulator", "IRF8"),
    # Immune checkpoint & immunoregulatory ligands (100% excluded from ALL IFN signatures)
    ("Checkpoint_Immunoreg", "LGALS9"), ("Checkpoint_Immunoreg", "CD47"), ("Checkpoint_Immunoreg", "CD80"),
    ("Checkpoint_Immunoreg", "CD86"), ("Checkpoint_Immunoreg", "NECTIN2"), ("Checkpoint_Immunoreg", "PVR"),
    ("Checkpoint_Immunoreg", "CD274"), ("Checkpoint_Immunoreg", "PDCD1LG2"), ("Checkpoint_Immunoreg", "VSIR"),
    ("Checkpoint_Immunoreg", "CD276"), ("Checkpoint_Immunoreg", "TGFB1"),
    # Intracellular metabolic enzyme (explicitly NOT a surface ligand)
    ("Metabolic_Enzyme", "IDO1"),
    # MHC-I / MHC-II antigen presentation machinery (100% excluded from ALL IFN signatures)
    ("Antigen_Presentation", "HLA-A"), ("Antigen_Presentation", "HLA-B"), ("Antigen_Presentation", "HLA-C"),
    ("Antigen_Presentation", "B2M"), ("Antigen_Presentation", "TAP1"), ("Antigen_Presentation", "TAP2"),
    ("Antigen_Presentation", "HLA-DRA"), ("Antigen_Presentation", "HLA-DRB1"), ("Antigen_Presentation", "CD74"),
    # Chemokines (independent of IFN-10)
    ("Chemokine", "CXCL9"), ("Chemokine", "CXCL10"), ("Chemokine", "CXCL11"), ("Chemokine", "CCL2"), ("Chemokine", "CCL5"),
    # Alternative TAM lineage markers & ambient-RNA soup negative controls
    ("Resident_Hypoxic_TAM", "C1QA"), ("Resident_Hypoxic_TAM", "MRC1"), ("Resident_Hypoxic_TAM", "CD163"), ("Resident_Hypoxic_TAM", "SPP1"),
    ("Ambient_Soup_Control", "KRT19"), ("Ambient_Soup_Control", "COL1A1"), ("Ambient_Soup_Control", "COL3A1"), ("Ambient_Soup_Control", "PRSS1"),
]

def normalize_cp10k_log1p(counts_csr):
    lib = np.asarray(counts_csr.sum(axis=1)).ravel().astype(np.float64)
    lib_safe = np.where(lib > 0, lib, 1.0)
    inv = sp.diags(1e4 / lib_safe)
    norm = inv @ counts_csr
    norm_dense = norm.toarray().astype(np.float32)
    return np.log1p(norm_dense)

def score_genes_zmean(X_log, gene_to_idx, genes):
    avail = [g for g in genes if g in gene_to_idx]
    sub = X_log[:, [gene_to_idx[g] for g in avail]]
    mu = sub.mean(axis=0)
    sd = sub.std(axis=0) + 1e-8
    return ((sub - mu) / sd).mean(axis=1), sub.mean(axis=1), avail

def assign_q1_global(vals):
    q1, q2 = np.percentile(vals, [33.3333, 66.6667])
    out = np.array(["mid"] * len(vals), dtype=object)
    out[vals <= q1] = "low"
    out[vals >= q2] = "high"
    return out, float(q1), float(q2)

def assign_q2_within_donor(vals, donor_vec, min_donor_cells=15):
    out = np.array(["unassigned"] * len(vals), dtype=object)
    df_tmp = pd.DataFrame({"v": vals, "d": donor_vec})
    for d, idx in df_tmp.groupby("d").groups.items():
        if len(idx) < min_donor_cells:
            continue
        sub_v = vals[idx]
        q1, q2 = np.percentile(sub_v, [33.3333, 66.6667])
        if q2 <= q1:
            # Handle ties cleanly: median split if tertiles tie
            med = np.median(sub_v)
            lbl = np.array(["mid"] * len(idx), dtype=object)
            lbl[sub_v < med] = "low"
            lbl[sub_v > med] = "high"
        else:
            lbl = np.array(["mid"] * len(idx), dtype=object)
            lbl[sub_v <= q1] = "low"
            lbl[sub_v >= q2] = "high"
        out[idx] = lbl
    return out

def bootstrap_hl_diff(x, y, n_boot=2000, seed=42):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Exact Hodges-Lehmann estimator = median of all pairwise differences (x_i - y_j)
    diffs = np.subtract.outer(x, y).ravel()
    hl = float(np.median(diffs))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        xb = rng.choice(x, size=len(x), replace=True)
        yb = rng.choice(y, size=len(y), replace=True)
        boots.append(np.median(np.subtract.outer(xb, yb).ravel()))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return hl, float(lo), float(hi)

# =========================================================================
# 1. COHORT 1: PENG ET AL. 2019 (DISCOVERY)
# =========================================================================
print("[1/4] Processing Peng 2019 Discovery Macrophages (5,361 TAMs across 35 donors)...", flush=True)
peng_raw = ad.read_h5ad(V1_ROOT / "02_data/peng2019/peng_raw.h5ad")
peng_ann = ad.read_h5ad(V1_ROOT / "02_data/peng2019/peng_annotated.h5ad")
common_c = peng_ann.obs_names.intersection(peng_raw.obs_names)
peng_all = peng_raw[common_c].copy()
for col in peng_ann.obs.columns:
    peng_all.obs[col] = peng_ann.obs.loc[common_c, col].values

total_cells_per_donor = peng_all.obs.groupby("Patient", observed=True).size().to_dict()
mac_p = peng_all[peng_all.obs["Cell_type"] == "Macrophage cell"].copy()
counts_p = mac_p.X.tocsr() if sp.issparse(mac_p.X) else sp.csr_matrix(mac_p.X)
Xlog_p = normalize_cp10k_log1p(counts_p)
g2i_p = {g: i for i, g in enumerate(mac_p.var_names)}

z14_p, rawmean14_p, av14_p = score_genes_zmean(Xlog_p, g2i_p, IFN14_PRIMARY)
z10_p, rawmean10_p, av10_p = score_genes_zmean(Xlog_p, g2i_p, IFN10_STRICT)
z15_p, rawmean15_p, av15_p = score_genes_zmean(Xlog_p, g2i_p, IFN15_LEGACY)
stat1_p = Xlog_p[:, g2i_p["STAT1"]]

q1_14_p, cut_lo14_p, cut_hi14_p = assign_q1_global(z14_p)
q1_10_p, _, _ = assign_q1_global(z10_p)
q1_15_p, _, _ = assign_q1_global(z15_p)
# Fixed transferable unstandardized log1p mean cutoff from Peng discovery (66.67th percentile of rawmean14_p)
fixed_rawmean14_hi_cut = float(np.percentile(rawmean14_p, 66.6667))
fixed_rawmean14_lo_cut = float(np.percentile(rawmean14_p, 33.3333))

donors_p_vec = mac_p.obs["Patient"].astype(str).values
cond_p_vec = mac_p.obs["CONDITION"].astype(str).values
q2_14_p = assign_q2_within_donor(z14_p, donors_p_vec, min_donor_cells=15)
q2_10_p = assign_q2_within_donor(z10_p, donors_p_vec, min_donor_cells=15)

# Merge UMAP coordinates from v1 for consistent visualization
v1_umap_p = pd.read_csv(V1_ROOT / "04_results/fig6_macrophage/peng_mac_umap_and_scores.csv").set_index("barcode")
mac_p_cells_df = pd.DataFrame({
    "barcode": mac_p.obs_names,
    "patient": donors_p_vec,
    "condition": np.where(cond_p_vec == "T", "PDAC", "Control"),
    "UMAP1": v1_umap_p.reindex(mac_p.obs_names)["UMAP1"].values,
    "UMAP2": v1_umap_p.reindex(mac_p.obs_names)["UMAP2"].values,
    "stat1_log1p": stat1_p,
    "ifn14_zscore": z14_p,
    "ifn14_rawmean_log1p": rawmean14_p,
    "ifn10_zscore": z10_p,
    "ifn15_legacy_zscore": z15_p,
    "Q1_global_IFN14": q1_14_p,
    "Q2_within_donor_IFN14": q2_14_p,
    "Q1_global_IFN10": q1_10_p,
    "Q2_within_donor_IFN10": q2_10_p,
})
mac_p_cells_df.to_csv(ST3_DIR / "peng2019_macrophage_cell_states.tsv", sep="\t", index=False)

# Between-donor abundance summary (Rule Q1 ONLY)
peng_donor_rows = []
for pat, grp in mac_p_cells_df.groupby("patient"):
    n_mac = len(grp)
    n_tot = int(total_cells_per_donor.get(pat, n_mac))
    n_hi_q1 = int((grp["Q1_global_IFN14"] == "high").sum())
    n_hi_10 = int((grp["Q1_global_IFN10"] == "high").sum())
    peng_donor_rows.append({
        "cohort": "Peng2019",
        "patient": pat,
        "condition": grp["condition"].iloc[0],
        "n_total_cells": n_tot,
        "n_macrophages": n_mac,
        "mac_frac_of_all_cells": n_mac / n_tot,
        "ifn14_high_frac_of_mac_Q1": n_hi_q1 / n_mac,
        "ifn14_high_frac_of_all_cells_Q1": n_hi_q1 / n_tot,
        "ifn10_high_frac_of_mac_Q1": n_hi_10 / n_mac,
        "donor_mean_ifn14_z": float(grp["ifn14_zscore"].mean()),
        "donor_mean_ifn10_z": float(grp["ifn10_zscore"].mean()),
        "donor_mean_stat1_log1p": float(grp["stat1_log1p"].mean()),
        "donor_pct_stat1_pos": float((grp["stat1_log1p"] > 0).mean() * 100.0),
    })
# Check if any donor in Peng had 0 macrophages!
all_peng_pats = peng_all.obs[["Patient", "CONDITION"]].drop_duplicates()
present_pats = set(r["patient"] for r in peng_donor_rows)
for _, r in all_peng_pats.iterrows():
    if str(r["Patient"]) not in present_pats:
        n_tot = int(total_cells_per_donor[r["Patient"]])
        peng_donor_rows.append({
            "cohort": "Peng2019", "patient": str(r["Patient"]),
            "condition": "PDAC" if r["CONDITION"] == "T" else "Control",
            "n_total_cells": n_tot, "n_macrophages": 0, "mac_frac_of_all_cells": 0.0,
            "ifn14_high_frac_of_mac_Q1": np.nan, "ifn14_high_frac_of_all_cells_Q1": 0.0,
            "ifn10_high_frac_of_mac_Q1": np.nan, "donor_mean_ifn14_z": np.nan,
            "donor_mean_ifn10_z": np.nan, "donor_mean_stat1_log1p": np.nan, "donor_pct_stat1_pos": np.nan,
        })
df_peng_donors = pd.DataFrame(peng_donor_rows).sort_values(["condition", "patient"])
df_peng_donors.to_csv(ST3_DIR / "peng2019_donor_abundance_Q1.tsv", sep="\t", index=False)

# Statistical tests on between-donor abundance (all donors with >=1 mac & stability >=10 macs)
abundance_tests = []
for min_m in [1, 10]:
    sub = df_peng_donors[df_peng_donors["n_macrophages"] >= min_m]
    tum = sub[sub["condition"] == "PDAC"]
    ctl = sub[sub["condition"] == "Control"]
    for metric in [
        "ifn14_high_frac_of_mac_Q1",
        "ifn14_high_frac_of_all_cells_Q1",
        "ifn10_high_frac_of_mac_Q1",
        "donor_mean_ifn14_z",
        "donor_mean_stat1_log1p",
    ]:
        x = tum[metric].dropna().values
        y = ctl[metric].dropna().values
        u_stat, p_val = stats.mannwhitneyu(x, y, alternative="two-sided")
        hl, ci_lo, ci_hi = bootstrap_hl_diff(x, y)
        abundance_tests.append({
            "cohort": "Peng2019",
            "contrast": "PDAC_vs_Control_Pancreas",
            "min_macrophages_per_donor": min_m,
            "metric": metric,
            "n_group1": len(x),
            "n_group2": len(y),
            "median_group1": float(np.median(x)),
            "iqr_group1": float(np.percentile(x, 75) - np.percentile(x, 25)),
            "median_group2": float(np.median(y)),
            "iqr_group2": float(np.percentile(y, 75) - np.percentile(y, 25)),
            "hodges_lehmann_diff": hl,
            "ci95_low": ci_lo,
            "ci95_high": ci_hi,
            "mwu_p": float(p_val),
        })

# Export Peng within-donor paired pseudobulk count matrices for edgeR QL:
# Primary: Q2_within_donor_IFN14 (min_cells=5), plus targeted sensitivities: Q2 min_cells=3, Q2 min_cells=10, Q1_global_IFN14 min_cells=5, Q2_within_donor_IFN10 min_cells=5
def export_paired_pseudobulk(counts_csr, gene_names, donor_vec, grp_vec, extra_meta_dict, tag, min_cells=5):
    cols = {}
    mrows = []
    for d in sorted(set(donor_vec)):
        m_hi = (donor_vec == d) & (grp_vec == "high")
        m_lo = (donor_vec == d) & (grp_vec == "low")
        if m_hi.sum() >= min_cells and m_lo.sum() >= min_cells:
            for lbl, mask in [("high", m_hi), ("low", m_lo)]:
                sid = f"{d}__{lbl}"
                cols[sid] = np.asarray(counts_csr[mask].sum(axis=0)).ravel().astype(np.int64)
                row = {"sample_id": sid, "patient": d, "group": lbl, "n_cells": int(mask.sum())}
                if d in extra_meta_dict:
                    row.update(extra_meta_dict[d])
                mrows.append(row)
    pd.DataFrame(cols, index=gene_names).to_csv(PB_DIR / f"counts_{tag}.csv")
    mdf = pd.DataFrame(mrows)
    mdf.to_csv(PB_DIR / f"meta_{tag}.csv", index=False)
    print(f"  Exported {tag}: {len(mdf)//2} paired donors (min_cells={min_cells})", flush=True)
    return len(mdf) // 2

peng_extra = {r["patient"]: {"condition": r["condition"]} for _, r in df_peng_donors.iterrows()}
export_paired_pseudobulk(counts_p, mac_p.var_names, donors_p_vec, q2_14_p, peng_extra, "Peng_Q2_IFN14_min5", min_cells=5)
export_paired_pseudobulk(counts_p, mac_p.var_names, donors_p_vec, q2_14_p, peng_extra, "Peng_Q2_IFN14_min3", min_cells=3)
export_paired_pseudobulk(counts_p, mac_p.var_names, donors_p_vec, q2_14_p, peng_extra, "Peng_Q2_IFN14_min10", min_cells=10)
export_paired_pseudobulk(counts_p, mac_p.var_names, donors_p_vec, q1_14_p, peng_extra, "Peng_Q1_IFN14_min5", min_cells=5)
export_paired_pseudobulk(counts_p, mac_p.var_names, donors_p_vec, q2_10_p, peng_extra, "Peng_Q2_IFN10_min5", min_cells=5)
# Also PDAC-tumor-only donors in Peng (excluding normal pancreas donors) to verify effect is purely intratumoral!
tum_mask_p = cond_p_vec == "T"
export_paired_pseudobulk(
    counts_p[tum_mask_p], mac_p.var_names, donors_p_vec[tum_mask_p], q2_14_p[tum_mask_p],
    peng_extra, "Peng_PDAConly_Q2_IFN14_min5", min_cells=5
)

# =========================================================================
# 2. COHORT 3: GSE202051 (EXTERNAL snRNA-seq VALIDATION, 43 PATIENTS)
# =========================================================================
print("[2/4] Processing GSE202051 External Validation Macrophages (10,140 TAMs across 43 PDAC patients)...", flush=True)
g202 = ad.read_h5ad(V1_ROOT / "02_data/gse202051/gse202051_macrophages.h5ad")
counts_202 = g202.X.tocsr() if sp.issparse(g202.X) else sp.csr_matrix(g202.X)
Xlog_202 = normalize_cp10k_log1p(counts_202)
g2i_202 = {g: i for i, g in enumerate(g202.var_names)}

z14_202, rawmean14_202, _ = score_genes_zmean(Xlog_202, g2i_202, IFN14_PRIMARY)
z10_202, rawmean10_202, _ = score_genes_zmean(Xlog_202, g2i_202, IFN10_STRICT)
stat1_202 = Xlog_202[:, g2i_202["STAT1"]]

q1_14_202, _, _ = assign_q1_global(z14_202)
q1_10_202, _, _ = assign_q1_global(z10_202)
# Fixed discovery-calibrated threshold transfer (without cohort re-quantile!)
fixed_cal_202 = np.array(["mid"] * len(rawmean14_202), dtype=object)
fixed_cal_202[rawmean14_202 <= fixed_rawmean14_lo_cut] = "low"
fixed_cal_202[rawmean14_202 >= fixed_rawmean14_hi_cut] = "high"

pids_202 = g202.obs["pid"].astype(str).values
trt_raw_202 = g202.obs["new_treatment"].astype(str).values
trt_bin_202 = np.where(trt_raw_202 == "Untreated", "Untreated", "Neoadjuvant")
resp_202 = g202.obs["response"].astype(str).values

q2_14_202 = assign_q2_within_donor(z14_202, pids_202, min_donor_cells=15)
q2_10_202 = assign_q2_within_donor(z10_202, pids_202, min_donor_cells=15)

v1_umap_202 = pd.read_csv(V1_ROOT / "04_results/gse202051/gse202051_mac_umap_and_scores.csv").set_index("barcode")
mac_202_cells_df = pd.DataFrame({
    "barcode": g202.obs_names,
    "patient": pids_202,
    "treatment_binary": trt_bin_202,
    "treatment_arm": trt_raw_202,
    "pathologic_response": resp_202,
    "UMAP1": v1_umap_202.reindex(g202.obs_names)["UMAP1"].values,
    "UMAP2": v1_umap_202.reindex(g202.obs_names)["UMAP2"].values,
    "stat1_log1p": stat1_202,
    "ifn14_zscore": z14_202,
    "ifn14_rawmean_log1p": rawmean14_202,
    "ifn10_zscore": z10_202,
    "Q1_global_IFN14": q1_14_202,
    "Q1_fixed_PengCalibrated_IFN14": fixed_cal_202,
    "Q2_within_donor_IFN14": q2_14_202,
    "Q2_within_donor_IFN10": q2_10_202,
})
mac_202_cells_df.to_csv(ST3_DIR / "gse202051_macrophage_cell_states.tsv", sep="\t", index=False)

g202_donor_rows = []
for pid, grp in mac_202_cells_df.groupby("patient"):
    n_mac = len(grp)
    g202_donor_rows.append({
        "cohort": "GSE202051",
        "patient": pid,
        "treatment_binary": grp["treatment_binary"].iloc[0],
        "treatment_arm": grp["treatment_arm"].iloc[0],
        "pathologic_response": grp["pathologic_response"].iloc[0],
        "n_macrophages": n_mac,
        "ifn14_high_frac_of_mac_Q1": float((grp["Q1_global_IFN14"] == "high").mean()),
        "ifn14_high_frac_fixed_PengCal": float((grp["Q1_fixed_PengCalibrated_IFN14"] == "high").mean()),
        "ifn10_high_frac_of_mac_Q1": float((grp["Q1_global_IFN14"] == "high").mean()),
        "donor_mean_ifn14_z": float(grp["ifn14_zscore"].mean()),
        "donor_mean_stat1_log1p": float(grp["stat1_log1p"].mean()),
        "donor_pct_stat1_pos": float((grp["stat1_log1p"] > 0).mean() * 100.0),
    })
df_g202_donors = pd.DataFrame(g202_donor_rows).sort_values(["treatment_binary", "patient"])
df_g202_donors.to_csv(ST3_DIR / "gse202051_donor_abundance_and_treatment_Q1.tsv", sep="\t", index=False)

for min_m in [1, 10]:
    sub = df_g202_donors[df_g202_donors["n_macrophages"] >= min_m]
    neo = sub[sub["treatment_binary"] == "Neoadjuvant"]
    unt = sub[sub["treatment_binary"] == "Untreated"]
    for metric in ["ifn14_high_frac_of_mac_Q1", "ifn14_high_frac_fixed_PengCal", "donor_mean_ifn14_z", "donor_mean_stat1_log1p"]:
        x = neo[metric].dropna().values
        y = unt[metric].dropna().values
        _, p_val = stats.mannwhitneyu(x, y, alternative="two-sided")
        hl, ci_lo, ci_hi = bootstrap_hl_diff(x, y)
        abundance_tests.append({
            "cohort": "GSE202051",
            "contrast": "Neoadjuvant_vs_Untreated_PDAC",
            "min_macrophages_per_donor": min_m,
            "metric": metric,
            "n_group1": len(x),
            "n_group2": len(y),
            "median_group1": float(np.median(x)),
            "iqr_group1": float(np.percentile(x, 75) - np.percentile(x, 25)),
            "median_group2": float(np.median(y)),
            "iqr_group2": float(np.percentile(y, 75) - np.percentile(y, 25)),
            "hodges_lehmann_diff": hl,
            "ci95_low": ci_lo,
            "ci95_high": ci_hi,
            "mwu_p": float(p_val),
        })

g202_extra = {r["patient"]: {"treatment_binary": r["treatment_binary"], "treatment_arm": r["treatment_arm"]} for _, r in df_g202_donors.iterrows()}
export_paired_pseudobulk(counts_202, g202.var_names, pids_202, q2_14_202, g202_extra, "GSE202051_All43_Q2_IFN14_min5", min_cells=5)
unt_mask = trt_bin_202 == "Untreated"
neo_mask = trt_bin_202 == "Neoadjuvant"
export_paired_pseudobulk(counts_202[unt_mask], g202.var_names, pids_202[unt_mask], q2_14_202[unt_mask], g202_extra, "GSE202051_Untreated18_Q2_IFN14_min5", min_cells=5)
export_paired_pseudobulk(counts_202[neo_mask], g202.var_names, pids_202[neo_mask], q2_14_202[neo_mask], g202_extra, "GSE202051_Neoadjuvant25_Q2_IFN14_min5", min_cells=5)
export_paired_pseudobulk(counts_202, g202.var_names, pids_202, q1_14_202, g202_extra, "GSE202051_All43_Q1_IFN14_min5", min_cells=5)
export_paired_pseudobulk(counts_202, g202.var_names, pids_202, q2_10_202, g202_extra, "GSE202051_All43_Q2_IFN10_min5", min_cells=5)

# =========================================================================
# 3. COHORT 2: GSE154778 (16 PATIENTS: 10 PRIMARY + 6 METASTASIS)
# =========================================================================
print("[3/4] Processing GSE154778 Macrophage_Myeloid (3,318 cells across 16 PDAC patients)...", flush=True)
g154_full = ad.read_h5ad(V1_ROOT / "02_data/gse154778/gse154778_all16_velocity.h5ad")
mac_154 = g154_full[g154_full.obs["cell_type"] == "Macrophage_Myeloid"].copy()
counts_154 = mac_154.layers["counts"].tocsr() if sp.issparse(mac_154.layers["counts"]) else sp.csr_matrix(mac_154.layers["counts"])
Xlog_154 = normalize_cp10k_log1p(counts_154)
g2i_154 = {g: i for i, g in enumerate(mac_154.var_names)}

z14_154, rawmean14_154, _ = score_genes_zmean(Xlog_154, g2i_154, IFN14_PRIMARY)
z10_154, _, _ = score_genes_zmean(Xlog_154, g2i_154, IFN10_STRICT)
stat1_154 = Xlog_154[:, g2i_154["STAT1"]]
q1_14_154, _, _ = assign_q1_global(z14_154)
smp_154 = mac_154.obs["sample"].astype(str).values
site_154 = np.where(pd.Series(smp_154).str.startswith("MET"), "Metastasis", "Primary")
q2_14_154 = assign_q2_within_donor(z14_154, smp_154, min_donor_cells=15)

g154_donor_rows = []
tot_154 = g154_full.obs.groupby("sample", observed=True).size().to_dict()
for s, grp_idx in pd.Series(smp_154).groupby(smp_154).groups.items():
    n_m = len(grp_idx)
    n_t = int(tot_154.get(s, n_m))
    g154_donor_rows.append({
        "cohort": "GSE154778",
        "patient": s,
        "lesion_site": "Metastasis" if s.startswith("MET") else "Primary",
        "n_total_cells": n_t,
        "n_macrophages": n_m,
        "ifn14_high_frac_of_mac_Q1": float((q1_14_154[grp_idx] == "high").mean()),
        "ifn14_high_frac_of_all_cells_Q1": float((q1_14_154[grp_idx] == "high").sum() / n_t),
        "donor_mean_ifn14_z": float(z14_154[grp_idx].mean()),
        "donor_mean_stat1_log1p": float(stat1_154[grp_idx].mean()),
    })
df_g154_donors = pd.DataFrame(g154_donor_rows).sort_values(["lesion_site", "patient"])
df_g154_donors.to_csv(ST3_DIR / "gse154778_donor_abundance_primary_vs_met_Q1.tsv", sep="\t", index=False)

for min_m in [1, 10]:
    sub = df_g154_donors[df_g154_donors["n_macrophages"] >= min_m]
    met = sub[sub["lesion_site"] == "Metastasis"]
    pri = sub[sub["lesion_site"] == "Primary"]
    if len(met) >= 2 and len(pri) >= 2:
        for metric in ["ifn14_high_frac_of_mac_Q1", "donor_mean_ifn14_z", "donor_mean_stat1_log1p"]:
            x = met[metric].dropna().values
            y = pri[metric].dropna().values
            _, p_val = stats.mannwhitneyu(x, y, alternative="two-sided")
            hl, ci_lo, ci_hi = bootstrap_hl_diff(x, y)
            abundance_tests.append({
                "cohort": "GSE154778",
                "contrast": "Metastasis_vs_Primary_PDAC",
                "min_macrophages_per_donor": min_m,
                "metric": metric,
                "n_group1": len(x),
                "n_group2": len(y),
                "median_group1": float(np.median(x)),
                "iqr_group1": float(np.percentile(x, 75) - np.percentile(x, 25)),
                "median_group2": float(np.median(y)),
                "iqr_group2": float(np.percentile(y, 75) - np.percentile(y, 25)),
                "hodges_lehmann_diff": hl,
                "ci95_low": ci_lo,
                "ci95_high": ci_hi,
                "mwu_p": float(p_val),
            })

g154_extra = {r["patient"]: {"lesion_site": r["lesion_site"]} for _, r in df_g154_donors.iterrows()}
export_paired_pseudobulk(counts_154, mac_154.var_names, smp_154, q2_14_154, g154_extra, "GSE154778_All16_Q2_IFN14_min5", min_cells=5)

df_ab_tests = pd.DataFrame(abundance_tests)
df_ab_tests["mwu_fdr_within_contrast"] = np.nan
for (c_id, c_tr, mm), idx in df_ab_tests.groupby(["cohort", "contrast", "min_macrophages_per_donor"]).groups.items():
    df_ab_tests.loc[idx, "mwu_fdr_within_contrast"] = multipletests(df_ab_tests.loc[idx, "mwu_p"], method="fdr_bh")[1]
df_ab_tests.to_csv(ST3_DIR / "between_donor_abundance_tests_Q1.tsv", sep="\t", index=False)

# =========================================================================
# 4. PRE-SPECIFIED 39-GENE PANEL DONOR-PAIRED NONPARAMETRIC WILCOXON ACROSS ALL 3 COHORTS
# =========================================================================
print("[4/4] Computing donor-paired Wilcoxon effect sizes & 95% CIs on pre-specified functional panel across all 3 cohorts...", flush=True)
def paired_panel_stats(X_log, g2i, donor_vec, grp_vec, cohort_label, scheme_label, min_cells=5):
    rows = []
    donors = sorted(set(donor_vec))
    for cat, g in PRESPECIFIED_PANEL:
        if g not in g2i:
            continue
        vals = X_log[:, g2i[g]]
        hi_m, lo_m, hi_p, lo_p, used_d = [], [], [], [], []
        for d in donors:
            m_hi = (donor_vec == d) & (grp_vec == "high")
            m_lo = (donor_vec == d) & (grp_vec == "low")
            if m_hi.sum() >= min_cells and m_lo.sum() >= min_cells:
                hi_m.append(float(vals[m_hi].mean()))
                lo_m.append(float(vals[m_lo].mean()))
                hi_p.append(float((vals[m_hi] > 0).mean() * 100.0))
                lo_p.append(float((vals[m_lo] > 0).mean() * 100.0))
                used_d.append(d)
        hi_m = np.array(hi_m); lo_m = np.array(lo_m)
        hi_p = np.array(hi_p); lo_p = np.array(lo_p)
        diff_m = hi_m - lo_m
        diff_p = hi_p - lo_p
        if len(diff_m) >= 5 and np.any(np.abs(diff_m) > 1e-9):
            _, p_m = stats.wilcoxon(hi_m, lo_m, alternative="two-sided")
            rng = np.random.default_rng(42)
            b_med = [np.median(rng.choice(diff_m, size=len(diff_m), replace=True)) for _ in range(2000)]
            ci_lo, ci_hi = float(np.percentile(b_med, 2.5)), float(np.percentile(b_med, 97.5))
        else:
            p_m, ci_lo, ci_hi = 1.0, 0.0, 0.0
        if len(diff_p) >= 5 and np.any(np.abs(diff_p) > 1e-9):
            _, p_pct = stats.wilcoxon(hi_p, lo_p, alternative="two-sided")
        else:
            p_pct = 1.0
        rows.append({
            "cohort": cohort_label,
            "scheme": scheme_label,
            "category": cat,
            "gene": g,
            "in_IFN14": g in IFN14_PRIMARY,
            "in_IFN10": g in IFN10_STRICT,
            "n_paired_donors": len(diff_m),
            "donor_mean_high": float(hi_m.mean()) if len(hi_m) else np.nan,
            "donor_mean_low": float(lo_m.mean()) if len(lo_m) else np.nan,
            "median_paired_delta_log1p": float(np.median(diff_m)) if len(diff_m) else np.nan,
            "mean_paired_delta_log1p": float(diff_m.mean()) if len(diff_m) else np.nan,
            "ci95_delta_low": ci_lo,
            "ci95_delta_high": ci_hi,
            "pct_donors_positive_delta": float((diff_m > 0).mean() * 100.0) if len(diff_m) else np.nan,
            "donor_pct_pos_high": float(hi_p.mean()) if len(hi_p) else np.nan,
            "donor_pct_pos_low": float(lo_p.mean()) if len(lo_p) else np.nan,
            "wilcoxon_mean_p": float(p_m),
            "wilcoxon_pct_p": float(p_pct),
        })
    df_out = pd.DataFrame(rows)
    # Panel FDR across non-defining test genes!
    ind_mask = (~df_out["in_IFN14"]) if "IFN14" in scheme_label else (~df_out["in_IFN10"])
    df_out["wilcoxon_mean_padj_indep"] = np.nan
    df_out.loc[ind_mask, "wilcoxon_mean_padj_indep"] = multipletests(df_out.loc[ind_mask, "wilcoxon_mean_p"], method="fdr_bh")[1]
    df_out["wilcoxon_pct_padj_indep"] = np.nan
    df_out.loc[ind_mask, "wilcoxon_pct_padj_indep"] = multipletests(df_out.loc[ind_mask, "wilcoxon_pct_p"], method="fdr_bh")[1]
    return df_out

panel_dfs = [
    paired_panel_stats(Xlog_p, g2i_p, donors_p_vec, q2_14_p, "Peng2019_All", "Q2_IFN14_min5", 5),
    paired_panel_stats(Xlog_p[tum_mask_p], g2i_p, donors_p_vec[tum_mask_p], q2_14_p[tum_mask_p], "Peng2019_PDAConly", "Q2_IFN14_min5", 5),
    paired_panel_stats(Xlog_p, g2i_p, donors_p_vec, q2_10_p, "Peng2019_All", "Q2_IFN10_min5", 5),
    paired_panel_stats(Xlog_p, g2i_p, donors_p_vec, q1_14_p, "Peng2019_All", "Q1_IFN14_min5", 5),
    paired_panel_stats(Xlog_202, g2i_202, pids_202, q2_14_202, "GSE202051_All43", "Q2_IFN14_min5", 5),
    paired_panel_stats(Xlog_202[unt_mask], g2i_202, pids_202[unt_mask], q2_14_202[unt_mask], "GSE202051_Untreated18", "Q2_IFN14_min5", 5),
    paired_panel_stats(Xlog_202[neo_mask], g2i_202, pids_202[neo_mask], q2_14_202[neo_mask], "GSE202051_Neoadjuvant25", "Q2_IFN14_min5", 5),
    paired_panel_stats(Xlog_202, g2i_202, pids_202, q2_10_202, "GSE202051_All43", "Q2_IFN10_min5", 5),
    paired_panel_stats(Xlog_154, g2i_154, smp_154, q2_14_154, "GSE154778_All16", "Q2_IFN14_min5", 5),
]
df_panel_all = pd.concat(panel_dfs, ignore_index=True)
df_panel_all.to_csv(ST3_DIR / "paired_donor_functional_panel_all_cohorts.tsv", sep="\t", index=False)

# Also save donor-by-gene paired matrices for Peng & GSE202051 for transparent slopegraph / boxplot plotting!
slope_genes = ["STAT1", "IRF1", "LGALS9", "CD47", "CD80", "NECTIN2", "CD274", "PDCD1LG2", "IDO1", "HLA-A", "B2M", "TAP1", "CXCL9", "CXCL10", "TGFB1", "KRT19", "COL1A1"]
donor_pair_records = []
for coh_lbl, X_m, g2i_m, d_v, g_v, extra_v in [
    ("Peng2019", Xlog_p, g2i_p, donors_p_vec, q2_14_p, peng_extra),
    ("GSE202051", Xlog_202, g2i_202, pids_202, q2_14_202, g202_extra),
    ("GSE154778", Xlog_154, g2i_154, smp_154, q2_14_154, g154_extra),
]:
    for d in sorted(set(d_v)):
        m_hi = (d_v == d) & (g_v == "high")
        m_lo = (d_v == d) & (g_v == "low")
        if m_hi.sum() >= 5 and m_lo.sum() >= 5:
            sub_info = extra_v.get(d, {})
            sub_tag = sub_info.get("condition") or sub_info.get("treatment_binary") or sub_info.get("lesion_site") or ""
            for sg in slope_genes:
                if sg in g2i_m:
                    vh = X_m[m_hi, g2i_m[sg]]
                    vl = X_m[m_lo, g2i_m[sg]]
                    donor_pair_records.append({
                        "cohort": coh_lbl,
                        "patient": d,
                        "stratum": sub_tag,
                        "gene": sg,
                        "n_high": int(m_hi.sum()),
                        "n_low": int(m_lo.sum()),
                        "mean_high": float(vh.mean()),
                        "mean_low": float(vl.mean()),
                        "delta_mean": float(vh.mean() - vl.mean()),
                        "pct_pos_high": float((vh > 0).mean() * 100.0),
                        "pct_pos_low": float((vl > 0).mean() * 100.0),
                    })
pd.DataFrame(donor_pair_records).to_csv(ST3_DIR / "donor_level_paired_gene_values.tsv", sep="\t", index=False)

print("\n=== SUMMARY: Peng2019 Q1 Between-Donor Abundance ===", flush=True)
print(df_ab_tests[df_ab_tests["cohort"] == "Peng2019"][["min_macrophages_per_donor", "metric", "median_group1", "median_group2", "hodges_lehmann_diff", "mwu_p"]].to_string(index=False), flush=True)
print("\n=== SUMMARY: Key Checkpoints & Controls Across Cohorts (Q2_IFN14_min5) ===", flush=True)
sub_show = df_panel_all[(df_panel_all["scheme"] == "Q2_IFN14_min5") & (df_panel_all["gene"].isin(["STAT1", "IRF1", "LGALS9", "CD47", "CD80", "NECTIN2", "CD274", "PDCD1LG2", "IDO1", "B2M", "TAP1", "KRT19", "COL1A1"]))]
print(sub_show[["cohort", "gene", "n_paired_donors", "donor_mean_high", "donor_mean_low", "pct_donors_positive_delta", "wilcoxon_mean_p", "wilcoxon_mean_padj_indep"]].to_string(index=False), flush=True)
print("[Stage 3 Python Part Complete]", flush=True)
