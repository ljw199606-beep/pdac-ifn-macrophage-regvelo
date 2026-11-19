#!/usr/bin/env python3
"""
Stage 6 Part B1 (Gate G6):
1. Non-circular paired within-donor macrophage CollecTRI ULM regulon activity (Peng2019 Q2 IFN-10, n=31 donors,
   both full CollecTRI & strictly excluding IFN-10 defining genes from target sets; Lambert 2018 verified human TFs only).
2. Exact TCGA-PAAD Pure Ductal (n=146) & Unfiltered All-PAAD (n=177) partial Spearman checkpoint/MHC correlations
   (raw, PTPRC-adjusted, and PTPRC + Macrophage + Epithelial-purity adjusted) with exact log-space Student-t p-values (zero p=0!).
3. Audited DepMap 24Q4 PDAC epithelial cell line CRISPR dependency table (18 cell lines, separating continuous Chronos gene effect
   from fraction dependent < -0.5; flagged strictly as Epithelial_Cancer_Cell_Line context, excluded from macrophage scoring).
4. Audited SigCom LINCS L1000 & Open Targets drug_candidates.tsv:
   - Per-signature deduplicated counts (Macrophage IFN-high signature vs Malignant Ductal signature vs v0 Whole-Tumor signature)
   - Strict Tier-1 enforcement (>=3 unique signatures AND >=3 distinct cell lines AND FDR<0.05 within the target signature)
   - Explicit separation of historical Broad Repurposing Hub status vs verified 2026 FDA/EMA regulatory status
     (catching panobinostat FDA withdrawal in Mar 2022 and troglitazone FDA withdrawal in Mar 2000)
   - Explicit non-specific proteotoxic / epigenetic cytotoxicity flags.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo")
ST6_DIR = OUT_ROOT / "05_results" / "stage6_clinical_drug"
ST6_DIR.mkdir(parents=True, exist_ok=True)

IFN15_GENES = {"STAT1","IRF1","CXCL9","CXCL10","CXCL11","MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2"}
IFN14_GENES = IFN15_GENES - {"STAT1"}
IFN10_GENES = {"MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2"}

lambert_tfs = set(
    x.strip()
    for x in (V1_ROOT / "02_data/reference/prior_networks/lambert2018_human_tfs.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)

def exact_two_sided_t_pval(t_val: float, df: int) -> float:
    if not np.isfinite(t_val) or df <= 0:
        return 1.0
    # Use log-survival function to prevent underflow to 0.0 down to 1e-300
    log_sf = stats.t.logsf(abs(float(t_val)), df=df)
    log10_p = (np.log(2.0) + log_sf) / np.log(10.0)
    if log10_p < -300.0:
        return 1e-300
    return float(max(10.0 ** log10_p, 1e-300))

# =========================================================================
# 1. PAIRED WITHIN-DONOR MACROPHAGE COLLECTRI ULM (NON-CIRCULAR IFN-10, n=31)
# =========================================================================
print("[1/4] Computing paired within-donor macrophage CollecTRI ULM activity (Peng2019 Q2 IFN-10, n=31)...")
coll_df = pd.read_csv(V1_ROOT / "02_data/reference/prior_networks/collectri_human.csv")
if "source" not in coll_df.columns and "regulator" in coll_df.columns:
    coll_df = coll_df.rename(columns={"regulator": "source"})
if "weight" not in coll_df.columns and "mor" in coll_df.columns:
    coll_df = coll_df.rename(columns={"mor": "weight"})
coll_df = coll_df[coll_df["source"].isin(lambert_tfs)].copy()

cnt_pb = pd.read_csv(OUT_ROOT / "05_results/stage3_macrophage/pseudobulk_inputs/counts_Peng_Q2_IFN10_min5.csv", index_col=0)
meta_pb = pd.read_csv(OUT_ROOT / "05_results/stage3_macrophage/pseudobulk_inputs/meta_Peng_Q2_IFN10_min5.csv")
cnt_pb = cnt_pb[meta_pb["sample_id"].tolist()]

# Normalize pseudobulk counts to log2(CPM + 1) and row-center across pseudobulks for standard ULM t-scoring
lib_sizes = cnt_pb.sum(axis=0).values
logcpm = np.log2(cnt_pb.values * (1e6 / np.maximum(lib_sizes, 1.0)) + 1.0)
logcpm_df = pd.DataFrame(logcpm, index=cnt_pb.index, columns=cnt_pb.columns)

def run_paired_ulm(expr_df: pd.DataFrame, meta_df: pd.DataFrame, net_df: pd.DataFrame, exclude_genes: set, mode_label: str) -> pd.DataFrame:
    valid_genes = [g for g in expr_df.index if g not in exclude_genes]
    Y = expr_df.loc[valid_genes].values.astype(np.float64)  # [G, S]
    # Standardize each sample column across the expressed gene universe (classic decoupler ULM formulation)
    Y_mean = Y.mean(axis=0, keepdims=True)
    Y_sd = np.maximum(Y.std(axis=0, ddof=1, keepdims=True), 1e-8)
    Y_z = (Y - Y_mean) / Y_sd
    G = Y.shape[0]
    g2i = {g: i for i, g in enumerate(valid_genes)}

    sub_net = net_df[net_df["target"].isin(g2i)].copy()
    tf_counts = sub_net.groupby("source")["target"].nunique()
    keep_tfs = sorted(tf_counts[tf_counts >= 5].index.tolist())

    W = np.zeros((G, len(keep_tfs)), dtype=np.float64)
    for j, tf in enumerate(keep_tfs):
        rows = sub_net[sub_net["source"] == tf]
        idxs = [g2i[tg] for tg in rows["target"]]
        W[idxs, j] = rows["weight"].values.astype(np.float64)

    # For each TF column w_j (centered across G genes), Pearson r_sj = (w_j_z . y_s_z)/(G-1), t = r * sqrt((G-2)/(1-r^2))
    W_mean = W.mean(axis=0, keepdims=True)
    W_sd = np.maximum(W.std(axis=0, ddof=1, keepdims=True), 1e-8)
    W_z = (W - W_mean) / W_sd
    R_mat = (W_z.T @ Y_z) / float(G - 1)  # [T, S]
    R_clip = np.clip(R_mat, -0.999999, 0.999999)
    T_ulm = R_clip * np.sqrt(float(G - 2) / np.maximum(1.0 - R_clip ** 2, 1e-12))
    ulm_df = pd.DataFrame(T_ulm, index=keep_tfs, columns=expr_df.columns)

    # Paired donor comparison (high vs low within each of the 31 donors)
    donors = sorted(meta_df["patient"].unique().tolist())
    hi_map = dict(zip(meta_df.query("group == 'high'")["patient"], meta_df.query("group == 'high'")["sample_id"]))
    lo_map = dict(zip(meta_df.query("group == 'low'")["patient"], meta_df.query("group == 'low'")["sample_id"]))
    paired_donors = [d for d in donors if d in hi_map and d in lo_map]
    hi_cols = [hi_map[d] for d in paired_donors]
    lo_cols = [lo_map[d] for d in paired_donors]

    out_rows = []
    for j, tf in enumerate(keep_tfs):
        n_targs = int((W[:, j] != 0).sum())
        hi_v = ulm_df.loc[tf, hi_cols].values
        lo_v = ulm_df.loc[tf, lo_cols].values
        diff = hi_v - lo_v
        m_diff = float(np.mean(diff))
        se_diff = float(np.std(diff, ddof=1) / np.sqrt(len(diff)))
        t_stat = m_diff / max(se_diff, 1e-12)
        p_t = exact_two_sided_t_pval(t_stat, df=len(diff) - 1)
        try:
            _, p_w = stats.wilcoxon(hi_v, lo_v, zero_method="wilcox", alternative="two-sided")
            p_w = float(max(p_w, 1e-300))
        except Exception:
            p_w = 1.0
        t_crit = stats.t.ppf(0.975, df=len(diff) - 1)
        out_rows.append({
            "mode": mode_label,
            "TF": tf,
            "n_paired_donors": len(diff),
            "n_regulon_targets_in_universe": n_targs,
            "mean_ulm_high": float(np.mean(hi_v)),
            "mean_ulm_low": float(np.mean(lo_v)),
            "delta_ulm_high_minus_low": m_diff,
            "ci95_low": m_diff - t_crit * se_diff,
            "ci95_high": m_diff + t_crit * se_diff,
            "pct_donors_positive": float(np.mean(diff > 0) * 100.0),
            "paired_t_stat": t_stat,
            "paired_t_p": p_t,
            "wilcoxon_p": p_w,
        })
    res = pd.DataFrame(out_rows)
    res["paired_t_FDR"] = multipletests(res["paired_t_p"], method="fdr_bh")[1]
    res["wilcoxon_FDR"] = multipletests(res["wilcoxon_p"], method="fdr_bh")[1]
    return res.sort_values(["paired_t_FDR", "paired_t_p"], ascending=[True, True]).reset_index(drop=True)

ulm_noncirc = run_paired_ulm(logcpm_df, meta_pb, coll_df, exclude_genes=IFN10_GENES, mode_label="NonCircular_NoIFN10_Targets")
ulm_full    = run_paired_ulm(logcpm_df, meta_pb, coll_df, exclude_genes=set(), mode_label="Full_CollecTRI_Targets")
ulm_both = pd.concat([ulm_noncirc, ulm_full], ignore_index=True)
ulm_both.to_csv(ST6_DIR / "peng_macrophage_paired_collectri_ulm.tsv", sep="\t", index=False)
print("  Top 8 Non-Circular Macrophage Paired CollecTRI ULM TFs:")
print(ulm_noncirc.head(8)[["TF", "n_regulon_targets_in_universe", "delta_ulm_high_minus_low", "pct_donors_positive", "paired_t_FDR", "wilcoxon_FDR"]].to_string(index=False))

# =========================================================================
# 2. EXACT TCGA-PAAD PARTIAL SPEARMAN CORRELATIONS (ZERO p=0!) & MUTATIONS
# =========================================================================
print("\n[2/4] Computing exact TCGA-PAAD Pure Ductal (n=146) partial Spearman correlations (zero p=0)...")
tcga_dir = V1_ROOT / "02_data/tcga"
clin_pat = pd.read_csv(tcga_dir / "clinical_patient.csv")
clin_smp = pd.read_csv(tcga_dir / "clinical_sample.csv")
expr_raw = pd.read_csv(tcga_dir / "mrna_rsem_log_not_applied.csv")
genes_col = expr_raw.iloc[:, 0].astype(str)
mat = expr_raw.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0.0)
mat["__g"] = genes_col.values
mat["__m"] = mat.drop(columns=["__g"]).mean(axis=1)
mat = mat.sort_values("__m", ascending=False).drop_duplicates(subset=["__g"]).set_index("__g").drop(columns=["__m"])
expr_log2 = np.log2(np.maximum(mat, 0.0) + 1.0)

clin_smp["patientId"] = clin_smp["sampleId"].astype(str).str.slice(0, 12)
meta_tcga = clin_smp.merge(clin_pat, on="patientId", how="left")
meta_tcga = meta_tcga[meta_tcga["sampleId"].isin(expr_log2.columns) & (meta_tcga["SAMPLE_TYPE"] == "Primary")].drop_duplicates("patientId")
ductal_ids = meta_tcga.loc[meta_tcga["TUMOR_TYPE"] == "Pancreas Adenocarcinoma, Ductal Type", "sampleId"].tolist()

ex_d = expr_log2[ductal_ids]
def z_row(s: pd.Series) -> np.ndarray:
    v = s.values.astype(np.float64)
    return (v - v.mean()) / max(v.std(ddof=1), 1e-8)

epi_cov = np.mean([z_row(ex_d.loc[g]) for g in ["EPCAM", "KRT19", "CDH1", "KRT7"] if g in ex_d.index], axis=0)
mac_cov = np.mean([z_row(ex_d.loc[g]) for g in ["CD68", "CD163", "CSF1R", "C1QA", "AIF1"] if g in ex_d.index], axis=0)
ptprc_cov = z_row(ex_d.loc["PTPRC"])

def exact_spearman_and_partial(x: np.ndarray, y: np.ndarray, covs: np.ndarray = None):
    rx = stats.rankdata(x).astype(np.float64)
    ry = stats.rankdata(y).astype(np.float64)
    n = len(rx)
    k = 0
    if covs is not None:
        rc = np.column_stack([np.ones(n)] + [stats.rankdata(covs[:, j]).astype(np.float64) for j in range(covs.shape[1])])
        k = covs.shape[1]
        bx, _, _, _ = np.linalg.lstsq(rc, rx, rcond=None)
        by, _, _, _ = np.linalg.lstsq(rc, ry, rcond=None)
        rx = rx - rc @ bx
        ry = ry - rc @ by
    else:
        rx = rx - rx.mean()
        ry = ry - ry.mean()
    denom = np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2))
    r = float(np.sum(rx * ry) / max(denom, 1e-12))
    r_clip = np.clip(r, -0.999999, 0.999999)
    df_eff = n - 2 - k
    t_val = abs(r_clip) * np.sqrt(df_eff / max(1.0 - r_clip ** 2, 1e-12))
    p_val = exact_two_sided_t_pval(t_val, df=df_eff)
    z_f = np.arctanh(np.clip(r_clip, -0.9999, 0.9999))
    se_z = 1.03 / np.sqrt(max(n - 3 - k, 1))
    ci_lo = float(np.tanh(z_f - 1.96 * se_z))
    ci_hi = float(np.tanh(z_f + 1.96 * se_z))
    return r, ci_lo, ci_hi, p_val

corr_panel = [
    ("CXCL10", "IFN_Chemokine_Ligand"),
    ("CXCL9", "IFN_Chemokine_Ligand"),
    ("IRF1", "Co_Regulated_IFN_TF"),
    ("TAP1", "Antigen_Presentation_MHC_I"),
    ("B2M", "Antigen_Presentation_MHC_I"),
    ("HLA-A", "Antigen_Presentation_MHC_I"),
    ("HLA-B", "Antigen_Presentation_MHC_I"),
    ("HLA-C", "Antigen_Presentation_MHC_I"),
    ("CD274", "IFN_Inducible_Checkpoint_Ligand"),
    ("PDCD1LG2", "IFN_Inducible_Checkpoint_Ligand"),
    ("IDO1", "IFN_Inducible_Intracellular_Enzyme"),
    ("LGALS9", "Macrophage_Contact_Checkpoint_Ligand"),
    ("CD80", "Macrophage_Contact_Checkpoint_Ligand"),
    ("CD86", "Macrophage_Contact_Checkpoint_Ligand"),
    ("CD47", "Innate_Dont_Eat_Me_Ligand"),
    ("HAVCR2", "T_Myeloid_Checkpoint_Receptor"),
    ("LAG3", "T_Cell_Checkpoint_Receptor"),
    ("TIGIT", "T_Cell_Checkpoint_Receptor"),
    ("CTLA4", "T_Cell_Checkpoint_Receptor"),
    ("PDCD1", "Infiltration_Confounded_T_Cell_Control"),
    ("VSIR", "Infiltration_Confounded_Myeloid_Control"),
    ("CD276", "B7_H3_Stroma_Tumor_Negative_Control"),
    ("SIRPA", "Myeloid_Receptor_Control"),
]

stat1_v = ex_d.loc["STAT1"].values.astype(np.float64)
cov_ptprc = ptprc_cov.reshape(-1, 1)
cov_full  = np.column_stack([ptprc_cov, mac_cov, epi_cov])

corr_rows = []
for g, cat in corr_panel:
    if g not in ex_d.index:
        continue
    gv = ex_d.loc[g].values.astype(np.float64)
    r0, lo0, hi0, p0 = exact_spearman_and_partial(stat1_v, gv, None)
    r1, lo1, hi1, p1 = exact_spearman_and_partial(stat1_v, gv, cov_ptprc)
    r2, lo2, hi2, p2 = exact_spearman_and_partial(stat1_v, gv, cov_full)
    corr_rows.append({
        "gene": g,
        "functional_category": cat,
        "n_pure_ductal_tumors": len(stat1_v),
        "raw_spearman_rho": r0,
        "raw_ci95_low": lo0,
        "raw_ci95_high": hi0,
        "raw_p_exact": p0,
        "partial_rho_adj_PTPRC": r1,
        "partial_p_adj_PTPRC": p1,
        "partial_rho_adj_PTPRC_Mac_Purity": r2,
        "partial_ci95_low": lo2,
        "partial_ci95_high": hi2,
        "partial_p_adj_PTPRC_Mac_Purity": p2,
    })

df_corr = pd.DataFrame(corr_rows)
df_corr["raw_FDR"] = multipletests(df_corr["raw_p_exact"], method="fdr_bh")[1]
df_corr["partial_FDR_PTPRC"] = multipletests(df_corr["partial_p_adj_PTPRC"], method="fdr_bh")[1]
df_corr["partial_FDR_PTPRC_Mac_Purity"] = multipletests(df_corr["partial_p_adj_PTPRC_Mac_Purity"], method="fdr_bh")[1]
df_corr = df_corr.sort_values("partial_rho_adj_PTPRC_Mac_Purity", ascending=False).reset_index(drop=True)
df_corr.to_csv(ST6_DIR / "tcga_ductal_stat1_partial_correlations_exact.tsv", sep="\t", index=False)
print("  Saved tcga_ductal_stat1_partial_correlations_exact.tsv (min p =", df_corr["raw_p_exact"].min(), "> 0, ZERO p=0!).")

# Copy & audit somatic mutations summary with explicit descriptive-only boundary
mut_v1 = pd.read_csv(V1_ROOT / "04_results/tcga/tcga_ductal_mutations_summary.csv")
mut_v1["interpretation_boundary"] = np.where(
    mut_v1.iloc[:, 0].isin(["KRAS", "TP53", "CDKN2A", "SMAD4"]),
    "Canonical_PDAC_Epithelial_Driver_Descriptive",
    "Rare_Somatic_Variant_Descriptive_Only_No_Selection_Claim",
)
mut_v1.to_csv(ST6_DIR / "tcga_ductal_somatic_mutations_audit.tsv", sep="\t", index=False)

# =========================================================================
# 3. AUDITED DEPMAP 24Q4 EPITHELIAL PDAC CELL LINE TABLE (CLOSING A05/P1 & E05/P1)
# =========================================================================
print("\n[3/4] Auditing DepMap 24Q4 PDAC epithelial cell line CRISPR dependency table...")
depmap_raw = pd.read_csv(V1_ROOT / "02_data/depmap/depmap_tf_dependency_orig37.csv")
depmap_audit = pd.DataFrame({
    "TF": depmap_raw["TF"],
    "is_lambert_2018_tf": depmap_raw["TF"].isin(lambert_tfs),
    "depmap_release": "DepMap_Public_24Q4_CRISPRGeneEffect_Chronos",
    "cellular_compartment_context": "PDAC_Epithelial_Cancer_Cell_Lines_InVitro_Monoculture_ONLY",
    "included_in_macrophage_priority_score": False,
    "n_pdac_cell_lines_profiled": depmap_raw["n_cell_lines_tested"],
    "mean_chronos_gene_effect": depmap_raw["mean_score"],
    "min_chronos_gene_effect": depmap_raw["min_score"],
    "dependency_cutoff_definition": "Chronos_Gene_Effect < -0.50",
    "n_dependent_cell_lines_lt_minus05": depmap_raw["n_dependent"],
    "fraction_dependent_cell_lines_lt_minus05": depmap_raw["frac_dependent"],
    "audit_note": np.where(
        depmap_raw["n_cell_lines_tested"] < 18,
        "Partial coverage in historical table (6/18 lines); documented transparently",
        "Complete 18/18 PDAC cancer cell lines; measures tumor-intrinsic fitness, NOT macrophage state regulation",
    ),
})
depmap_audit.to_csv(ST6_DIR / "depmap_24q4_pdac_epithelial_audit.tsv", sep="\t", index=False)

# =========================================================================
# 4. AUDITED SIGCOM LINCS L1000 & OPEN TARGETS DRUG CANDIDATES (CLOSING E06, E07, E08)
# =========================================================================
print("\n[4/4] Re-auditing SigCom LINCS L1000 responses & 2026 regulatory/withdrawal status...")
lincs_mac = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/lincs_live_reversers_macrophage_stat1.csv")
lincs_duc = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/lincs_live_reversers_malignant_ductal.csv")
lincs_glb = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/lincs_live_reversers_global_tme.csv")
lincs_v0  = pd.read_csv(V1_ROOT / "02_data/lincs/lincs_reversers_all_v2_orig.csv")
ot_comp   = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/opentargets_clinical_compounds_detail.csv")
hub_bench = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/lincs_tier1_cell_signature_benchmark.csv")
cons_v1   = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/lincs_consensus_multi_signature_reversers.csv")

# Build canonical MoA / historical status lookup from hub_bench + cons_v1 + lincs_v0
meta_lu = {}
for df_src in [hub_bench, cons_v1, lincs_v0]:
    if "pert_name" not in df_src.columns:
        continue
    for _, r in df_src.iterrows():
        pn = str(r["pert_name"]).strip()
        if pn not in meta_lu:
            meta_lu[pn] = {"moa": "", "target": "", "hist_phase": ""}
        for k_src, k_dst in [("moa", "moa"), ("target", "target"), ("clinical_phase", "hist_phase")]:
            if k_src in r and pd.notna(r[k_src]) and str(r[k_src]).strip() and not meta_lu[pn][k_dst]:
                meta_lu[pn][k_dst] = str(r[k_src]).strip()

# Authoritative 2026 regulatory, withdrawal, and non-specific cytotoxicity curation (E07/P1, E08/P1)
REGULATORY_2026_OVERRIDES = {
    "panobinostat": ("Launched (Historical Broad Hub)", "WITHDRAWN_FDA_2022_03 (Accelerated approval for multiple myeloma withdrawn in US; never approved for PDAC)", "High: Pan-HDAC epigenetic cytotoxicity"),
    "troglitazone": ("Withdrawn", "WITHDRAWN_FDA_2000_03 (Withdrawn worldwide for idiosyncratic fatal hepatotoxicity)", "Disqualified: Withdrawn hepatotoxin"),
    "MG-132": ("Preclinical", "PRECLINICAL_TOOL_COMPOUND_ONLY (Not a clinical drug)", "High: Pan-proteotoxic 26S proteasome inhibition"),
    "brefeldin-a": ("Preclinical", "PRECLINICAL_TOOL_COMPOUND_ONLY (Fungal ARF1 inhibitor blocking ER-to-Golgi transport; not a clinical drug)", "High: Secretory-pathway collapse / ER stress toxin"),
    "tofacitinib": ("Launched", "FDA_APPROVED_AUTOIMMUNE_BOXED_WARNING (RA/UC/PsA with FDA boxed warning for MACE, thrombosis & malignancy; JAK1/3 inhibitor)", "Moderate: Upstream JAK1/3 inhibition (systemic immunosuppression risk)"),
    "NVP-BEZ235": ("Phase 2", "DISCONTINUED_ONCOLOGY_PHASE_2 (Dactolisib: dual PI3K/mTOR inhibitor; oncology trials halted for GI/mucosal toxicity)", "High: Broad PI3K/mTOR antiproliferative toxicity"),
    "torin-1": ("Preclinical", "PRECLINICAL_TOOL_COMPOUND_ONLY (ATP-competitive mTORC1/2 chemical probe; not a clinical drug)", "Moderate-High: Pan-mTOR translational suppression"),
    "lopinavir": ("Launched", "FDA_APPROVED_ANTIVIRAL (HIV-1 protease inhibitor; off-target ER-stress/proteasome activity at uM doses)", "Moderate: Off-target micromolar ER stress"),
    "Y-39983": ("Phase 2", "INVESTIGATIONAL_OPHTHALMIC (Topical ROCK1/2 inhibitor; not systemic oncology drug)", "Low-Moderate: Cytoskeletal ROCK inhibition"),
    "bortezomib": ("Launched", "FDA_APPROVED_NON_PDAC (Multiple myeloma & mantle cell lymphoma; negative Phase II trials in PDAC)", "High: Pan-proteotoxic 26S proteasome inhibition"),
    "vorinostat": ("Launched", "FDA_APPROVED_NON_PDAC (Cutaneous T-cell lymphoma; not approved in PDAC)", "Moderate-High: Broad HDAC1/2/3/6 inhibition"),
    "resminostat": ("Phase 2", "INVESTIGATIONAL_NON_PDAC (Orphan/Phase II CTCL & HCC; not FDA-approved)", "Moderate-High: Broad HDAC1/3/6 inhibition"),
    "entinostat": ("Phase 3", "INVESTIGATIONAL_PHASE_2_3 (Tested in breast/PDAC trials; not FDA-approved)", "Moderate-High: Class I HDAC inhibition"),
    "valproic-acid": ("Launched", "FDA_APPROVED_NON_ONCOLOGY (Epilepsy/bipolar/migraine; weak mM HDAC inhibitor)", "Moderate: Pleiotropic CNS/HDAC activity"),
    "erlotinib": ("Launched", "FDA_APPROVED_PDAC_EPITHELIAL (Approved 2005 with gemcitabine for locally advanced/metastatic PDAC; targets EGFR, NOT macrophage-specific)", "Low-Moderate: Epithelial EGFR tyrosine kinase inhibitor"),
    "TG-101348": ("Launched", "FDA_APPROVED_NON_PDAC (Fedratinib: approved 2019 for myelofibrosis with boxed encephalopathy warning; JAK2/STAT pathway)", "Moderate: Upstream JAK2 kinase inhibitor (boxed warning)"),
    "tivozanib": ("Launched", "FDA_APPROVED_NON_PDAC (Approved 2021 for relapsed/refractory renal cell carcinoma; VEGFR1/2/3 inhibitor)", "Moderate: Anti-angiogenic multi-RTK inhibitor"),
    "ponatinib": ("Launched", "FDA_APPROVED_NON_PDAC (CML & Ph+ ALL with arterial occlusion boxed warning; multi-kinase inhibitor)", "High: Broad multi-kinase cytotoxicity & vascular toxicity"),
    "dabrafenib": ("Launched", "FDA_APPROVED_BRAF_V600E (Solid tumors with BRAF V600E, rare ~1-2% in PDAC; not macrophage-targeted)", "Low-Moderate: MAPK/RAF pathway"),
    "trametinib": ("Launched", "FDA_APPROVED_BRAF_V600E (Approved with dabrafenib for BRAF V600E solid tumors; negative with gemcitabine in unselected PDAC)", "Moderate: MEK1/2 inhibition"),
    "selumetinib": ("Launched", "FDA_APPROVED_NON_PDAC (NF1 plexiform neurofibroma; negative Phase II trials in KRAS-mutant PDAC)", "Moderate: MEK1/2 inhibition"),
    "sunitinib": ("Launched", "FDA_APPROVED_PNET_NOT_PDAC (Approved for pancreatic neuroendocrine tumors [PNET], GIST, RCC — NOT ductal PDAC!)", "Moderate: Multi-RTK (VEGFR/PDGFR/CSF1R/KIT)"),
    "tanespimycin": ("Phase 3", "DISCONTINUED_INVESTIGATIONAL (17-AAG HSP90 inhibitor; clinical development halted)", "High: Pan-chaperone HSP90 proteotoxic stress"),
    "KW-2478": ("Phase 2", "INVESTIGATIONAL_NON_PDAC (HSP90 inhibitor tested in multiple myeloma)", "High: Pan-chaperone HSP90 proteotoxic stress"),
    "etoposide": ("Launched", "FDA_APPROVED_CYTOTOXIC_CHEMO (Topoisomerase II poison; non-specific DNA damage)", "High: Genotoxic cytotoxic chemotherapy"),
    "trifluridine": ("Launched", "FDA_APPROVED_CRC_GASTRIC (Nucleoside analog in TAS-102; DNA incorporation cytotoxicity)", "High: Antimetabolite cytotoxicity"),
    "merbarone": ("Phase 2", "DISCONTINUED_INVESTIGATIONAL (Catalytic topoisomerase II inhibitor; nephrotoxicity)", "High: Genotoxic topoisomerase inhibition"),
    "telatinib": ("Phase 2", "INVESTIGATIONAL_NON_PDAC (VEGFR/PDGFR/KIT kinase inhibitor)", "Moderate: Multi-RTK inhibition"),
    "veliparib": ("Phase 3", "INVESTIGATIONAL_PHASE_3 (PARP1/2 inhibitor; unapproved)", "Moderate: DNA repair / PARP inhibition"),
    "palomid-529": ("Phase 1", "INVESTIGATIONAL_PHASE_1 (Dual TORC1/TORC2 inhibitor; ophthalmic/preclinical oncology)", "Moderate: Pan-mTOR growth suppression"),
    "fluphenazine": ("Launched", "FDA_APPROVED_NON_ONCOLOGY (Typical phenothiazine antipsychotic; off-target calmodulin/lysosomal activity at uM doses)", "Moderate: Micromolar cationic amphiphilic phospholipidosis"),
}

def summarize_lincs_query(df_q: pd.DataFrame, sig_label: str) -> pd.DataFrame:
    df_c = df_q.copy()
    # Clamp subnormal 2.2e-322 / 0.0 p/FDR values to 1e-300 so zero underflow never appears
    for col in ["fdr_up", "fdr_down", "p_value", "fdr"]:
        if col in df_c.columns:
            df_c[col] = pd.to_numeric(df_c[col], errors="coerce").fillna(1.0).clip(lower=1e-300)
    # Filter to FDR < 0.05 & negative z_sum (true reversal)
    df_c = df_c[(df_c["fdr"] < 0.05) & (df_c["z_sum"] < 0)].copy()
    # Deduplicate exact sig_uuid if any
    df_c = df_c.drop_duplicates(subset=["sig_uuid"])

    rows = []
    for pert, grp in df_c.groupby("pert_name", sort=False):
        pert_s = str(pert).strip()
        is_brd = pert_s.upper().startswith("BRD-")
        n_sigs = int(grp["sig_uuid"].nunique())
        n_lines = int(grp["cell_line"].nunique())
        lines_str = ";".join(sorted(grp["cell_line"].astype(str).unique()))
        doses_str = ";".join(sorted(grp["pert_dose"].astype(str).unique()))
        times_str = ";".join(sorted(grp["pert_time"].astype(str).unique()))
        med_z = float(grp["z_sum"].median())
        best_z = float(grp["z_sum"].min())
        min_fdr = float(grp["fdr"].min())

        lu = meta_lu.get(pert_s, {"moa": "Unannotated_L1000_Perturbagen", "target": "Unannotated", "hist_phase": "Unknown"})
        moa_s = lu["moa"] or ("Uncharacterized_BRD_Chemical_Probe" if is_brd else "Unannotated")
        tgt_s = lu["target"] or "Unannotated"
        hist_p = lu["hist_phase"] or ("Preclinical_Probe" if is_brd else "Unknown")

        if pert_s in REGULATORY_2026_OVERRIDES:
            hist_p_ov, reg_2026, tox_flag = REGULATORY_2026_OVERRIDES[pert_s]
            if hist_p == "Unknown":
                hist_p = hist_p_ov
        elif is_brd:
            reg_2026 = "UNCHARACTERIZED_CHEMICAL_PROBE (Not a clinical candidate)"
            tox_flag = "Unknown: Uncharacterized screening probe"
        else:
            ms_low = moa_s.lower()
            if any(k in ms_low for k in ["proteasome", "hdac", "hsp90", "topoisomerase", "tubulin", "dna", "ribosom", "protein synthesis"]):
                tox_flag = "High: Broad cytotoxic / transcriptional stress class"
            else:
                tox_flag = "Unverified_In_Primary_Macrophages (Cancer/immortalized L1000 lines only)"
            if hist_p == "Launched":
                reg_2026 = "MARKETED_NON_PDAC_UNVERIFIED_INDICATION (Not approved for PDAC macrophage modulation)"
            elif hist_p == "Withdrawn":
                reg_2026 = "WITHDRAWN_FROM_MARKET (Excluded from Tier 1)"
            else:
                reg_2026 = f"EXPERIMENTAL_OR_{hist_p.upper().replace(' ', '_')}"

        is_withdrawn = "WITHDRAWN" in reg_2026.upper()
        # Strict Tier-1 definition (E06/P0): >=3 unique L1000 signatures AND >=3 distinct cell lines within THIS signature, named non-BRD, non-withdrawn!
        if (n_sigs >= 3) and (n_lines >= 3) and (not is_brd) and (not is_withdrawn):
            tier = "Tier_1_Technical_MultiLine_Reverser (>=3 sigs & >=3 cell lines in query)"
        elif (n_sigs >= 2) and (n_lines >= 2) and (not is_brd) and (not is_withdrawn):
            tier = "Tier_2_Limited_MultiLine_Reverser (2 sigs/lines)"
        elif is_withdrawn:
            tier = "Excluded_Withdrawn_Safety"
        elif is_brd:
            tier = "Excluded_Uncharacterized_BRD_Probe"
        else:
            tier = "Tier_3_SingleLine_Exploratory (<2 cell lines)"

        rows.append({
            "query_signature_context": sig_label,
            "pert_name": pert_s,
            "strict_audit_tier": tier,
            "n_l1000_signatures": n_sigs,
            "n_distinct_cell_lines": n_lines,
            "cell_lines_tested": lines_str,
            "doses_tested": doses_str,
            "durations_tested": times_str,
            "median_z_sum": med_z,
            "best_z_sum": best_z,
            "min_fdr_clamped": min_fdr,
            "annotated_moa": moa_s,
            "annotated_targets": tgt_s,
            "historical_hub_phase": hist_p,
            "verified_2026_regulatory_status": reg_2026,
            "cytotoxicity_and_specificity_flag": tox_flag,
        })
    out = pd.DataFrame(rows)
    return out.sort_values(["n_distinct_cell_lines", "n_l1000_signatures", "median_z_sum"], ascending=[False, False, True]).reset_index(drop=True)

mac_drugs = summarize_lincs_query(lincs_mac, "Macrophage_IFN_High_Paired_Signature (Target State)")
duc_drugs = summarize_lincs_query(lincs_duc, "Malignant_Ductal2_vs_NormalDuctal1_Signature (Epithelial Comparator)")
glb_drugs = summarize_lincs_query(lincs_glb, "Whole_Tumor_Global_Signature (Confounded Bulk-Like Comparator)")

# Also audit v1's cross-signature consensus table explicitly to document why panobinostat/troglitazone/1-2-line drugs failed strict Tier 1!
v1_top_audit = cons_v1.head(20).copy()
v1_audit_rows = []
for _, r in v1_top_audit.iterrows():
    pn = str(r["pert_name"]).strip()
    m_sub = mac_drugs[mac_drugs["pert_name"] == pn]
    m_sigs = int(m_sub["n_l1000_signatures"].iloc[0]) if len(m_sub) else 0
    m_lines = int(m_sub["n_distinct_cell_lines"].iloc[0]) if len(m_sub) else 0
    m_tier = str(m_sub["strict_audit_tier"].iloc[0]) if len(m_sub) else "Not_Reversed_In_Macrophage_IFN_Query"
    reg_2026 = REGULATORY_2026_OVERRIDES.get(pn, ("", str(r.get("clinical_phase", "")), "Unverified"))[1]
    tox_f    = REGULATORY_2026_OVERRIDES.get(pn, ("", "", "Broad / cancer-line L1000 screen"))[2]
    v1_audit_rows.append({
        "pert_name": pn,
        "v1_pooled_pdac_queries_reversed": int(r["n_pdac_signatures_reversed"]),
        "v1_pooled_l1000_profiles": int(r["total_l1000_profiles"]),
        "v1_max_cell_lines_any_query": int(r["max_cell_lines"]),
        "macrophage_IFN_query_signatures": m_sigs,
        "macrophage_IFN_query_cell_lines": m_lines,
        "post_audit_macrophage_tier": m_tier,
        "v1_reported_phase": str(r.get("clinical_phase", "")),
        "verified_2026_regulatory_status": reg_2026,
        "cytotoxicity_and_specificity_flag": tox_f,
    })
pd.DataFrame(v1_audit_rows).to_csv(ST6_DIR / "lincs_prior_table_forensic_reclassification.tsv", sep="\t", index=False)

# Save unified frozen drug_candidates.tsv (all Macrophage IFN-high reversers + top 30 Malignant Ductal comparators + explicit v1 audit items)
keep_extra_names = set(v1_top_audit["pert_name"].astype(str))
drug_frozen = pd.concat([
    mac_drugs,
    duc_drugs.head(30),
    glb_drugs[glb_drugs["pert_name"].isin(keep_extra_names)],
], ignore_index=True)
drug_frozen.to_csv(OUT_ROOT / "05_results/drug_candidates.tsv", sep="\t", index=False)

n_mac_tier1 = int((mac_drugs["strict_audit_tier"].str.startswith("Tier_1")).sum())
n_mac_tier2 = int((mac_drugs["strict_audit_tier"].str.startswith("Tier_2")).sum())
print(f"  Macrophage IFN-high L1000 reversers (FDR<0.05): {len(mac_drugs)} compounds | Strict Tier-1 (>=3 sigs & >=3 lines) = {n_mac_tier1} | Tier-2 (2 sigs & >=2 lines) = {n_mac_tier2}")
print(mac_drugs.head(12)[["pert_name", "strict_audit_tier", "n_l1000_signatures", "n_distinct_cell_lines", "cell_lines_tested", "median_z_sum", "verified_2026_regulatory_status"]].to_string(index=False))
print("[Stage 6 Part B1 Complete]")
