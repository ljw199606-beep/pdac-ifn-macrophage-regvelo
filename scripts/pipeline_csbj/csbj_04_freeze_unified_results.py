#!/usr/bin/env python3
"""
CSBJ Step 6 Unified Frozen Results Layer (csbj_04_freeze_unified_results.py):
Consolidates all verified Phase 1 & Phase 2 outputs (including 247-subject reconciliation,
97,020->91,496 Visium spot audit, exact edgeR QL LODO, voom donor-bootstrap CIs, N_perm=1,000
degree+expression-matched RegVelo permutations with BH q, and pure-ductal TCGA BH/PH audit)
into a single programmatic JSON + TSV source of truth:
  - CSBJ_Final_Submission/04_Reproducibility/data_manifest/csbj_frozen_numbers.json
  - CSBJ_Final_Submission/06_Quality_Audit/phase2_tables/target_evidence_matrix_csbj.tsv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

BASE_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT   = Path(r"D:\Research\2026_PDAC_RegVelo")
CSBJ_ROOT = BASE_ROOT / "CSBJ_Final_Submission"
MAN_DIR   = CSBJ_ROOT / "04_Reproducibility" / "data_manifest"
P2_DIR    = CSBJ_ROOT / "06_Quality_Audit" / "phase2_tables"
ST3_DIR   = BASE_ROOT / "05_results" / "stage3_macrophage"
ST6_DIR   = BASE_ROOT / "05_results" / "stage6_clinical_drug"

counts_json = json.loads((MAN_DIR / "frozen_cohort_counts.json").read_text(encoding="utf-8"))
comp_stat   = pd.read_csv(BASE_ROOT / "06_figures/source_data/Figure1d_donor_composition_stats_source.tsv", sep="\t").set_index("Cell_type")
q1_tests    = pd.read_csv(ST3_DIR / "between_donor_abundance_tests_Q1.tsv", sep="\t")
de_peng14   = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_Peng_Q2_IFN14_min5.tsv", sep="\t")
de_peng10   = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_Peng_Q2_IFN10_min5.tsv", sep="\t")
de_202_10   = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_GSE202051_All43_Q2_IFN10_min5.tsv", sep="\t")
fgsea_df    = pd.read_csv(ST3_DIR / "noncircular_hallmark_fgsea_all52_pathways.tsv", sep="\t")
edger_lodo  = pd.read_csv(P2_DIR / "edgeR_exact_LODO_and_voom_bootstrap_audit.tsv", sep="\t")
amb_audit   = pd.read_csv(P2_DIR / "ambient_rna_and_background_audit.tsv", sep="\t").set_index("gene")
conc_voom   = pd.read_csv(P2_DIR / "cross_cohort_voom_donor_bootstrap_and_LODO_concordance.tsv", sep="\t")
rgv_1000    = pd.read_csv(P2_DIR / "regvelo_1000perm_macrophage_summary.tsv", sep="\t")
sp_pat      = pd.read_csv(MAN_DIR / "GSE272362_patient_level_morans_summary.tsv", sep="\t")
surv_bh     = pd.read_csv(P2_DIR / "tcga_survival_3models_with_BH_and_PH_audit.tsv", sep="\t")
corr_tcga   = pd.read_csv(ST6_DIR / "tcga_ductal_stat1_partial_correlations_exact.tsv", sep="\t").set_index("gene")
paired_pan  = pd.read_csv(ST3_DIR / "paired_donor_functional_panel_all_cohorts.tsv", sep="\t")

# =========================================================================
# 1. REBUILD COMPARTMENT-SEPARATED EVIDENCE MATRIX WITH UPGRADED 1,000-PERM REGVELO
# =========================================================================
lambert_tfs = set(
    x.strip()
    for x in (V1_ROOT / "02_data/reference/prior_networks/lambert2018_human_tfs.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)

def pos_percentile_score(df_in: pd.DataFrame, key_col: str, stat_col: str, dir_mask: pd.Series, universe_keys: list) -> pd.Series:
    sub = df_in[[key_col, stat_col]].copy().loc[dir_mask.values].dropna(subset=[stat_col]).drop_duplicates(subset=[key_col])
    s_out = pd.Series(0.0, index=universe_keys, dtype=float)
    if len(sub) == 0:
        return s_out
    ranks = sub[stat_col].rank(method="average", pct=True)
    for k, rk in zip(sub[key_col], ranks):
        if k in s_out.index:
            s_out.loc[k] = float(rk)
    return s_out

ulm_nc     = pd.read_csv(ST6_DIR / "peng_macrophage_paired_collectri_ulm.tsv", sep="\t").query("mode == 'NonCircular_NoIFN10_Targets'").copy()
scenic_c14 = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_macrophage_paired_C14_aucell.csv")
scenic_sum = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_motif_regulons_summary.csv")
ot_tract   = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/opentargets_tractability_and_drugs.csv")

mac_tf_universe = sorted([g for g in de_peng10["gene"].unique() if g in lambert_tfs])

p10_tf = de_peng10[de_peng10["gene"].isin(mac_tf_universe) & (~de_peng10["in_defining_signature"])].copy()
p10_tf["neglog10_fdr"] = -np.log10(np.clip(p10_tf["FDR_edgeR_nonDefining"].values, 1e-50, 1.0))
M1 = pos_percentile_score(p10_tf, "gene", "neglog10_fdr", p10_tf["log2FC_edgeR"] > 0, mac_tf_universe)

ulm_tf = ulm_nc[ulm_nc["TF"].isin(mac_tf_universe)].copy()
ulm_tf["neglog10_fdr"] = -np.log10(np.clip(ulm_tf["paired_t_FDR"].values, 1e-50, 1.0))
M2 = pos_percentile_score(ulm_tf, "TF", "neglog10_fdr", ulm_tf["delta_ulm_high_minus_low"] > 0, mac_tf_universe)

sc_tf = scenic_c14[scenic_c14["TF"].isin(mac_tf_universe)].copy()
sc_tf["neglog10_fdr"] = -np.log10(np.clip(sc_tf["wilcoxon_padj"].values, 1e-50, 1.0))
M3 = pos_percentile_score(sc_tf, "TF", "neglog10_fdr", sc_tf["delta_auc_high_vs_low"] > 0, mac_tf_universe)

g202_tf = de_202_10[de_202_10["gene"].isin(mac_tf_universe) & (~de_202_10["in_defining_signature"])].copy()
g202_tf["neglog10_fdr"] = -np.log10(np.clip(g202_tf["FDR_edgeR_nonDefining"].values, 1e-50, 1.0))
M4 = pos_percentile_score(g202_tf, "gene", "neglog10_fdr", g202_tf["log2FC_edgeR"] > 0, mac_tf_universe)

rgv_single = rgv_1000[(rgv_1000["perturbation_type"] == "Single_TF_KO") & (rgv_1000["perturbation"].isin(mac_tf_universe))].copy()
rgv_single["abs_z_null"] = rgv_single["standardized_z_vs_null"].abs()
# Score M5 by empirical percentile of degree/expression-standardized |z_vs_null| (fair across degrees!)
M5 = pos_percentile_score(rgv_single, "perturbation", "abs_z_null", rgv_single["abs_z_null"] > 0, mac_tf_universe)

mac_mat = pd.DataFrame({"TF": mac_tf_universe, "compartment_context": "Macrophage_IFN_State_Primary"})
mac_mat = mac_mat.merge(p10_tf[["gene", "log2FC_edgeR", "FDR_edgeR_nonDefining"]].rename(columns={"gene": "TF", "log2FC_edgeR": "M1_peng_IFN10_log2FC", "FDR_edgeR_nonDefining": "M1_peng_IFN10_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(ulm_tf[["TF", "n_regulon_targets_in_universe", "delta_ulm_high_minus_low", "paired_t_FDR"]].rename(columns={"n_regulon_targets_in_universe": "M2_collectri_n_targets", "delta_ulm_high_minus_low": "M2_collectri_delta_ULM", "paired_t_FDR": "M2_collectri_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(sc_tf[["TF", "delta_auc_high_vs_low", "wilcoxon_padj"]].rename(columns={"delta_auc_high_vs_low": "M3_scenic_delta_AUC", "wilcoxon_padj": "M3_scenic_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(scenic_sum[["TF", "n_targets"]].drop_duplicates("TF").rename(columns={"n_targets": "M3_scenic_motif_targets"}), on="TF", how="left")
mac_mat = mac_mat.merge(g202_tf[["gene", "log2FC_edgeR", "FDR_edgeR_nonDefining"]].rename(columns={"gene": "TF", "log2FC_edgeR": "M4_gse202051_IFN10_log2FC", "FDR_edgeR_nonDefining": "M4_gse202051_IFN10_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(rgv_single[["perturbation", "n_muted_prior_edges", "mean_delta_pp", "cross_seed_sd_pp", "standardized_z_vs_null", "empirical_perm_p", "BH_adjusted_q"]].rename(columns={"perturbation": "TF", "n_muted_prior_edges": "M5_regvelo_muted_edges", "mean_delta_pp": "M5_regvelo_hard_delta_pp", "cross_seed_sd_pp": "M5_regvelo_seed_sd_pp", "standardized_z_vs_null": "M5_regvelo_z_null", "empirical_perm_p": "M5_regvelo_perm1000_p", "BH_adjusted_q": "M5_regvelo_BH_q"}), on="TF", how="left")

mac_mat["score_M1_Peng_IFN10_DEG"] = [M1[t] for t in mac_mat["TF"]]
mac_mat["score_M2_Peng_CollecTRI_ULM"] = [M2[t] for t in mac_mat["TF"]]
mac_mat["score_M3_Peng_pySCENIC_Motif"] = [M3[t] for t in mac_mat["TF"]]
mac_mat["score_M4_GSE202051_Replication"] = [M4[t] for t in mac_mat["TF"]]
mac_mat["score_M5_GSE154778_RegVelo_1000Perm"] = [M5[t] for t in mac_mat["TF"]]

S_mat = mac_mat[["score_M1_Peng_IFN10_DEG", "score_M2_Peng_CollecTRI_ULM", "score_M3_Peng_pySCENIC_Motif", "score_M4_GSE202051_Replication", "score_M5_GSE154778_RegVelo_1000Perm"]].values
peng_block = S_mat[:, 0:3].mean(axis=1)
mac_mat["cohort_balanced_score"] = (peng_block + S_mat[:, 3] + S_mat[:, 4]) / 3.0
mac_mat["equal_5layer_score"]    = S_mat.mean(axis=1)
mac_mat["primary_rank"]          = mac_mat["cohort_balanced_score"].rank(ascending=False, method="min").astype(int)
mac_mat["equal_5layer_rank"]     = mac_mat["equal_5layer_score"].rank(ascending=False, method="min").astype(int)

for drop_idx, col_nm in [(0, "rank_drop_M1"), (1, "rank_drop_M2"), (2, "rank_drop_M3"), (3, "rank_drop_M4_GSE202051"), (4, "rank_drop_M5_RegVelo")]:
    if drop_idx in (0, 1, 2):
        rem = [i for i in (0, 1, 2) if i != drop_idx]
        sc_tmp = (S_mat[:, rem].mean(axis=1) + S_mat[:, 3] + S_mat[:, 4]) / 3.0
    elif drop_idx == 3:
        sc_tmp = (peng_block + S_mat[:, 4]) / 2.0
    else:
        sc_tmp = (peng_block + S_mat[:, 3]) / 2.0
    mac_mat[col_nm] = pd.Series(sc_tmp).rank(ascending=False, method="min").astype(int).values

rng = np.random.default_rng(20261118)
W_dir = rng.dirichlet([1.0, 1.0, 1.0], size=1000)
C3 = np.column_stack([peng_block, S_mat[:, 3], S_mat[:, 4]])
scores_dir = C3 @ W_dir.T
ranks_dir = np.apply_along_axis(lambda x: pd.Series(x).rank(ascending=False, method="min").values, 0, scores_dir)
mac_mat["dirichlet_median_rank"] = np.median(ranks_dir, axis=1)
mac_mat["dirichlet_iqr_rank"]    = np.percentile(ranks_dir, 75, axis=1) - np.percentile(ranks_dir, 25, axis=1)
mac_mat["dirichlet_top5_freq"]   = np.mean(ranks_dir <= 5, axis=1)

mac_mat = mac_mat.merge(
    ot_tract[["gene", "druggability_tier", "max_clinical_phase", "sm_tractable", "protac_tractable"]].rename(columns={"gene": "TF"}),
    on="TF", how="left",
)
mac_mat = mac_mat.sort_values(["primary_rank", "equal_5layer_rank"]).reset_index(drop=True)
mac_mat.to_csv(P2_DIR / "macrophage_ifn_tf_evidence_matrix_csbj_full.tsv", sep="\t", index=False)

epi_prev = pd.read_csv(ST6_DIR / "epithelial_malignant_ductal_tf_reference_matrix.tsv", sep="\t")
pd.concat([mac_mat.head(80), epi_prev.head(40)], ignore_index=True, sort=False).to_csv(P2_DIR / "target_evidence_matrix_csbj.tsv", sep="\t", index=False)

print("Top 10 Macrophage IFN-State Regulators (CSBJ 1,000-Perm Upgraded Matrix):")
print(mac_mat.head(10)[["primary_rank", "TF", "cohort_balanced_score", "rank_drop_M5_RegVelo", "dirichlet_median_rank", "M5_regvelo_hard_delta_pp", "M5_regvelo_perm1000_p", "M5_regvelo_BH_q"]].to_string(index=False))

# =========================================================================
# 2. ASSEMBLE MASTER csbj_frozen_numbers.json FOR ZERO-TRANSCRIPTION MANUSCRIPT IMPORT
# =========================================================================
q1_p14 = q1_tests[(q1_tests["cohort"] == "Peng2019") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == "ifn14_high_frac_of_mac_Q1")].iloc[0]
q1_p10 = q1_tests[(q1_tests["cohort"] == "Peng2019") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == "ifn10_high_frac_of_mac_Q1")].iloc[0]
q1_pall = q1_tests[(q1_tests["cohort"] == "Peng2019") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == "ifn14_high_frac_of_all_cells_Q1")].iloc[0]
q1_p14_m10 = q1_tests[(q1_tests["cohort"] == "Peng2019") & (q1_tests["min_macrophages_per_donor"] == 10) & (q1_tests["metric"] == "ifn14_high_frac_of_mac_Q1")].iloc[0]
q1_202 = q1_tests[(q1_tests["cohort"] == "GSE202051") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == "ifn14_high_frac_of_mac_Q1")].iloc[0]
q1_154 = q1_tests[(q1_tests["cohort"] == "GSE154778") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == "ifn14_high_frac_of_mac_Q1")].iloc[0]

fg_col = "tag" if "tag" in fgsea_df.columns else "model"
fg_peng = fgsea_df[fgsea_df[fg_col] == "Peng_Q2_IFN14_min5"].set_index("pathway")
lodo_p14 = edger_lodo[edger_lodo["cohort_model"] == "Peng2019_Q2_IFN14_min5"].set_index("gene")
lodo_p10 = edger_lodo[edger_lodo["cohort_model"] == "Peng2019_Q2_IFN10_min5"].set_index("gene")
surv_d   = surv_bh[(surv_bh["cohort"] == "Pure_Ductal_PDAC_n146") & (surv_bh["endpoint"] == "OS")].set_index("feature")

def get_conc(vt: str, gs_prefix: str) -> dict:
    r = conc_voom[(conc_voom["validation_cohort_model"] == vt) & (conc_voom["gene_set_definition"].str.startswith(gs_prefix))].iloc[0]
    return {
        "n_genes": int(r["n_genes"]),
        "edgeR_rho": round(float(r["edgeR_spearman_rho"]), 3),
        "voom_rho": round(float(r["voom_precision_spearman_rho"]), 3),
        "boot_ci_low": round(float(r["donor_bootstrap_95CI_low"]), 3),
        "boot_ci_high": round(float(r["donor_bootstrap_95CI_high"]), 3),
        "lodo_min": round(float(r["LODO_min_spearman_rho"]), 3),
        "lodo_max": round(float(r["LODO_max_spearman_rho"]), 3),
        "sign_pct": round(float(r["sign_agreement_pct"]), 1),
    }

rgv_s_sorted = rgv_1000[rgv_1000["perturbation_type"] == "Single_TF_KO"].sort_values(["empirical_perm_p", "abs_mean_delta_pp"], ascending=[True, False]).reset_index(drop=True)
rgv_by_name  = rgv_1000.set_index("perturbation")

def get_rgv_tf(nm: str) -> dict:
    r = rgv_by_name.loc[nm]
    return {
        "edges": int(r["n_muted_prior_edges"]),
        "delta_pp": round(float(r["mean_delta_pp"]), 2),
        "seed_sd_pp": round(float(r["cross_seed_sd_pp"]), 2),
        "donor_ci_low_pp": round(float(r["donor16_bootstrap_95CI_low_pp"]), 2),
        "donor_ci_high_pp": round(float(r["donor16_bootstrap_95CI_high_pp"]), 2),
        "z_null": round(float(r["standardized_z_vs_null"]), 2),
        "perm_p": round(float(r["empirical_perm_p"]), 4),
        "bh_q": round(float(r["BH_adjusted_q"]), 4) if pd.notna(r["BH_adjusted_q"]) else None,
    }

sp_pt = sp_pat[sp_pat["anatomical_group"] == "Primary_Tumor"].set_index("Feature Name")
sp_hm = sp_pat[sp_pat["anatomical_group"] == "Hepatic_Metastasis"].set_index("Feature Name")
sp_lnm = sp_pat[sp_pat["anatomical_group"] == "Lymph_Node_Metastasis"].set_index("Feature Name")

counts_flat = dict(counts_json)
counts_flat.update({
    "primary_5cohorts_unique_subjects": int(counts_json["n_unique_subjects_primary_247"]),
    "primary_5cohorts_pdac_patients": int(counts_json["n_pdac_patients_primary_236"]),
    "primary_5cohorts_control_donors": int(counts_json["n_control_donors_11"]),
    "primary_5cohorts_specimens_or_sections": int(counts_json["n_specimens_primary_270"]),
    "cra001160_qc_cells": int(counts_json["Peng2019"]["n_cells_qc"]),
    "cra001160_tams": int(counts_json["Peng2019"]["n_mac"]),
    "gse202051_tams": int(counts_json["GSE202051"]["n_mac"]),
    "gse154778_srr_runs": int(counts_json["GSE154778"]["n_srr_runs"]),
    "gse154778_starsolo_raw_cells": int(counts_json["GSE154778"]["n_emptydrops"]),
    "gse154778_post_qc_all_cells": int(counts_json["GSE154778"]["n_cells_qc"]),
    "gse154778_tams": int(counts_json["GSE154778"]["n_mac"]),
    "gse154778_malignant_ductal_subset": int(counts_json["GSE154778"]["n_ductal_malignant_subset"]),
    "gse272362_array_barcodes": int(counts_json["GSE272362"]["n_array_barcodes"]),
    "gse272362_in_tissue_spots_30_sections": int(counts_json["GSE272362"]["n_in_tissue_spots_spaceranger"]),
    "gse272362_in_tissue_spots_29_unique_lesions": int(counts_json["GSE272362"]["n_in_tissue_spots_29sec"]),
    "gse272362_khaliq_post_qc_spots": int(counts_json["GSE272362"]["n_post_qc_spots_khaliq_series"]),
    "tcga_pure_ductal_patients": int(counts_json["TCGA_PAAD"]["n_pure_ductal"]),
    "tcga_all_primary_paad_sensitivity": int(counts_json["TCGA_PAAD"]["n_all_paad"]),
})

frozen = {
    "target_journal": "Computational and Structural Biotechnology Journal (CSBJ, Science Partner Journals)",
    "article_type": "Research Article (General Section)",
    "counts": counts_flat,
    "peng_composition": {
        "ductal2_med_pdac_pct": round(float(comp_stat.loc["Ductal cell type 2", "median_pct_tumor"]), 1),
        "ductal2_med_ctrl_pct": round(float(comp_stat.loc["Ductal cell type 2", "median_pct_normal"]), 1),
        "ductal2_bh_q": float(comp_stat.loc["Ductal cell type 2", "mwu_padj"]),
        "mac_med_pdac_pct": round(float(comp_stat.loc["Macrophage cell", "median_pct_tumor"]), 1),
        "mac_med_ctrl_pct": round(float(comp_stat.loc["Macrophage cell", "median_pct_normal"]), 1),
        "mac_fold_mean": round(float(comp_stat.loc["Macrophage cell", "fold_change_mean"]), 2),
        "mac_mwu_p": float(comp_stat.loc["Macrophage cell", "mwu_p"]),
        "mac_bh_q": float(comp_stat.loc["Macrophage cell", "mwu_padj"]),
        "tcell_med_pdac_pct": round(float(comp_stat.loc["T cell", "median_pct_tumor"]), 1),
        "tcell_med_ctrl_pct": round(float(comp_stat.loc["T cell", "median_pct_normal"]), 2),
        "tcell_bh_q": float(comp_stat.loc["T cell", "mwu_padj"]),
        "stellate_med_pdac_pct": round(float(comp_stat.loc["Stellate cell", "median_pct_tumor"]), 1),
        "stellate_med_ctrl_pct": round(float(comp_stat.loc["Stellate cell", "median_pct_normal"]), 1),
        "stellate_bh_q": float(comp_stat.loc["Stellate cell", "mwu_padj"]),
    },
    "q1_abundance": {
        "peng_ifn14_med_pdac_pct": round(float(q1_p14["median_group1"]) * 100, 1),
        "peng_ifn14_med_ctrl_pct": round(float(q1_p14["median_group2"]) * 100, 1),
        "peng_ifn14_hl_pct": round(float(q1_p14["hodges_lehmann_diff"]) * 100, 1),
        "peng_ifn14_ci_low_pct": round(float(q1_p14["ci95_low"]) * 100, 1),
        "peng_ifn14_ci_high_pct": round(float(q1_p14["ci95_high"]) * 100, 1),
        "peng_ifn14_mwu_p": float(q1_p14["mwu_p"]),
        "peng_ifn14_min10_mwu_p": float(q1_p14_m10["mwu_p"]),
        "peng_ifn10_med_pdac_pct": round(float(q1_p10["median_group1"]) * 100, 1),
        "peng_ifn10_med_ctrl_pct": round(float(q1_p10["median_group2"]) * 100, 1),
        "peng_ifn10_mwu_p": float(q1_p10["mwu_p"]),
        "peng_allcells_med_pdac_pct": round(float(q1_pall["median_group1"]) * 100, 2),
        "peng_allcells_med_ctrl_pct": round(float(q1_pall["median_group2"]) * 100, 2),
        "peng_allcells_mwu_p": float(q1_pall["mwu_p"]),
        "gse202051_neoadj_med_pct": round(float(q1_202["median_group1"]) * 100, 1),
        "gse202051_untr_med_pct": round(float(q1_202["median_group2"]) * 100, 1),
        "gse202051_hl_pct": round(float(q1_202["hodges_lehmann_diff"]) * 100, 1),
        "gse202051_ci_low_pct": round(float(q1_202["ci95_low"]) * 100, 1),
        "gse202051_ci_high_pct": round(float(q1_202["ci95_high"]) * 100, 1),
        "gse202051_mwu_p": round(float(q1_202["mwu_p"]), 3),
        "gse154778_met_med_pct": round(float(q1_154["median_group1"]) * 100, 1),
        "gse154778_pt_med_pct": round(float(q1_154["median_group2"]) * 100, 1),
        "gse154778_hl_pct": round(float(q1_154["hodges_lehmann_diff"]) * 100, 1),
        "gse154778_ci_low_pct": round(float(q1_154["ci95_low"]) * 100, 1),
        "gse154778_ci_high_pct": round(float(q1_154["ci95_high"]) * 100, 1),
        "gse154778_mwu_p": round(float(q1_154["mwu_p"]), 4),
        "gse154778_mwu_fdr": round(float(q1_154["mwu_fdr_within_contrast"]), 3),
    },
    "paired_pseudobulk_and_lodo": {
        "n_tested_genes_peng": int(len(de_peng14)),
        "n_nondefining_DEGs_peng_FDR005": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05)).sum()),
        "n_nondefining_DEGs_peng_FDR005_up": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05) & (de_peng14["log2FC_edgeR"] > 0)).sum()),
        "n_nondefining_DEGs_peng_FDR005_down": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05) & (de_peng14["log2FC_edgeR"] < 0)).sum()),
        "n_nondefining_DEGs_peng_FDR005_absFC025": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05) & (de_peng14["log2FC_edgeR"].abs() > 0.25)).sum()),
        "n_nondefining_DEGs_peng_FDR005_up_FC025": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05) & (de_peng14["log2FC_edgeR"] > 0.25)).sum()),
        "n_nondefining_DEGs_peng_FDR005_down_FC025": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05) & (de_peng14["log2FC_edgeR"] < -0.25)).sum()),
        "STAT1_log2FC": round(float(lodo_p14.loc["STAT1", "full_edgeR_log2FC"]), 3),
        "STAT1_ci_low": round(float(lodo_p14.loc["STAT1", "full_edgeR_ci95_low"]), 2),
        "STAT1_ci_high": round(float(lodo_p14.loc["STAT1", "full_edgeR_ci95_high"]), 2),
        "STAT1_FDR": float(lodo_p14.loc["STAT1", "full_edgeR_nonDefining_FDR"]),
        "STAT1_LODO_FDR_runs": str(lodo_p14.loc["STAT1", "edgeR_LODO_runs_FDR_lt_005"]),
        "IRF1_IFN10_log2FC": round(float(lodo_p10.loc["IRF1", "full_edgeR_log2FC"]), 3),
        "IRF1_IFN10_FDR": float(lodo_p10.loc["IRF1", "full_edgeR_nonDefining_FDR"]),
        "IRF1_IFN10_LODO_FDR_runs": str(lodo_p10.loc["IRF1", "edgeR_LODO_runs_FDR_lt_005"]),
        "CD47_log2FC": round(float(lodo_p14.loc["CD47", "full_edgeR_log2FC"]), 3),
        "CD47_ci_low": round(float(lodo_p14.loc["CD47", "full_edgeR_ci95_low"]), 2),
        "CD47_ci_high": round(float(lodo_p14.loc["CD47", "full_edgeR_ci95_high"]), 2),
        "CD47_FDR": float(lodo_p14.loc["CD47", "full_edgeR_nonDefining_FDR"]),
        "CD47_voom_boot_low": round(float(lodo_p14.loc["CD47", "voom_donor_bootstrap_95CI_low"]), 3),
        "CD47_voom_boot_high": round(float(lodo_p14.loc["CD47", "voom_donor_bootstrap_95CI_high"]), 3),
        "CD47_LODO_p_runs": str(lodo_p14.loc["CD47", "edgeR_LODO_runs_p_lt_005"]),
        "CD47_LODO_FDR_runs": str(lodo_p14.loc["CD47", "edgeR_LODO_runs_FDR_lt_005"]),
        "IRF7_log2FC": round(float(lodo_p14.loc["IRF7", "full_edgeR_log2FC"]), 3),
        "IRF7_FDR": float(lodo_p14.loc["IRF7", "full_edgeR_nonDefining_FDR"]),
        "STAT2_log2FC": round(float(lodo_p14.loc["STAT2", "full_edgeR_log2FC"]), 3),
        "STAT2_FDR": float(lodo_p14.loc["STAT2", "full_edgeR_nonDefining_FDR"]),
        "TAP1_log2FC": round(float(lodo_p14.loc["TAP1", "full_edgeR_log2FC"]), 3),
        "TAP1_FDR": float(lodo_p14.loc["TAP1", "full_edgeR_nonDefining_FDR"]),
        "HLAA_log2FC": round(float(lodo_p14.loc["HLA-A", "full_edgeR_log2FC"]), 3),
        "HLAA_FDR": float(lodo_p14.loc["HLA-A", "full_edgeR_nonDefining_FDR"]),
        "B2M_log2FC": round(float(lodo_p14.loc["B2M", "full_edgeR_log2FC"]), 3),
        "B2M_FDR": float(lodo_p14.loc["B2M", "full_edgeR_nonDefining_FDR"]),
        "CD80_log2FC": round(float(lodo_p14.loc["CD80", "full_edgeR_log2FC"]), 3),
        "CD80_PValue": float(lodo_p14.loc["CD80", "full_edgeR_PValue"]),
        "CD80_FDR": round(float(lodo_p14.loc["CD80", "full_edgeR_nonDefining_FDR"]), 3),
        "CD80_LODO_p_runs": str(lodo_p14.loc["CD80", "edgeR_LODO_runs_p_lt_005"]),
        "LGALS9_log2FC": round(float(lodo_p14.loc["LGALS9", "full_edgeR_log2FC"]), 3),
        "LGALS9_PValue": float(lodo_p14.loc["LGALS9", "full_edgeR_PValue"]),
        "LGALS9_FDR": round(float(lodo_p14.loc["LGALS9", "full_edgeR_nonDefining_FDR"]), 3),
        "LGALS9_LODO_p_runs": str(lodo_p14.loc["LGALS9", "edgeR_LODO_runs_p_lt_005"]),
        "COL1A1_unpaired_FDR": float(amb_audit.loc["COL1A1", "unpaired_Tumor_vs_Control_TAM_FDR"]),
        "COL1A1_paired_log2FC": round(float(lodo_p14.loc["COL1A1", "full_edgeR_log2FC"]), 3),
        "COL1A1_paired_FDR": round(float(lodo_p14.loc["COL1A1", "full_edgeR_nonDefining_FDR"]), 3),
        "KRT19_unpaired_FDR": float(amb_audit.loc["KRT19", "unpaired_Tumor_vs_Control_TAM_FDR"]),
        "KRT19_paired_log2FC": round(float(lodo_p14.loc["KRT19", "full_edgeR_log2FC"]), 3),
        "KRT19_paired_FDR": round(float(lodo_p14.loc["KRT19", "full_edgeR_nonDefining_FDR"]), 3),
        "fgsea_MHCI_NES": round(float(fg_peng.loc["CUSTOM_ANTIGEN_PRESENTATION_MHC_NO_IFN14", "NES"]), 4),
        "fgsea_MHCI_padj": float(fg_peng.loc["CUSTOM_ANTIGEN_PRESENTATION_MHC_NO_IFN14", "padj"]),
        "fgsea_Checkpoint_NES": round(float(fg_peng.loc["CUSTOM_ADAPTIVE_CHECKPOINT_LIGANDS_NO_IFN14", "NES"]), 4),
        "fgsea_Checkpoint_padj": float(fg_peng.loc["CUSTOM_ADAPTIVE_CHECKPOINT_LIGANDS_NO_IFN14", "padj"]),
    },
    "cross_cohort_concordance": {
        "GSE202051_All43_AllGenes": get_conc("GSE202051_All43_Q2_IFN14_min5", "1_"),
        "GSE202051_All43_Panel23": get_conc("GSE202051_All43_Q2_IFN14_min5", "2_"),
        "GSE202051_All43_DiscDEGs": get_conc("GSE202051_All43_Q2_IFN14_min5", "3_"),
        "GSE202051_Untreated18_DiscDEGs": get_conc("GSE202051_Untreated18_Q2_IFN14_min5", "3_"),
        "GSE202051_Neoadjuvant25_DiscDEGs": get_conc("GSE202051_Neoadjuvant25_Q2_IFN14_min5", "3_"),
        "GSE154778_All16_AllGenes": get_conc("GSE154778_All16_Q2_IFN14_min5", "1_"),
        "GSE154778_All16_Panel23": get_conc("GSE154778_All16_Q2_IFN14_min5", "2_"),
        "GSE154778_All16_DiscDEGs": get_conc("GSE154778_All16_Q2_IFN14_min5", "3_"),
    },
    "regvelo_1000perm": {
        "n_permutations": 1000,
        "n_seeds": 3,
        "n_single_tfs_tested": int(len(rgv_s_sorted)),
        "n_emp_p_lt_001": int((rgv_s_sorted["empirical_perm_p"] < 0.01).sum()),
        "n_emp_p_lt_005": int((rgv_s_sorted["empirical_perm_p"] < 0.05).sum()),
        "n_bh_q_lt_005": int((rgv_s_sorted["BH_adjusted_q"] < 0.05).sum()),
        "n_bh_q_lt_010": int((rgv_s_sorted["BH_adjusted_q"] < 0.10).sum()),
        "top5_single_tfs": rgv_s_sorted.head(5)["perturbation"].tolist(),
        "NFKB1": get_rgv_tf("NFKB1"),
        "REL": get_rgv_tf("REL"),
        "IRF4": get_rgv_tf("IRF4"),
        "NR4A3": get_rgv_tf("NR4A3"),
        "RUNX1": get_rgv_tf("RUNX1"),
        "FOS": get_rgv_tf("FOS"),
        "STAT1": get_rgv_tf("STAT1"),
        "IRF1": get_rgv_tf("IRF1"),
        "STAT2": get_rgv_tf("STAT2"),
        "IRF2": get_rgv_tf("IRF2"),
        "STAT1_plus_IRF1": get_rgv_tf("STAT1+IRF1"),
        "GAPDH": get_rgv_tf("GAPDH"),
    },
    "evidence_matrix": {
        "top5_macrophage_tfs": mac_mat.head(5)["TF"].tolist(),
        "rank1_tf": str(mac_mat.iloc[0]["TF"]),
        "rank1_score": round(float(mac_mat.iloc[0]["cohort_balanced_score"]), 3),
        "rank2_tf": str(mac_mat.iloc[1]["TF"]),
        "rank2_score": round(float(mac_mat.iloc[1]["cohort_balanced_score"]), 3),
        "STAT1_rank": int(mac_mat.loc[mac_mat["TF"] == "STAT1", "primary_rank"].iloc[0]),
        "STAT1_score": round(float(mac_mat.loc[mac_mat["TF"] == "STAT1", "cohort_balanced_score"].iloc[0]), 3),
        "STAT1_drop_M5_rank": int(mac_mat.loc[mac_mat["TF"] == "STAT1", "rank_drop_M5_RegVelo"].iloc[0]),
        "IRF1_rank": int(mac_mat.loc[mac_mat["TF"] == "IRF1", "primary_rank"].iloc[0]),
        "IRF1_score": round(float(mac_mat.loc[mac_mat["TF"] == "IRF1", "cohort_balanced_score"].iloc[0]), 3),
        "IRF2_rank": int(mac_mat.loc[mac_mat["TF"] == "IRF2", "primary_rank"].iloc[0]),
        "STAT2_rank": int(mac_mat.loc[mac_mat["TF"] == "STAT2", "primary_rank"].iloc[0]),
    },
    "spatial_gse272362": {
        "STAT1_PT_pat_med_I": round(float(sp_pt.loc["STAT1", "patient_level_median_I"]), 3),
        "STAT1_PT_sec30_med_I": round(float(sp_pt.loc["STAT1", "section30_median_I"]), 3),
        "STAT1_PT_sec29_med_I": round(float(sp_pt.loc["STAT1", "section29_median_I"]), 3),
        "STAT1_PT_n_pat_sig": f"{int(sp_pt.loc['STAT1', 'n_patients_with_FDR005_section'])}/{int(sp_pt.loc['STAT1', 'n_unique_patients'])}",
        "STAT1_HM_pat_med_I": round(float(sp_hm.loc["STAT1", "patient_level_median_I"]), 3),
        "STAT1_LNM_pat_med_I": round(float(sp_lnm.loc["STAT1", "patient_level_median_I"]), 3),
        "B2M_PT_pat_med_I": round(float(sp_pt.loc["B2M", "patient_level_median_I"]), 3),
        "CD47_PT_pat_med_I": round(float(sp_pt.loc["CD47", "patient_level_median_I"]), 3),
    },
    "tcga_clinical": {
        "STAT1_HR_multi": round(float(surv_d.loc["STAT1", "HR_multi"]), 2),
        "STAT1_ci_low": round(float(surv_d.loc["STAT1", "CI95_low_multi"]), 2),
        "STAT1_ci_high": round(float(surv_d.loc["STAT1", "CI95_high_multi"]), 2),
        "STAT1_p_multi": round(float(surv_d.loc["STAT1", "p_multi"]), 3),
        "STAT1_bh_q_multi": round(float(surv_d.loc["STAT1", "BH_q_multi"]), 3),
        "IRF1_HR_multi": round(float(surv_d.loc["IRF1", "HR_multi"]), 2),
        "IRF1_p_multi": round(float(surv_d.loc["IRF1", "p_multi"]), 3),
        "FOXM1_HR_multi": round(float(surv_d.loc["FOXM1", "HR_multi"]), 2),
        "FOXM1_ci_low": round(float(surv_d.loc["FOXM1", "CI95_low_multi"]), 2),
        "FOXM1_ci_high": round(float(surv_d.loc["FOXM1", "CI95_high_multi"]), 2),
        "FOXM1_p_multi": round(float(surv_d.loc["FOXM1", "p_multi"]), 4),
        "FOXM1_bh_q_multi": round(float(surv_d.loc["FOXM1", "BH_q_multi"]), 3),
        "HMGA2_HR_multi": round(float(surv_d.loc["HMGA2", "HR_multi"]), 2),
        "HMGA2_ci_low": round(float(surv_d.loc["HMGA2", "CI95_low_multi"]), 2),
        "HMGA2_ci_high": round(float(surv_d.loc["HMGA2", "CI95_high_multi"]), 2),
        "HMGA2_p_multi": round(float(surv_d.loc["HMGA2", "p_multi"]), 4),
        "GATA6_HR_multi": round(float(surv_d.loc["GATA6", "HR_multi"]), 2),
        "GATA6_ci_low": round(float(surv_d.loc["GATA6", "CI95_low_multi"]), 2),
        "GATA6_ci_high": round(float(surv_d.loc["GATA6", "CI95_high_multi"]), 2),
        "GATA6_p_multi": round(float(surv_d.loc["GATA6", "p_multi"]), 3),
        "partial_rho_CXCL10": round(float(corr_tcga.loc["CXCL10", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_rho_CD274": round(float(corr_tcga.loc["CD274", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_rho_IDO1": round(float(corr_tcga.loc["IDO1", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_rho_CD47": round(float(corr_tcga.loc["CD47", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_rho_PDCD1": round(float(corr_tcga.loc["PDCD1", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_p_PDCD1": round(float(corr_tcga.loc["PDCD1", "partial_p_adj_PTPRC_Mac_Purity"]), 3),
        "partial_rho_VSIR": round(float(corr_tcga.loc["VSIR", "partial_rho_adj_PTPRC_Mac_Purity"]), 2),
        "partial_p_VSIR": round(float(corr_tcga.loc["VSIR", "partial_p_adj_PTPRC_Mac_Purity"]), 3),
    },
}
(MAN_DIR / "csbj_frozen_numbers.json").write_text(json.dumps(frozen, indent=2), encoding="utf-8")
print("[CSBJ Step 6 Unified Frozen Numbers Saved to csbj_frozen_numbers.json!]")
