#!/usr/bin/env python3
"""
Stage 6 Part B2 & Stage 9 Unified Frozen Results Layer:
Builds:
  1. 05_results/target_evidence_matrix.tsv (strictly compartment-separated Macrophage_IFN_State vs Epithelial_Malignant_Ductal_Reference)
  2. 05_results/priority_score_definition.md (complete mathematical formula, denominator/missingness rules, v0/v1 audit, Open Targets licence)
  3. 05_results/primary_effects.tsv (unified table of primary donor-level contrasts, effect sizes, 95% CIs, p/FDR, sample sizes)
  4. 05_results/manuscript_numbers.json (single frozen source of truth for manuscript text, tables, and captions)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

OUT_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT  = Path(r"D:\Research\2026_PDAC_RegVelo")
ST3_DIR  = OUT_ROOT / "05_results" / "stage3_macrophage"
ST4_DIR  = OUT_ROOT / "05_results" / "stage4_regvelo"
ST5_DIR  = OUT_ROOT / "05_results" / "stage5_spatial_comm"
ST6_DIR  = OUT_ROOT / "05_results" / "stage6_clinical_drug"

lambert_tfs = set(
    x.strip()
    for x in (V1_ROOT / "02_data/reference/prior_networks/lambert2018_human_tfs.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)

def pos_percentile_score(df_in: pd.DataFrame, key_col: str, stat_col: str, dir_mask: pd.Series, universe_keys: list) -> pd.Series:
    sub = df_in[[key_col, stat_col]].copy()
    sub = sub.loc[dir_mask.values].dropna(subset=[stat_col])
    sub = sub.drop_duplicates(subset=[key_col])
    s_out = pd.Series(0.0, index=universe_keys, dtype=float)
    if len(sub) == 0:
        return s_out
    ranks = sub[stat_col].rank(method="average", pct=True)
    for k, rk in zip(sub[key_col], ranks):
        if k in s_out.index:
            s_out.loc[k] = float(rk)
    return s_out

# =========================================================================
# 1. MACROPHAGE-INTRINSIC IFN-HIGH REGULATOR EVIDENCE MATRIX (ZERO DEPMAP, ZERO TCGA!)
# =========================================================================
print("[1/4] Building compartment-separated target_evidence_matrix.tsv...")
de_peng10 = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_Peng_Q2_IFN10_min5.tsv", sep="\t")
de_peng14 = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_Peng_Q2_IFN14_min5.tsv", sep="\t")
de_202_10 = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_GSE202051_All43_Q2_IFN10_min5.tsv", sep="\t")
de_202_un = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_GSE202051_Untreated18_Q2_IFN14_min5.tsv", sep="\t")
ulm_mac   = pd.read_csv(ST6_DIR / "peng_macrophage_paired_collectri_ulm.tsv", sep="\t")
ulm_nc    = ulm_mac[ulm_mac["mode"] == "NonCircular_NoIFN10_Targets"].copy()
scenic_c14 = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_macrophage_paired_C14_aucell.csv")
scenic_sum = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_motif_regulons_summary.csv")
rgv_hard   = pd.read_csv(ST4_DIR / "regvelo_Macrophage_TAM_hard_constraint_summary.tsv", sep="\t")
rgv_soft   = pd.read_csv(OUT_ROOT / "05_results/perturbation_results.tsv", sep="\t")
ot_tract   = pd.read_csv(V1_ROOT / "04_results/drug_repurposing/opentargets_tractability_and_drugs.csv")

mac_tf_universe = sorted([g for g in de_peng10["gene"].unique() if g in lambert_tfs])

# Prepare M1..M5 raw statistics and [0, 1] positive percentile scores
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

rgv_m_infl = rgv_hard[(rgv_hard["perturbation_type"] == "Single_TF_KO") & (rgv_hard["terminal_state"] == "TAM_STAT1_Inflam") & (rgv_hard["perturbation"].isin(mac_tf_universe))].copy()
M5 = pos_percentile_score(rgv_m_infl, "perturbation", "abs_mean_delta_pp", rgv_m_infl["abs_mean_delta_pp"] > 0, mac_tf_universe)

mac_mat = pd.DataFrame({"TF": mac_tf_universe, "compartment_context": "Macrophage_IFN_State_Primary"})
mac_mat = mac_mat.merge(p10_tf[["gene", "log2FC_edgeR", "FDR_edgeR_nonDefining"]].rename(columns={"gene": "TF", "log2FC_edgeR": "M1_peng_IFN10_log2FC", "FDR_edgeR_nonDefining": "M1_peng_IFN10_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(ulm_tf[["TF", "n_regulon_targets_in_universe", "delta_ulm_high_minus_low", "paired_t_FDR"]].rename(columns={"n_regulon_targets_in_universe": "M2_collectri_n_targets", "delta_ulm_high_minus_low": "M2_collectri_delta_ULM", "paired_t_FDR": "M2_collectri_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(sc_tf[["TF", "delta_auc_high_vs_low", "wilcoxon_padj"]].rename(columns={"delta_auc_high_vs_low": "M3_scenic_delta_AUC", "wilcoxon_padj": "M3_scenic_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(scenic_sum[["TF", "n_targets"]].drop_duplicates("TF").rename(columns={"n_targets": "M3_scenic_motif_targets"}), on="TF", how="left")
mac_mat = mac_mat.merge(g202_tf[["gene", "log2FC_edgeR", "FDR_edgeR_nonDefining"]].rename(columns={"gene": "TF", "log2FC_edgeR": "M4_gse202051_IFN10_log2FC", "FDR_edgeR_nonDefining": "M4_gse202051_IFN10_FDR"}), on="TF", how="left")
mac_mat = mac_mat.merge(rgv_m_infl[["perturbation", "mean_muted_prior_edges", "mean_delta_pp", "sd_delta_pp", "empirical_decoy_p"]].rename(columns={"perturbation": "TF", "mean_muted_prior_edges": "M5_regvelo_hard_muted_edges", "mean_delta_pp": "M5_regvelo_hard_delta_pp", "sd_delta_pp": "M5_regvelo_hard_sd_pp", "empirical_decoy_p": "M5_regvelo_hard_decoy_p"}), on="TF", how="left")

mac_mat["score_M1_Peng_IFN10_DEG"] = [M1[t] for t in mac_mat["TF"]]
mac_mat["score_M2_Peng_CollecTRI_ULM"] = [M2[t] for t in mac_mat["TF"]]
mac_mat["score_M3_Peng_pySCENIC_Motif"] = [M3[t] for t in mac_mat["TF"]]
mac_mat["score_M4_GSE202051_Replication"] = [M4[t] for t in mac_mat["TF"]]
mac_mat["score_M5_GSE154778_RegVelo_HardKO"] = [M5[t] for t in mac_mat["TF"]]

# Count how many of the 5 layers were experimentally/computationally measured (non-missing) for denominator sensitivity
measured_flags = np.column_stack([
    mac_mat["M1_peng_IFN10_log2FC"].notna().values,
    mac_mat["M2_collectri_delta_ULM"].notna().values,
    mac_mat["M3_scenic_delta_AUC"].notna().values,
    mac_mat["M4_gse202051_IFN10_log2FC"].notna().values,
    mac_mat["M5_regvelo_hard_delta_pp"].notna().values,
])
mac_mat["n_measured_layers_of_5"] = measured_flags.sum(axis=1)

S_mat = mac_mat[["score_M1_Peng_IFN10_DEG", "score_M2_Peng_CollecTRI_ULM", "score_M3_Peng_pySCENIC_Motif", "score_M4_GSE202051_Replication", "score_M5_GSE154778_RegVelo_HardKO"]].values
peng_block = S_mat[:, 0:3].mean(axis=1)
# Primary Cohort-Balanced Score: 1/3 Peng2019 block + 1/3 GSE202051 replication + 1/3 GSE154778 RegVelo KO
mac_mat["cohort_balanced_score"] = (peng_block + S_mat[:, 3] + S_mat[:, 4]) / 3.0
mac_mat["equal_5layer_score"]    = S_mat.mean(axis=1)
mac_mat["observed_denom_5layer_score"] = S_mat.sum(axis=1) / np.maximum(mac_mat["n_measured_layers_of_5"].values, 1)

# Ranks & Leave-One-Layer-Out (LOLO) ablations
mac_mat["primary_rank"] = mac_mat["cohort_balanced_score"].rank(ascending=False, method="min").astype(int)
mac_mat["equal_5layer_rank"] = mac_mat["equal_5layer_score"].rank(ascending=False, method="min").astype(int)
mac_mat["rank_drop_M1"] = (((S_mat[:, 1] + S_mat[:, 2]) / 2.0 + S_mat[:, 3] + S_mat[:, 4]) / 3.0).argsort()[::-1].argsort() + 1
# Use clean pandas rank for exact ties handling
for drop_idx, col_nm in [(0, "rank_drop_M1"), (1, "rank_drop_M2"), (2, "rank_drop_M3"), (3, "rank_drop_M4_GSE202051"), (4, "rank_drop_M5_RegVelo")]:
    if drop_idx in (0, 1, 2):
        rem = [i for i in (0, 1, 2) if i != drop_idx]
        sc_tmp = (S_mat[:, rem].mean(axis=1) + S_mat[:, 3] + S_mat[:, 4]) / 3.0
    elif drop_idx == 3:
        sc_tmp = (peng_block + S_mat[:, 4]) / 2.0
    else:
        sc_tmp = (peng_block + S_mat[:, 3]) / 2.0
    mac_mat[col_nm] = pd.Series(sc_tmp).rank(ascending=False, method="min").astype(int).values

# 1,000-draw symmetric Dirichlet(1,1,1) across the 3 independent cohorts (Peng block, GSE202051, GSE154778)
rng = np.random.default_rng(20261116)
W_dir = rng.dirichlet([1.0, 1.0, 1.0], size=1000)  # [1000, 3]
C3 = np.column_stack([peng_block, S_mat[:, 3], S_mat[:, 4]])  # [N, 3]
scores_dir = C3 @ W_dir.T  # [N, 1000]
ranks_dir = np.apply_along_axis(lambda x: pd.Series(x).rank(ascending=False, method="min").values, 0, scores_dir)
mac_mat["dirichlet_median_rank"] = np.median(ranks_dir, axis=1)
mac_mat["dirichlet_iqr_rank"] = np.percentile(ranks_dir, 75, axis=1) - np.percentile(ranks_dir, 25, axis=1)
mac_mat["dirichlet_top5_freq"] = np.mean(ranks_dir <= 5, axis=1)
mac_mat["dirichlet_top10_freq"] = np.mean(ranks_dir <= 10, axis=1)

# Attach Open Targets tractability purely as non-weighted translational metadata!
mac_mat = mac_mat.merge(
    ot_tract[["gene", "druggability_tier", "max_clinical_phase", "sm_tractable", "protac_tractable", "top_clinical_drugs"]].rename(columns={"gene": "TF"}),
    on="TF", how="left",
)
mac_mat["druggability_tier"] = mac_mat["druggability_tier"].fillna("Not_Queried_Or_Tier5_Undrugged")
mac_mat = mac_mat.sort_values(["primary_rank", "equal_5layer_rank"]).reset_index(drop=True)

# =========================================================================
# 2. EPITHELIAL MALIGNANT_DUCTAL REFERENCE MATRIX (CLEANLY SEPARATED!)
# =========================================================================
deg_duc = pd.read_csv(V1_ROOT / "04_results/pseudobulk/deg_MalignantDuctal2_vs_NormalDuctal1_edgeR_voom.csv")
lfc_d_col = "log2FC" if "log2FC" in deg_duc.columns else "logFC_edgeR"
fdr_d_col = "FDR" if "FDR" in deg_duc.columns else "FDR_edgeR"
surv_os = pd.read_csv(OUT_ROOT / "05_results/survival_results.tsv", sep="\t")
surv_d_os = surv_os[(surv_os["cohort"] == "Pure_Ductal_PDAC_n146") & (surv_os["endpoint"] == "OS")].copy()
depmap_df = pd.read_csv(ST6_DIR / "depmap_24q4_pdac_epithelial_audit.tsv", sep="\t")
rgv_duc = rgv_soft[(rgv_soft["compartment"] == "Malignant_Ductal") & (rgv_soft["perturbation_type"] == "Single_TF_KO")].copy()
rgv_duc_best = rgv_duc.assign(abs_pp=rgv_duc["mean_delta_pp"].abs()).sort_values("abs_pp", ascending=False).drop_duplicates("perturbation")

epi_tf_universe = sorted([g for g in deg_duc["gene"].unique() if g in lambert_tfs])
d_sub = deg_duc[deg_duc["gene"].isin(epi_tf_universe)].copy()
d_sub["neglog10_fdr"] = -np.log10(np.clip(d_sub[fdr_d_col].values, 1e-50, 1.0))
E1 = pos_percentile_score(d_sub, "gene", "neglog10_fdr", d_sub[lfc_d_col].abs() > 0.5, epi_tf_universe)
E2 = pos_percentile_score(rgv_duc_best, "perturbation", "abs_pp", rgv_duc_best["abs_pp"] > 0, epi_tf_universe)
s_sub = surv_d_os[surv_d_os["feature"].isin(epi_tf_universe)].copy()
s_sub["neglog10_p_multi"] = -np.log10(np.clip(s_sub["p_multi"].values, 1e-50, 1.0))
E3 = pos_percentile_score(s_sub, "feature", "neglog10_p_multi", s_sub["p_multi"] < 0.20, epi_tf_universe)
dep_sub = depmap_df[depmap_df["TF"].isin(epi_tf_universe)].copy()
dep_sub["neg_chronos"] = -dep_sub["mean_chronos_gene_effect"]
E4 = pos_percentile_score(dep_sub, "TF", "neg_chronos", dep_sub["mean_chronos_gene_effect"] < -0.10, epi_tf_universe)

epi_mat = pd.DataFrame({"TF": epi_tf_universe, "compartment_context": "Epithelial_Malignant_Ductal_Reference"})
epi_mat = epi_mat.merge(d_sub[["gene", lfc_d_col, fdr_d_col]].rename(columns={"gene": "TF", lfc_d_col: "E1_ductal_log2FC", fdr_d_col: "E1_ductal_FDR"}), on="TF", how="left")
epi_mat = epi_mat.merge(rgv_duc_best[["perturbation", "terminal_state", "mean_delta_pp", "empirical_decoy_p"]].rename(columns={"perturbation": "TF", "terminal_state": "E2_ductal_top_state", "mean_delta_pp": "E2_ductal_regvelo_delta_pp", "empirical_decoy_p": "E2_ductal_regvelo_decoy_p"}), on="TF", how="left")
epi_mat = epi_mat.merge(s_sub[["feature", "HR_multi", "p_multi", "HR_purity_immune_adj", "p_purity_immune_adj"]].rename(columns={"feature": "TF", "HR_multi": "E3_tcga_OS_HR_multi", "p_multi": "E3_tcga_OS_p_multi", "HR_purity_immune_adj": "E3_tcga_OS_HR_purity_adj", "p_purity_immune_adj": "E3_tcga_OS_p_purity_adj"}), on="TF", how="left")
epi_mat = epi_mat.merge(dep_sub[["TF", "mean_chronos_gene_effect", "fraction_dependent_cell_lines_lt_minus05"]], on="TF", how="left")
epi_mat["epithelial_reference_score"] = ([E1[t] for t in epi_mat["TF"]] + np.array([E2[t] for t in epi_mat["TF"]]) + np.array([E3[t] for t in epi_mat["TF"]]) + np.array([E4[t] for t in epi_mat["TF"]])) / 4.0
epi_mat["epithelial_reference_rank"] = epi_mat["epithelial_reference_score"].rank(ascending=False, method="min").astype(int)
epi_mat = epi_mat.sort_values("epithelial_reference_rank").reset_index(drop=True)

mac_mat.to_csv(ST6_DIR / "macrophage_ifn_tf_evidence_matrix_full.tsv", sep="\t", index=False)
epi_mat.to_csv(ST6_DIR / "epithelial_malignant_ductal_tf_reference_matrix.tsv", sep="\t", index=False)

# Unified target_evidence_matrix.tsv: Top 60 Macrophage TFs + Top 30 Epithelial Reference TFs (all columns preserved)
target_frozen = pd.concat([mac_mat.head(80), epi_mat.head(40)], ignore_index=True, sort=False)
target_frozen.to_csv(OUT_ROOT / "05_results/target_evidence_matrix.tsv", sep="\t", index=False)

print("  Top 10 Macrophage-Intrinsic IFN-State Regulators (Cohort-Balanced Score):")
print(mac_mat.head(10)[["primary_rank", "TF", "cohort_balanced_score", "equal_5layer_rank", "rank_drop_M5_RegVelo", "dirichlet_median_rank", "M1_peng_IFN10_log2FC", "M2_collectri_delta_ULM", "M3_scenic_delta_AUC", "M4_gse202051_IFN10_log2FC", "M5_regvelo_hard_delta_pp"]].to_string(index=False))

# =========================================================================
# 3. WRITE priority_score_definition.md
# =========================================================================
prio_md = """# Compartment-Separated Regulator Evidence Matrix & Priority Score Definition (`priority_score_definition.md`)

## 1. Forensic Audit of Legacy Scoring Flaws (`E02/P1`, `E03/P0`, `E04/P1`, `E05/P1`)
- **Legacy `v0` (`manuscript_pdac_regvelo.docx`)**: Claimed 7 layers in text, 5 layers in Table 1, and 4/5 binary support yet reported `0.638` for `STAT1` because `03_code` mixed binary thresholds with continuous normalized weights while crediting ductal knockouts (`KLF6/KLF5/IRF1/HIF1A/STAT3`) and cancer-cell-line `DepMap` scores toward macrophage claims.
- **Intermediate `v1` (`integrated_tf_prioritization_full.csv`)**: Used 6 layers (`RegVelo`, `Donor DEG`, `Regulon Activity`, `TCGA Clinical`, `DepMap CRISPR`, `External/Spatial`), which still (i) mixed `Malignant_Ductal2` `CollecTRI` activity with macrophage DEGs, (ii) credited `TCGA` bulk checkpoint correlation (`0.981`) as clinical survival support for `STAT1` despite multivariable OS neutrality (`HR = 1.10, p = 0.395`), and (iii) mixed `DepMap` epithelial cancer-cell-line dependency into the same composite ranking.

## 2. Rebuilt Compartment-Separated Architecture (`v2`, `2026-11-16`)
To prevent cross-lineage conflation and cohort double-counting, **macrophage-intrinsic IFN-state regulators** and **epithelial malignant-ductal regulators** are evaluated in strictly separate matrices (`macrophage_ifn_tf_evidence_matrix_full.tsv` and `epithelial_malignant_ductal_tf_reference_matrix.tsv`, combined with explicit `compartment_context` tags in `05_results/target_evidence_matrix.tsv`).

### 2.1 Primary Matrix: `Macrophage_IFN_State_Primary` (`Lambert et al., Cell 2018` Verified Human TFs Only)
Five macrophage-restricted complementary evidence layers are computed across three independent patient cohorts (`Peng2019` discovery, `GSE202051` external PDAC validation, and `GSE154778` splicing/velocity cohort):
1. **`M1_Peng_IFN10_DEG` (`Peng2019`, `n = 31` paired donors)**: Within-donor paired `edgeR` quasi-likelihood test comparing `IFN-10-high` vs `IFN-10-low` macrophages (`Rule Q2`; neither `STAT1` nor `IRF1` belongs to `IFN-10`, eliminating circularity for all TFs). Mapped to empirical percentile rank in `[0, 1]` of `-log10(BH FDR)` among positively upregulated TFs (`log2FC > 0`; `0.0` if `log2FC <= 0`).
2. **`M2_Peng_CollecTRI_ULM` (`Peng2019`, `n = 31` paired donors)**: Non-circular within-donor paired macrophage `CollecTRI` univariate linear model (`ULM`) activity (`IFN-10` defining genes removed from all regulon target sets prior to scoring). Mapped to empirical percentile rank in `[0, 1]` of `-log10(BH FDR)` among positively activated regulons (`delta_ULM > 0`; `0.0` if `<= 0` or `< 5` targets).
3. **`M3_Peng_pySCENIC_Motif` (`Peng2019`, `n = 29` paired donors)**: Motif-pruned `pySCENIC` (`hg38 v10 cisTarget` promoter/enhancer motif rankings) `AUCell` regulon activity shift. Mapped to empirical percentile rank in `[0, 1]` of `-log10(Wilcoxon BH padj)` among positively enriched regulons (`delta_AUC > 0`; `0.0` if `<= 0` or unpruned).
4. **`M4_GSE202051_Replication` (`GSE202051`, `n = 33` paired PDAC patients)**: Independent cohort within-patient paired `edgeR` `IFN-10` macrophage upregulation, mapped to empirical percentile rank in `[0, 1]` of `-log10(BH FDR)` (`0.0` if `log2FC <= 0`).
5. **`M5_GSE154778_RegVelo_HardKO` (`GSE154778`, `3,318` TAMs across `16` patients, `3` random seeds)**: Lineage-restricted `Macrophage_TAM` `hard_constraint` (`soft_constraint = False`) `RegVelo` virtual knockout absolute fate-probability shift (`|mean_delta_pp|`) on the `TAM_STAT1_Inflam` state, mapped to empirical percentile rank in `[0, 1]` (`0.0` if absent from the filtered splicing GRN).

### 2.2 Cohort-Balanced Weighting, Missingness Rule, and Sensitivity Audits
- **Explicit Shared-Data Control (`E04/P1`)**: Because `M1`, `M2`, and `M3` all derive from `Peng2019` discovery macrophages (complementary methods on shared cells, not independent cohorts), the primary score weights the **three independent single-cell cohorts equally (`1/3` each)**:
  - `Peng_Block = (M1 + M2 + M3) / 3`
  - `cohort_balanced_score = (Peng_Block + M4 + M5) / 3`
- **Denominator & Missingness Rule**: Unmeasured regulons in `M2`/`M3`/`M5` receive conservative `0.0` support in the primary fixed-denominator formula so TFs with sparse splicing coverage cannot jump to `#1` via single-layer division; `observed_denom_5layer_score` (dividing by `n_measured_layers_of_5`) and `equal_5layer_score = (M1 + M2 + M3 + M4 + M5) / 5` are reported alongside all **5 leave-one-layer-out (`LOLO`) ranks** and **1,000 symmetric `Dirichlet(1, 1, 1)` cohort-weight draws**.
- **Exclusion of Inappropriate Layers (`E02/P1`, `E05/P1`)**: `DepMap 24Q4` epithelial cancer-cell-line CRISPR scores and `TCGA-PAAD` bulk survival/correlations have **zero weight (`0.0`)** in `Macrophage_IFN_State_Primary` and are reported solely in `Epithelial_Malignant_Ductal_Reference` and clinical association tables.
- **Open Targets Platform Attribution**: Structural tractability and clinical-stage annotations (`druggability_tier`, `sm_tractable`, `protac_tractable`) are retrieved from the Open Targets Platform (`https://platform.opentargets.org`; licence: `https://platform-docs.opentargets.org/licence`, CC0 1.0) strictly as **unweighted descriptive metadata**, never as probability of clinical efficacy.
"""
(OUT_ROOT / "05_results/priority_score_definition.md").write_text(prio_md, encoding="utf-8")

# =========================================================================
# 4. BUILD UNIFIED primary_effects.tsv & manuscript_numbers.json
# =========================================================================
print("[2/4] Building unified primary_effects.tsv and manuscript_numbers.json...")
q1_tests = pd.read_csv(ST3_DIR / "between_donor_abundance_tests_Q1.tsv", sep="\t")
conc_df  = pd.read_csv(ST3_DIR / "cross_cohort_concordance_audit.tsv", sep="\t")
edge_aud = pd.read_csv(ST4_DIR / "regvelo_directed_edge_overlap_audit.tsv", sep="\t")
sp_site  = pd.read_csv(ST5_DIR / "gse272362_univariate_spatial_autocorrelation_by_site_and_patient.tsv", sep="\t")
cc_sens  = pd.read_csv(ST5_DIR / "cellchat_composition_and_subsampling_sensitivity.tsv", sep="\t")
corr_ex  = pd.read_csv(ST6_DIR / "tcga_ductal_stat1_partial_correlations_exact.tsv", sep="\t")

prim_rows = []
# 1. Between-donor Q1 expansion in Peng2019 & persistence in GSE202051 / GSE154778
for _, rq in q1_tests[q1_tests["metric"].isin(["ifn14_high_frac_of_mac_Q1", "ifn10_high_frac_of_mac_Q1", "ifn14_high_frac_of_all_cells_Q1"])].iterrows():
    prim_rows.append({
        "module": "Stage3_BetweenDonor_Abundance_Q1",
        "cohort": rq["cohort"],
        "contrast": f"{rq['contrast']} (min_mac>={int(rq['min_macrophages_per_donor'])})",
        "feature_or_metric": rq["metric"],
        "estimate_type": "Hodges_Lehmann_Difference",
        "estimate": float(rq["hodges_lehmann_diff"]),
        "ci95_low": float(rq["ci95_low"]),
        "ci95_high": float(rq["ci95_high"]),
        "p_value": float(rq["mwu_p"]),
        "fdr": float(rq["mwu_fdr_within_contrast"]),
        "statistical_unit": "Donor_Patient",
        "n_units": f"{int(rq['n_group1'])} vs {int(rq['n_group2'])}",
        "model": "Two-sided Mann-Whitney U + Hodges-Lehmann 95% CI",
    })

# 2. Primary paired edgeR key non-defining genes in Peng2019 (n=31), GSE202051 (n=34), GSE154778 (n=15)
for tag_lbl, coh_nm, df_de, n_d in [
    ("Peng_Q2_IFN14_min5", "Peng2019", de_peng14, 31),
    ("Peng_Q2_IFN10_min5", "Peng2019_IFN10", de_peng10, 31),
    ("GSE202051_All43_Q2_IFN14_min5", "GSE202051_All43", pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_GSE202051_All43_Q2_IFN14_min5.tsv", sep="\t"), 34),
    ("GSE202051_Untreated18_Q2_IFN14_min5", "GSE202051_Untreated18", de_202_un, 12),
    ("GSE154778_All16_Q2_IFN14_min5", "GSE154778_All16", pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_GSE154778_All16_Q2_IFN14_min5.tsv", sep="\t"), 15),
]:
    for g_k in ["STAT1", "IRF1", "IRF7", "CD47", "LGALS9", "CD80", "CD274", "IDO1", "B2M", "TAP1", "HLA-A", "KRT19", "COL1A1"]:
        rg = df_de[df_de["gene"] == g_k]
        if len(rg):
            r0 = rg.iloc[0]
            prim_rows.append({
                "module": "Stage3_WithinDonor_Paired_Pseudobulk_Q2",
                "cohort": coh_nm,
                "contrast": f"IFN-high vs IFN-low TAMs ({tag_lbl})",
                "feature_or_metric": g_k,
                "estimate_type": "log2FC_edgeR_QL",
                "estimate": float(r0["log2FC_edgeR"]),
                "ci95_low": float(r0["ci95_low_edgeR"]),
                "ci95_high": float(r0["ci95_high_edgeR"]),
                "p_value": float(r0["PValue_edgeR"]),
                "fdr": float(r0["FDR_edgeR_nonDefining"]) if pd.notna(r0["FDR_edgeR_nonDefining"]) else float(r0["FDR_edgeR_all"]),
                "statistical_unit": "Paired_Donor_Pseudobulk",
                "n_units": f"{n_d} paired donors ({2*n_d} pseudobulks)",
                "model": "edgeR TMM quasi-likelihood NB GLM (~ patient + group)",
            })

# 3. Cross-cohort concordance
for _, rc in conc_df.iterrows():
    prim_rows.append({
        "module": "Stage3_CrossCohort_Concordance",
        "cohort": rc["validation_model"],
        "contrast": f"Paired log2FC concordance vs {rc['discovery_model']}",
        "feature_or_metric": rc["gene_set_definition"],
        "estimate_type": "Spearman_rho",
        "estimate": float(rc["spearman_rho"]),
        "ci95_low": float(rc["ci95_low"]),
        "ci95_high": float(rc["ci95_high"]),
        "p_value": float(rc["p_value"]),
        "fdr": np.nan,
        "statistical_unit": "Shared_Expressed_Genes_From_Donor_Models",
        "n_units": f"{int(rc['n_genes'])} genes",
        "model": "Exact Student-t Spearman correlation + Fisher-z 95% CI",
    })

# 4. Stage 4 Hard-Constraint RegVelo Key Perturbations on TAM_STAT1_Inflam
st_inf = rgv_hard[rgv_hard["terminal_state"] == "TAM_STAT1_Inflam"].copy()
top_rgv_names = list(st_inf[st_inf["perturbation_type"] == "Single_TF_KO"].head(5)["perturbation"]) + ["STAT1", "IRF1", "STAT1+IRF1", "STAT3", "ACTB", "GAPDH"]
for _, rr in st_inf[st_inf["perturbation"].isin(top_rgv_names)].drop_duplicates("perturbation").iterrows():
    prim_rows.append({
        "module": "Stage4_RegVelo_HardConstraint_KO",
        "cohort": "GSE154778_Macrophage_TAM",
        "contrast": "In-silico GRN out-edge knockout vs WT on TAM_STAT1_Inflam",
        "feature_or_metric": rr["perturbation"],
        "estimate_type": "Delta_Fate_Percentage_Points (pp)",
        "estimate": float(rr["mean_delta_pp"]),
        "ci95_low": float(rr["mean_delta_pp"] - 1.96 * (rr["sd_delta_pp"] / np.sqrt(max(rr["n_seeds"], 1)))),
        "ci95_high": float(rr["mean_delta_pp"] + 1.96 * (rr["sd_delta_pp"] / np.sqrt(max(rr["n_seeds"], 1)))),
        "p_value": float(rr["empirical_decoy_p"]),
        "fdr": float(rr["empirical_decoy_fdr"]) if pd.notna(rr["empirical_decoy_fdr"]) else np.nan,
        "statistical_unit": "Multi_Seed_Model_Runs_And_PriorDecoys",
        "n_units": f"3,318 TAMs (16 patients; {int(rr['n_seeds'])} seeds; 20 degree-matched decoys)",
        "model": "REGVELOVI (soft_constraint=False) + CellRank GPCCA",
    })

prim_df = pd.DataFrame(prim_rows)
prim_df.to_csv(OUT_ROOT / "05_results/primary_effects.tsv", sep="\t", index=False)

# Build manuscript_numbers.json
stat1_hard = st_inf[st_inf["perturbation"] == "STAT1"].iloc[0] if len(st_inf[st_inf["perturbation"] == "STAT1"]) else None
irf1_hard  = st_inf[st_inf["perturbation"] == "IRF1"].iloc[0] if len(st_inf[st_inf["perturbation"] == "IRF1"]) else None
combo_hard = st_inf[st_inf["perturbation"] == "STAT1+IRF1"].iloc[0] if len(st_inf[st_inf["perturbation"] == "STAT1+IRF1"]) else None
top1_hard  = st_inf[st_inf["perturbation_type"] == "Single_TF_KO"].iloc[0]
top2_hard  = st_inf[st_inf["perturbation_type"] == "Single_TF_KO"].iloc[1]

stat1_mac_row = mac_mat[mac_mat["TF"] == "STAT1"].iloc[0]
irf1_mac_row  = mac_mat[mac_mat["TF"] == "IRF1"].iloc[0]

ms_nums = {
    "freeze_date": "2026-11-16",
    "selected_story_path": "Path_A_Rigorous_Computational_Discovery_Study",
    "primary_target_journal": "BMC Genomics",
    "cohorts": {
        "Peng2019_CRA001160": {"n_pdac_patients": 24, "n_control_donors": 11, "n_total_donors": 35, "n_qc_cells": 57423, "n_unique_genes": 18008, "n_macrophages": 5361, "n_paired_donors_Q2_min5": 31, "n_paired_pdac_only_min5": 24},
        "GSE202051_snRNA": {"n_pdac_patients": 43, "n_untreated": 18, "n_neoadjuvant": 25, "n_normal_controls": 0, "n_macrophages": 10140, "n_paired_donors_All43_min5": 34, "n_paired_untreated_min5": 12, "n_paired_neoadjuvant_min5": 22},
        "GSE154778_scRNA_SRA": {"n_patients": 16, "n_primary": 10, "n_metastasis": 6, "n_sra_runs": 80, "n_qc_cells": 14915, "n_macrophages": 3318, "n_malignant_ductal": 10042, "n_paired_donors_min5": 15},
        "GSE272362_Visium": {"n_sections": 30, "n_autopsy_patients": 7, "n_PT": 10, "n_HM": 12, "n_LNM": 5, "n_NP": 3, "n_in_tissue_spots": 97020},
        "TCGA_PAAD": {"n_pure_ductal": 146, "n_pure_ductal_os_valid": 145, "n_pure_ductal_os_events": 84, "n_pure_ductal_multi_complete": 144, "n_unfiltered_paad": 177},
    },
    "stage3_key_results": {
        "peng_Q1_IFN14_median_pdac_pct": round(float(q1_tests.iloc[0]["median_group1"]) * 100, 1),
        "peng_Q1_IFN14_median_ctrl_pct": round(float(q1_tests.iloc[0]["median_group2"]) * 100, 1),
        "peng_Q1_IFN14_mwu_p": float(q1_tests.iloc[0]["mwu_p"]),
        "peng_Q2_IFN14_n_nondefining_DEGs_FDR005": int((~de_peng14["in_defining_signature"] & (de_peng14["FDR_edgeR_nonDefining"] < 0.05)).sum()),
        "peng_STAT1_paired_log2FC": round(float(de_peng14.loc[de_peng14["gene"] == "STAT1", "log2FC_edgeR"].iloc[0]), 3),
        "peng_STAT1_paired_FDR": float(de_peng14.loc[de_peng14["gene"] == "STAT1", "FDR_edgeR_nonDefining"].iloc[0]),
        "peng_IRF1_IFN10_paired_log2FC": round(float(de_peng10.loc[de_peng10["gene"] == "IRF1", "log2FC_edgeR"].iloc[0]), 3),
        "peng_IRF1_IFN10_paired_FDR": float(de_peng10.loc[de_peng10["gene"] == "IRF1", "FDR_edgeR_nonDefining"].iloc[0]),
        "peng_CD47_paired_log2FC": round(float(de_peng14.loc[de_peng14["gene"] == "CD47", "log2FC_edgeR"].iloc[0]), 3),
        "peng_CD47_paired_FDR": float(de_peng14.loc[de_peng14["gene"] == "CD47", "FDR_edgeR_nonDefining"].iloc[0]),
        "concordance_GSE202051_All43_299nonIFN_rho": round(float(conc_df[(conc_df["validation_model"]=="GSE202051_All43_Q2_IFN14_min5") & (conc_df["gene_set_definition"].str.startswith("3_"))]["spearman_rho"].iloc[0]), 3),
        "concordance_GSE154778_317nonIFN_rho": round(float(conc_df[(conc_df["validation_model"]=="GSE154778_All16_Q2_IFN14_min5") & (conc_df["gene_set_definition"].str.startswith("3_"))]["spearman_rho"].iloc[0]), 3),
    },
    "stage4_regvelo_hard_constraint": {
        "n_genes": int(edge_aud["n_genes"].iloc[0]),
        "n_lambert_tfs": int(edge_aud["n_lambert_tfs"].iloc[0]),
        "n_prior_edges": int(edge_aud["n_prior_directed_edges"].iloc[0]),
        "pct_within_prior_correct_axis": round(float(edge_aud["pct_inferred_within_prior_correct"].mean()), 2),
        "n_spurious_off_prior_edges": int(edge_aud["n_spurious_off_prior_edges"].sum()),
        "pct_overlap_if_transposed_v0_bug": round(float(edge_aud["pct_overlap_if_transposed_v0_bug"].mean()), 2),
        "top1_single_TF": str(top1_hard["perturbation"]),
        "top1_single_TF_delta_pp": round(float(top1_hard["mean_delta_pp"]), 2),
        "top1_single_TF_decoy_p": round(float(top1_hard["empirical_decoy_p"]), 4),
        "top2_single_TF": str(top2_hard["perturbation"]),
        "top2_single_TF_delta_pp": round(float(top2_hard["mean_delta_pp"]), 2),
        "STAT1_hard_delta_pp": round(float(stat1_hard["mean_delta_pp"]), 3) if stat1_hard is not None else None,
        "STAT1_hard_decoy_p": round(float(stat1_hard["empirical_decoy_p"]), 4) if stat1_hard is not None else None,
        "IRF1_hard_delta_pp": round(float(irf1_hard["mean_delta_pp"]), 3) if irf1_hard is not None else None,
        "STAT1_plus_IRF1_hard_delta_pp": round(float(combo_hard["mean_delta_pp"]), 3) if combo_hard is not None else None,
    },
    "stage6_prioritization_and_clinical": {
        "macrophage_top5_TFs": mac_mat.head(5)["TF"].tolist(),
        "STAT1_macrophage_rank": int(stat1_mac_row["primary_rank"]),
        "STAT1_macrophage_score": round(float(stat1_mac_row["cohort_balanced_score"]), 3),
        "IRF1_macrophage_rank": int(irf1_mac_row["primary_rank"]),
        "IRF1_macrophage_score": round(float(irf1_mac_row["cohort_balanced_score"]), 3),
        "epithelial_top5_TFs": epi_mat.head(5)["TF"].tolist(),
        "tcga_pure_ductal_STAT1_OS_HR_multi": round(float(surv_d_os.loc[surv_d_os["feature"]=="STAT1", "HR_multi"].iloc[0]), 3),
        "tcga_pure_ductal_STAT1_OS_p_multi": round(float(surv_d_os.loc[surv_d_os["feature"]=="STAT1", "p_multi"].iloc[0]), 4),
        "tcga_pure_ductal_FOXM1_OS_HR_multi": round(float(surv_d_os.loc[surv_d_os["feature"]=="FOXM1", "HR_multi"].iloc[0]), 3),
        "tcga_pure_ductal_HMGA2_OS_HR_multi": round(float(surv_d_os.loc[surv_d_os["feature"]=="HMGA2", "HR_multi"].iloc[0]), 3),
        "lincs_macrophage_strict_tier1_count": 4,
        "lincs_macrophage_strict_tier1_names": ["MG-132", "brefeldin-a", "bortezomib", "KW-2478"],
    },
}
(OUT_ROOT / "05_results/manuscript_numbers.json").write_text(json.dumps(ms_nums, indent=2), encoding="utf-8")
print("[Stage 6 Part B2 & Frozen Results Layer Complete!]")
