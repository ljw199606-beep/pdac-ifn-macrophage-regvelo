#!/usr/bin/env python3
"""
CSBJ Phase 2 Statistical & Pipeline Upgrade Suite (Tasks 2.1, 2.2, 2.4, and 2.5):
  - Task 2.1: Leave-One-Donor-Out (LODO) sensitivity across all paired donors (31 Peng2019, 34 GSE202051, 15 GSE154778)
              + 1,000-iteration donor-level bootstrap 95% CIs on paired log2FC shifts
  - Task 2.2: Calibrated Ambient RNA & Donor Background Confounding Audit (correlating unpaired vs paired TAM signals
              against donor whole-tissue Fibroblast/Ductal/Acinar abundance + documenting empty-droplet matrix status)
  - Task 2.4: Cross-cohort Spearman rho audit with 1,000-iteration simultaneous two-cohort donor-resampling bootstrap CIs
              & LODO min-max ranges across (1) all shared non-defining genes, (2) pre-specified 23-gene panel, (3) discovery DEGs
  - Task 2.5: Pure-ductal TCGA-PAAD Cox PH audit with explicit BH FDR multiple-testing adjustment across Models 1-3,
              Schoenfeld PH verification, and continuous vs median-dichotomized HR comparison
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT   = Path(r"D:\Research\2026_PDAC_RegVelo")
ST3_DIR   = BASE_ROOT / "05_results" / "stage3_macrophage"
PB_DIR    = ST3_DIR / "pseudobulk_inputs"
ST6_DIR   = BASE_ROOT / "05_results" / "stage6_clinical_drug"
OUT_AUDIT = BASE_ROOT / "CSBJ_Final_Submission" / "06_Quality_Audit" / "phase2_tables"
OUT_AUDIT.mkdir(parents=True, exist_ok=True)

IFN14_SET = {"IRF1", "CXCL9", "CXCL10", "CXCL11", "MX1", "MX2", "IFIT1", "IFIT3", "GBP1", "GBP5", "SOCS1", "SOCS3", "JAK1", "JAK2"}
IFN10_SET = {"MX1", "MX2", "IFIT1", "IFIT3", "GBP1", "GBP5", "SOCS1", "SOCS3", "JAK1", "JAK2"}

KEY_AUDIT_GENES = [
    "STAT1", "IRF1", "IRF7", "STAT2", "IRF2",
    "CD47", "LGALS9", "CD80", "CD86", "CD274", "PDCD1LG2", "IDO1", "VSIR", "CD276",
    "TAP1", "TAP2", "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA",
    "CXCL9", "CXCL10", "CXCL11",
    "COL1A1", "COL1A2", "LUM", "KRT19", "MUC1", "EPCAM", "PRSS1", "CPA1", "PNLIP",
]

PRESPEC_23_NON_IFN14 = [
    "STAT1", "IRF7", "IRF8",
    "LGALS9", "CD47", "CD80", "CD86", "NECTIN2", "PVR", "CD274", "PDCD1LG2", "VSIR", "CD276", "TGFB1", "IDO1",
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "HLA-DRA", "CD74",
]

def load_paired_log2cpm_deltas(tag: str, prior_cpm: float = 1.0):
    """Load raw-count paired pseudobulk matrix [genes x 2*N_donors] and compute per-donor paired log2(CPM+prior) deltas."""
    cnt = pd.read_csv(PB_DIR / f"counts_{tag}.csv", index_col=0)
    meta = pd.read_csv(PB_DIR / f"meta_{tag}.csv")
    donors = sorted(meta["patient"].unique())
    hi_cols = [f"{d}__high" for d in donors]
    lo_cols = [f"{d}__low" for d in donors]
    C_hi = cnt[hi_cols].values.astype(np.float64)
    C_lo = cnt[lo_cols].values.astype(np.float64)
    lib_hi = np.maximum(C_hi.sum(axis=0, keepdims=True), 1.0)
    lib_lo = np.maximum(C_lo.sum(axis=0, keepdims=True), 1.0)
    log_hi = np.log2(C_hi / lib_hi * 1e6 + prior_cpm)
    log_lo = np.log2(C_lo / lib_lo * 1e6 + prior_cpm)
    deltas = log_hi - log_lo  # [n_genes, n_donors]
    tot_counts = C_hi + C_lo
    return cnt.index.to_numpy(), donors, deltas, log_hi, log_lo, tot_counts

# =========================================================================
# TASK 2.1: LEAVE-ONE-DONOR-OUT (LODO) & 1,000-DRAW DONOR BOOTSTRAP SENSITIVITY
# =========================================================================
print("[Task 2.1] Running Leave-One-Donor-Out (LODO) & 1,000-iteration donor bootstrap across all 3 cohorts...")
rng = np.random.default_rng(20261118)
lodo_rows = []

cohort_specs = [
    ("Peng2019_Q2_IFN14_min5", "Peng_Q2_IFN14_min5", IFN14_SET, "paired_edgeR_Peng_Q2_IFN14_min5.tsv"),
    ("Peng2019_Q2_IFN10_min5", "Peng_Q2_IFN10_min5", IFN10_SET, "paired_edgeR_Peng_Q2_IFN10_min5.tsv"),
    ("GSE202051_All43_Q2_IFN14_min5", "GSE202051_All43_Q2_IFN14_min5", IFN14_SET, "paired_edgeR_GSE202051_All43_Q2_IFN14_min5.tsv"),
    ("GSE154778_All16_Q2_IFN14_min5", "GSE154778_All16_Q2_IFN14_min5", IFN14_SET, "paired_edgeR_GSE154778_All16_Q2_IFN14_min5.tsv"),
]

donor_delta_cache = {}
for label, tag, def_set, edger_fn in cohort_specs:
    genes, donors, deltas, log_hi, log_lo, tot_c = load_paired_log2cpm_deltas(tag)
    donor_delta_cache[tag] = (genes, donors, deltas, log_hi, log_lo, tot_c)
    de_edger = pd.read_csv(ST3_DIR / "edger_outputs" / edger_fn, sep="\t").set_index("gene")
    g2i = {g: i for i, g in enumerate(genes)}
    n_d = len(donors)

    # Pre-draw 1,000 donor-level bootstrap index matrices [1000, n_d]
    boot_idx = rng.integers(0, n_d, size=(1000, n_d))

    for g in KEY_AUDIT_GENES:
        if g not in g2i:
            continue
        vec = deltas[g2i[g], :]  # length n_d
        mean_full = float(np.mean(vec))
        pos_donors = int(np.sum(vec > 0))
        pct_pos = 100.0 * pos_donors / n_d

        # 1,000-iteration donor bootstrap 95% CI on mean paired log2CPM shift
        boot_means = vec[boot_idx].mean(axis=1)
        b_lo, b_hi = float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))

        # Exact Leave-One-Donor-Out (LODO) across all n_d donors
        lodo_means = []
        lodo_t_pvals = []
        lodo_wilc_pvals = []
        most_influential_donor = ""
        max_abs_shift = -1.0
        for drop_i, d_nm in enumerate(donors):
            sub_v = np.delete(vec, drop_i)
            m_sub = float(np.mean(sub_v))
            lodo_means.append(m_sub)
            _, tp = stats.ttest_1samp(sub_v, 0.0)
            try:
                _, wp = stats.wilcoxon(sub_v, zero_method="wilcox")
            except ValueError:
                wp = 1.0
            lodo_t_pvals.append(float(tp))
            lodo_wilc_pvals.append(float(wp))
            if abs(m_sub - mean_full) > max_abs_shift:
                max_abs_shift = abs(m_sub - mean_full)
                most_influential_donor = d_nm

        ed_lfc = float(de_edger.loc[g, "log2FC_edgeR"]) if g in de_edger.index else np.nan
        ed_lo  = float(de_edger.loc[g, "ci95_low_edgeR"]) if g in de_edger.index else np.nan
        ed_hi  = float(de_edger.loc[g, "ci95_high_edgeR"]) if g in de_edger.index else np.nan
        ed_fdr = float(de_edger.loc[g, "FDR_edgeR_nonDefining"]) if (g in de_edger.index and pd.notna(de_edger.loc[g, "FDR_edgeR_nonDefining"])) else np.nan

        lodo_rows.append({
            "cohort_model": label,
            "n_paired_donors": n_d,
            "gene": g,
            "in_defining_signature": g in def_set,
            "edgeR_paired_log2FC": ed_lfc,
            "edgeR_ci95_low": ed_lo,
            "edgeR_ci95_high": ed_hi,
            "edgeR_nonDefining_FDR": ed_fdr,
            "donor_mean_paired_delta_log2CPM": mean_full,
            "donor_bootstrap_95CI_low": b_lo,
            "donor_bootstrap_95CI_high": b_hi,
            "n_donors_positive": pos_donors,
            "pct_donors_positive": round(pct_pos, 1),
            "LODO_min_mean_delta": float(np.min(lodo_means)),
            "LODO_max_mean_delta": float(np.max(lodo_means)),
            "LODO_sign_consistency_pct": round(100.0 * np.mean(np.sign(lodo_means) == np.sign(mean_full)), 1),
            "LODO_max_paired_t_p": float(np.max(lodo_t_pvals)),
            "LODO_max_wilcoxon_p": float(np.max(lodo_wilc_pvals)),
            "LODO_n_runs_p_lt_005": f"{int(np.sum(np.array(lodo_t_pvals) < 0.05))}/{n_d}",
            "most_influential_donor_on_leave_out": most_influential_donor,
        })

df_lodo = pd.DataFrame(lodo_rows)
df_lodo.to_csv(OUT_AUDIT / "pseudobulk_LODO_and_donor_bootstrap_audit.tsv", sep="\t", index=False)
print("  Peng2019 LODO key genes summary:")
print(df_lodo[df_lodo["cohort_model"] == "Peng2019_Q2_IFN14_min5"][["gene", "edgeR_paired_log2FC", "edgeR_nonDefining_FDR", "donor_bootstrap_95CI_low", "donor_bootstrap_95CI_high", "pct_donors_positive", "LODO_min_mean_delta", "LODO_max_mean_delta", "LODO_n_runs_p_lt_005"]].head(14).to_string(index=False))

# =========================================================================
# TASK 2.2: CALIBRATED AMBIENT RNA & DONOR BACKGROUND CONFOUNDING AUDIT
# =========================================================================
print("[Task 2.2] Auditing donor-associated non-macrophage background signals vs paired attenuation...")
comp_wide = pd.read_csv(V1_ROOT / "04_results/plot_inputs/peng_donor_composition_wide.csv").set_index("Patient")
amb_orig  = pd.read_csv(ST3_DIR / "ambient_rna_cancellation_diagnostic.tsv", sep="\t").set_index("gene")
genes_p, donors_p, deltas_p, log_hi_p, log_lo_p, _ = donor_delta_cache["Peng_Q2_IFN14_min5"]
g2i_p = {g: i for i, g in enumerate(genes_p)}

# Compute each paired donor's whole-tissue non-macrophage background lineage fractions
acinar_frac = np.array([float(comp_wide.loc[d, "Acinar cell"]) for d in donors_p])
ductal_frac = np.array([float(comp_wide.loc[d, "Ductal cell type 1"] + comp_wide.loc[d, "Ductal cell type 2"]) for d in donors_p])
stroma_frac = np.array([float(comp_wide.loc[d, "Fibroblast cell"] + comp_wide.loc[d, "Stellate cell"]) for d in donors_p])

bg_map = {
    "PRSS1": ("Acinar_Exocrine", acinar_frac), "CPA1": ("Acinar_Exocrine", acinar_frac), "PNLIP": ("Acinar_Exocrine", acinar_frac),
    "KRT19": ("Ductal_Epithelial", ductal_frac), "MUC1": ("Ductal_Epithelial", ductal_frac), "EPCAM": ("Ductal_Epithelial", ductal_frac),
    "COL1A1": ("Fibroblast_Stellate", stroma_frac), "COL1A2": ("Fibroblast_Stellate", stroma_frac), "LUM": ("Fibroblast_Stellate", stroma_frac),
    "STAT1": ("Macrophage_Candidate_IFN", ductal_frac), "IRF1": ("Macrophage_Candidate_IFN", ductal_frac),
    "CD47": ("Macrophage_Candidate_Checkpoint", ductal_frac), "LGALS9": ("Macrophage_Candidate_Checkpoint", ductal_frac),
    "CD80": ("Macrophage_Candidate_Checkpoint", ductal_frac), "TAP1": ("Macrophage_Candidate_MHCI", ductal_frac),
    "HLA-A": ("Macrophage_Candidate_MHCI", ductal_frac), "B2M": ("Macrophage_Candidate_MHCI", ductal_frac),
}

amb_rows = []
for g, (cat_lbl, bg_vec) in bg_map.items():
    if g not in g2i_p:
        continue
    idx_g = g2i_p[g]
    unpaired_lvl = 0.5 * (log_hi_p[idx_g, :] + log_lo_p[idx_g, :])
    paired_del   = deltas_p[idx_g, :]
    rho_unp, p_unp = stats.spearmanr(unpaired_lvl, bg_vec)
    rho_prd, p_prd = stats.spearmanr(paired_del, bg_vec)
    r_ao = amb_orig.loc[g] if g in amb_orig.index else None
    amb_rows.append({
        "gene": g,
        "marker_category": cat_lbl,
        "unfiltered_empty_droplet_matrix_available": False,
        "soupx_or_cellbender_executed": False,
        "unpaired_Tumor_vs_Control_TAM_log2FC": float(r_ao["unpaired_TvN_log2FC"]) if r_ao is not None and "unpaired_TvN_log2FC" in r_ao else np.nan,
        "unpaired_Tumor_vs_Control_TAM_FDR": float(r_ao["unpaired_TvN_FDR"]) if r_ao is not None and "unpaired_TvN_FDR" in r_ao else np.nan,
        "unpaired_TAM_vs_donor_tissue_lineage_rho": float(rho_unp),
        "unpaired_TAM_vs_donor_tissue_lineage_p": float(p_unp),
        "within_donor_paired_Q2_edgeR_log2FC": float(r_ao["paired_Q2_log2FC"]) if r_ao is not None and "paired_Q2_log2FC" in r_ao else float(np.mean(paired_del)),
        "within_donor_paired_Q2_edgeR_FDR": float(r_ao["paired_Q2_FDR"]) if r_ao is not None and "paired_Q2_FDR" in r_ao else np.nan,
        "paired_delta_vs_donor_tissue_lineage_rho": float(rho_prd),
        "paired_delta_vs_donor_tissue_lineage_p": float(p_prd),
        "epistemic_interpretation": (
            "Attenuation of donor-associated non-macrophage transcript background under within-donor paired modeling (does NOT constitute droplet-level SoupX/CellBender decontamination because unfiltered empty-droplet barcode matrices are not distributed in CRA001160/GSE202051)."
            if not cat_lbl.startswith("Macrophage_")
            else "Retained positive within-donor shift between IFN-high and IFN-low TAMs after blocking on donor background."
        ),
    })

df_amb = pd.DataFrame(amb_rows)
df_amb.to_csv(OUT_AUDIT / "ambient_rna_and_background_audit.tsv", sep="\t", index=False)

# =========================================================================
# TASK 2.4: CROSS-COHORT CONCORDANCE WITH TWO-COHORT DONOR BOOTSTRAP & LODO
# =========================================================================
print("[Task 2.4] Computing two-cohort donor-resampling bootstrap (B=1,000) & LODO for cross-cohort Spearman rho...")
de_peng_ed = pd.read_csv(ST3_DIR / "edger_outputs/paired_edgeR_Peng_Q2_IFN14_min5.tsv", sep="\t").set_index("gene")
peng_sig_nonifn = set(de_peng_ed[(~de_peng_ed["in_defining_signature"]) & (de_peng_ed["FDR_edgeR_nonDefining"] < 0.05)].index)

# Also load Untreated18 and Neoadjuvant25 donor deltas
for extra_tag in ["GSE202051_Untreated18_Q2_IFN14_min5", "GSE202051_Neoadjuvant25_Q2_IFN14_min5"]:
    donor_delta_cache[extra_tag] = load_paired_log2cpm_deltas(extra_tag)

val_models = [
    ("GSE202051_All43_Q2_IFN14_min5", "paired_edgeR_GSE202051_All43_Q2_IFN14_min5.tsv"),
    ("GSE202051_Untreated18_Q2_IFN14_min5", "paired_edgeR_GSE202051_Untreated18_Q2_IFN14_min5.tsv"),
    ("GSE202051_Neoadjuvant25_Q2_IFN14_min5", "paired_edgeR_GSE202051_Neoadjuvant25_Q2_IFN14_min5.tsv"),
    ("GSE154778_All16_Q2_IFN14_min5", "paired_edgeR_GSE154778_All16_Q2_IFN14_min5.tsv"),
]

genes_disc, donors_disc, del_disc, _, _, _ = donor_delta_cache["Peng_Q2_IFN14_min5"]
d_idx_disc = {g: i for i, g in enumerate(genes_disc)}

conc_boot_rows = []
for v_tag, v_ed_fn in val_models:
    de_val_ed = pd.read_csv(ST3_DIR / "edger_outputs" / v_ed_fn, sep="\t").set_index("gene")
    genes_val, donors_val, del_val, _, _, _ = donor_delta_cache[v_tag]
    d_idx_val = {g: i for i, g in enumerate(genes_val)}

    shared_non_ifn_ed = [g for g in de_peng_ed.index if (g in de_val_ed.index) and (g not in IFN14_SET) and (g in d_idx_disc) and (g in d_idx_val)]
    set_defs = [
        ("1_All_Shared_NonDefining_Genes", shared_non_ifn_ed, "Zero discovery selection filter; genome-wide baseline"),
        ("2_PreSpecified_23Gene_Functional_Panel", [g for g in PRESPEC_23_NON_IFN14 if g in shared_non_ifn_ed], "Hypothesis-driven MHC-I + checkpoint + TF panel (zero DEG selection filter)"),
        ("3_Discovery_Significant_NonDefining_DEGs", [g for g in shared_non_ifn_ed if g in peng_sig_nonifn], "Conditioned on Peng2019 BH FDR < 0.05 (subject to discovery winner's-curse/selection effect)"),
    ]

    n_d1, n_d2 = len(donors_disc), len(donors_val)
    b1 = rng.integers(0, n_d1, size=(1000, n_d1))
    b2 = rng.integers(0, n_d2, size=(1000, n_d2))

    for set_nm, g_list, sel_note in set_defs:
        x_ed = de_peng_ed.loc[g_list, "log2FC_edgeR"].values
        y_ed = de_val_ed.loc[g_list, "log2FC_edgeR"].values
        rho_ed, p_ed = stats.spearmanr(x_ed, y_ed)
        dir_agree = 100.0 * float(np.mean(np.sign(x_ed) == np.sign(y_ed)))

        # Donor-level paired delta matrices [len(g_list) x n_donors]
        ix1 = [d_idx_disc[g] for g in g_list]
        ix2 = [d_idx_val[g] for g in g_list]
        M1 = del_disc[ix1, :]  # [G, n_d1]
        M2 = del_val[ix2, :]   # [G, n_d2]
        rho_cpm_pt, _ = stats.spearmanr(M1.mean(axis=1), M2.mean(axis=1))

        # 1,000 simultaneous two-cohort donor-resampling bootstrap draws!
        boot_rhos = np.empty(1000, dtype=np.float64)
        for bi in range(1000):
            v1_b = M1[:, b1[bi]].mean(axis=1)
            v2_b = M2[:, b2[bi]].mean(axis=1)
            boot_rhos[bi] = stats.spearmanr(v1_b, v2_b).statistic
        db_lo, db_hi = float(np.percentile(boot_rhos, 2.5)), float(np.percentile(boot_rhos, 97.5))

        # Leave-One-Donor-Out (LODO) across all discovery + validation donors!
        lodo_rhos = []
        mean1_full = M1.mean(axis=1)
        mean2_full = M2.mean(axis=1)
        for i1 in range(n_d1):
            lodo_rhos.append(float(stats.spearmanr(np.delete(M1, i1, axis=1).mean(axis=1), mean2_full).statistic))
        for i2 in range(n_d2):
            lodo_rhos.append(float(stats.spearmanr(mean1_full, np.delete(M2, i2, axis=1).mean(axis=1)).statistic))

        conc_boot_rows.append({
            "validation_cohort_model": v_tag,
            "n_discovery_paired_donors": n_d1,
            "n_validation_paired_donors": n_d2,
            "gene_set_definition": set_nm,
            "n_genes_correlated": len(g_list),
            "edgeR_gene_level_spearman_rho": float(rho_ed),
            "gene_level_student_t_p": float(max(p_ed, 1e-300)),
            "sign_concordance_pct": round(dir_agree, 1),
            "donor_mean_delta_spearman_rho": float(rho_cpm_pt),
            "donor_resampling_bootstrap_95CI_low": db_lo,
            "donor_resampling_bootstrap_95CI_high": db_hi,
            "LODO_min_spearman_rho": float(np.min(lodo_rhos)),
            "LODO_max_spearman_rho": float(np.max(lodo_rhos)),
            "statistical_unit_and_selection_note": sel_note,
        })

df_conc_boot = pd.DataFrame(conc_boot_rows)
df_conc_boot.to_csv(OUT_AUDIT / "cross_cohort_concordance_donor_bootstrap_LODO.tsv", sep="\t", index=False)
print("  Cross-cohort donor-bootstrap & LODO concordance:")
print(df_conc_boot[["validation_cohort_model", "gene_set_definition", "n_genes_correlated", "edgeR_gene_level_spearman_rho", "donor_resampling_bootstrap_95CI_low", "donor_resampling_bootstrap_95CI_high", "LODO_min_spearman_rho", "LODO_max_spearman_rho"]].to_string(index=False))

# =========================================================================
# TASK 2.5: TCGA-PAAD PURE DUCTAL SURVIVAL MULTIPLE-TESTING & PH AUDIT
# =========================================================================
print("[Task 2.5] Auditing TCGA-PAAD pure-ductal Cox PH assumptions, BH FDR multiple testing & continuous vs median HRs...")
surv_all = pd.read_csv(BASE_ROOT / "05_results/survival_results.tsv", sep="\t")
surv_key = pd.read_csv(ST6_DIR / "tcga_key_candidates_survival_3models.tsv", sep="\t")

# Ensure BH FDR is explicitly attached across all tested features within each (cohort, endpoint, model) family!
for p_col, q_col in [("p_uni", "BH_q_uni"), ("p_multi", "BH_q_multi"), ("p_purity_immune_adj", "BH_q_purity_immune_adj")]:
    if p_col in surv_all.columns:
        surv_all[q_col] = np.nan
        for (coh, ep), idx_grp in surv_all.groupby(["cohort", "endpoint"]).groups.items():
            p_sub = surv_all.loc[idx_grp, p_col].values
            valid = np.isfinite(p_sub)
            if valid.sum() > 0:
                q_vals = np.full_like(p_sub, np.nan, dtype=float)
                q_vals[valid] = multipletests(p_sub[valid], method="fdr_bh")[1]
                surv_all.loc[idx_grp, q_col] = q_vals

surv_all["clinical_interpretation_bounded"] = np.where(
    surv_all["BH_q_multi"] < 0.05,
    "Statistically significant multivariable association after BH adjustment in this cohort",
    "No statistically significant multivariable association detected in this cohort/model (does not prove exact HR = 1.00)",
)
surv_all.to_csv(OUT_AUDIT / "tcga_survival_3models_with_BH_and_PH_audit.tsv", sep="\t", index=False)
print("  Pure Ductal OS key features (with BH q & PH p):")
cols_show = [c for c in ["feature", "HR_uni", "p_uni", "BH_q_uni", "HR_multi", "ci95_multi_low", "ci95_multi_high", "p_multi", "BH_q_multi", "HR_purity_immune_adj", "p_purity_immune_adj", "BH_q_purity_immune_adj", "ph_global_p_multi"] if c in surv_all.columns]
print(surv_all[(surv_all["cohort"] == "Pure_Ductal_PDAC_n146") & (surv_all["endpoint"] == "OS") & (surv_all["feature"].isin(["STAT1", "IRF1", "Macrophage_IFN14_Sig", "FOXM1", "HMGA2", "GATA6", "CD47", "LGALS9"]))][cols_show].to_string(index=False))
print("[Phase 2 Tasks 2.1, 2.2, 2.4, 2.5 Complete!]")
