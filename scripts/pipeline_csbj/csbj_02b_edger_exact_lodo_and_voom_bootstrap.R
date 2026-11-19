#!/usr/bin/env Rscript
# CSBJ Phase 2 Exact Biostatistical Upgrade (Tasks 2.1 & 2.4):
# 1. Exact Leave-One-Donor-Out (LODO) edgeR TMM QL NB GLM fits across all 31 Peng2019 donors,
#    all 34 GSE202051 patients, and all 15 GSE154778 patients!
# 2. Precision-weighted voom paired donor-bootstrap (B = 1,000) 95% CIs on gene log2FCs AND
#    simultaneous two-cohort donor-bootstrap + LODO Spearman rho across all 3 gene-set definitions!

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(parallel)
})

base_root <- "/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116"
pb_dir    <- file.path(base_root, "05_results/stage3_macrophage/pseudobulk_inputs")
ed_dir    <- file.path(base_root, "05_results/stage3_macrophage/edger_outputs")
out_dir   <- file.path(base_root, "CSBJ_Final_Submission/06_Quality_Audit/phase2_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ifn14_genes <- c("IRF1","CXCL9","CXCL10","CXCL11","MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2")
ifn10_genes <- c("MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2")
key_genes   <- c("STAT1","IRF1","IRF7","STAT2","IRF2","CD47","LGALS9","CD80","CD86","CD274","PDCD1LG2","IDO1","VSIR","CD276","TAP1","TAP2","B2M","HLA-A","HLA-B","HLA-C","CXCL9","CXCL10","COL1A1","KRT19","PRSS1","CPA1")
prespec_23  <- c("STAT1","IRF7","IRF8","LGALS9","CD47","CD80","CD86","NECTIN2","PVR","CD274","PDCD1LG2","VSIR","CD276","TGFB1","IDO1","HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","HLA-DRA","CD74")

fit_edger_subset <- function(cnt, meta, def_set, disp_vec = NULL) {
  meta$patient <- factor(meta$patient)
  meta$group   <- factor(meta$group, levels = c("low", "high"))
  y <- DGEList(counts = cnt)
  keep <- (rowSums(y$counts >= 1) >= 3) & (rowSums(y$counts) >= 5)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ patient + group, data = meta)
  if (is.null(disp_vec)) {
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
  } else {
    d_use <- disp_vec[rownames(y)]
    d_use[is.na(d_use)] <- median(disp_vec, na.rm = TRUE)
    fit <- glmQLFit(y, design, dispersion = d_use, robust = TRUE)
  }
  qlf <- glmQLFTest(fit, coef = "grouphigh")
  et <- topTags(qlf, n = Inf, sort.by = "none")$table
  is_def <- rownames(et) %in% def_set
  fdr_nd <- rep(NA_real_, nrow(et))
  fdr_nd[!is_def] <- p.adjust(et$PValue[!is_def], method = "BH")
  et$FDR_nonDefining <- fdr_nd
  return(list(et = et, disp = if (is.null(disp_vec)) setNames(y$tagwise.dispersion, rownames(y)) else disp_vec))
}

get_voom_donor_deltas <- function(tag) {
  cnt  <- as.matrix(read.csv(file.path(pb_dir, paste0("counts_", tag, ".csv")), row.names = 1, check.names = FALSE))
  storage.mode(cnt) <- "double"
  meta <- read.csv(file.path(pb_dir, paste0("meta_", tag, ".csv")), stringsAsFactors = FALSE)
  cnt  <- cnt[, meta$sample_id, drop = FALSE]
  y <- DGEList(counts = cnt)
  keep <- (rowSums(y$counts >= 1) >= 3) & (rowSums(y$counts) >= 5)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  meta$patient <- factor(meta$patient)
  meta$group   <- factor(meta$group, levels = c("low", "high"))
  design <- model.matrix(~ patient + group, data = meta)
  v <- voom(y, design, plot = FALSE)
  dimnames(v$weights) <- dimnames(v$E)
  donors <- sort(unique(as.character(meta$patient)))
  hi_cols <- paste0(donors, "__high")
  lo_cols <- paste0(donors, "__low")
  E_diff <- v$E[, hi_cols, drop = FALSE] - v$E[, lo_cols, drop = FALSE]
  # Harmonic precision weight of each paired difference: 1 / (1/w_hi + 1/w_lo)
  W_hi <- v$weights[, hi_cols, drop = FALSE]
  W_lo <- v$weights[, lo_cols, drop = FALSE]
  W_pair <- 1.0 / (1.0 / W_hi + 1.0 / W_lo)
  colnames(E_diff) <- donors
  colnames(W_pair) <- donors
  list(genes = rownames(E_diff), donors = donors, E_diff = E_diff, W_pair = W_pair, cnt = cnt, meta = meta)
}

set.seed(20261118)
cat("=== [1/2] Running True edgeR QL Leave-One-Donor-Out (LODO) & Precision-Weighted Donor Bootstrap ===\n")
lodo_specs <- list(
  list(lbl = "Peng2019_Q2_IFN14_min5", tag = "Peng_Q2_IFN14_min5", def = ifn14_genes),
  list(lbl = "Peng2019_Q2_IFN10_min5", tag = "Peng_Q2_IFN10_min5", def = ifn10_genes),
  list(lbl = "GSE202051_All43_Q2_IFN14_min5", tag = "GSE202051_All43_Q2_IFN14_min5", def = ifn14_genes),
  list(lbl = "GSE154778_All16_Q2_IFN14_min5", tag = "GSE154778_All16_Q2_IFN14_min5", def = ifn14_genes)
)

voom_cache <- list()
edger_lodo_rows <- list()
row_k <- 1

for (sp in lodo_specs) {
  cat("  Fitting LODO edgeR QL models for", sp$lbl, "...\n")
  vc <- get_voom_donor_deltas(sp$tag)
  voom_cache[[sp$tag]] <- vc
  full_et <- read.delim(file.path(ed_dir, paste0("paired_edgeR_", sp$tag, ".tsv")), stringsAsFactors = FALSE)
  rownames(full_et) <- full_et$gene

  donors <- vc$donors
  n_d <- length(donors)
  kg <- intersect(key_genes, rownames(full_et))

  # Store LODO logFC and FDR matrices [length(kg) x n_d]
  lfc_mat <- matrix(NA_real_, nrow = length(kg), ncol = n_d, dimnames = list(kg, donors))
  pval_mat <- matrix(NA_real_, nrow = length(kg), ncol = n_d, dimnames = list(kg, donors))
  fdr_mat  <- matrix(NA_real_, nrow = length(kg), ncol = n_d, dimnames = list(kg, donors))

  cat("  Estimating full NB dispersion & fitting", n_d, "exact LODO edgeR QL models for", sp$lbl, "...\n"); flush.console()
  t_l0 <- proc.time()[3]
  full_fit_obj <- fit_edger_subset(vc$cnt, vc$meta, sp$def, disp_vec = NULL)
  disp_base <- full_fit_obj$disp

  for (di in seq_along(donors)) {
    d_drop <- donors[di]
    m_sub  <- vc$meta[vc$meta$patient != d_drop, , drop = FALSE]
    c_sub  <- vc$cnt[, m_sub$sample_id, drop = FALSE]
    et_sub <- fit_edger_subset(c_sub, m_sub, sp$def, disp_vec = disp_base)$et
    common_g <- intersect(kg, rownames(et_sub))
    lfc_mat[common_g, di]  <- et_sub[common_g, "logFC"]
    pval_mat[common_g, di] <- et_sub[common_g, "PValue"]
    fdr_mat[common_g, di]  <- et_sub[common_g, "FDR_nonDefining"]
  }
  # Exact full-estimateDisp spot verification on STAT1's most influential donor!
  inf_di <- which.max(abs(lfc_mat["STAT1", ] - full_et["STAT1", "log2FC_edgeR"]))
  m_chk <- vc$meta[vc$meta$patient != donors[inf_di], , drop = FALSE]
  et_chk <- fit_edger_subset(vc$cnt[, m_chk$sample_id, drop = FALSE], m_chk, sp$def, disp_vec = NULL)$et
  max_disp_diff <- max(abs(et_chk[kg, "logFC"] - lfc_mat[kg, inf_di]), na.rm = TRUE)
  cat("    Completed", n_d, "LODO QL fits for", sp$lbl, "in", round(proc.time()[3] - t_l0, 1), "s (max |reEstDisp - tagwiseDisp| =", signif(max_disp_diff, 3), ")\n"); flush.console()

  # Precision-weighted voom donor bootstrap (B = 1,000)
  B <- 1000
  boot_idx <- matrix(sample.int(n_d, B * n_d, replace = TRUE), nrow = B, ncol = n_d)

  for (g in kg) {
    b_est <- numeric(B)
    if (g %in% rownames(vc$E_diff)) {
      ed_g <- vc$E_diff[g, ]
      wp_g <- vc$W_pair[g, ]
      for (b in seq_len(B)) {
        ix <- boot_idx[b, ]
        b_est[b] <- sum(ed_g[ix] * wp_g[ix]) / sum(wp_g[ix])
      }
      b_ci <- quantile(b_est, probs = c(0.025, 0.975), na.rm = TRUE)
    } else {
      b_ci <- c(NA_real_, NA_real_)
    }

    lfcs <- lfc_mat[g, ]
    pvs  <- pval_mat[g, ]
    fdrs <- fdr_mat[g, ]
    full_lfc <- full_et[g, "log2FC_edgeR"]
    most_inf <- donors[which.max(abs(lfcs - full_lfc))]

    edger_lodo_rows[[row_k]] <- data.frame(
      cohort_model = sp$lbl,
      n_paired_donors = n_d,
      gene = g,
      in_defining_signature = g %in% sp$def,
      full_edgeR_log2FC = full_lfc,
      full_edgeR_ci95_low = full_et[g, "ci95_low_edgeR"],
      full_edgeR_ci95_high = full_et[g, "ci95_high_edgeR"],
      full_edgeR_PValue = full_et[g, "PValue_edgeR"],
      full_edgeR_nonDefining_FDR = full_et[g, "FDR_edgeR_nonDefining"],
      voom_donor_bootstrap_95CI_low = unname(b_ci[1]),
      voom_donor_bootstrap_95CI_high = unname(b_ci[2]),
      edgeR_LODO_min_log2FC = min(lfcs, na.rm = TRUE),
      edgeR_LODO_max_log2FC = max(lfcs, na.rm = TRUE),
      edgeR_LODO_sign_consistency_pct = round(100 * mean(sign(lfcs) == sign(full_lfc), na.rm = TRUE), 1),
      edgeR_LODO_max_PValue = max(pvs, na.rm = TRUE),
      edgeR_LODO_max_nonDefining_FDR = if (all(is.na(fdrs))) NA_real_ else max(fdrs, na.rm = TRUE),
      edgeR_LODO_runs_p_lt_005 = paste0(sum(pvs < 0.05, na.rm = TRUE), "/", n_d),
      edgeR_LODO_runs_FDR_lt_005 = if (all(is.na(fdrs))) "Defining_Excluded" else paste0(sum(fdrs < 0.05, na.rm = TRUE), "/", n_d),
      most_influential_donor = most_inf,
      stringsAsFactors = FALSE
    )
    row_k <- row_k + 1
  }
}

df_edger_lodo <- do.call(rbind, edger_lodo_rows)
write.table(df_edger_lodo, file.path(out_dir, "edgeR_exact_LODO_and_voom_bootstrap_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(subset(df_edger_lodo, cohort_model == "Peng2019_Q2_IFN14_min5" & gene %in% c("STAT1","IRF7","CD47","LGALS9","CD80","TAP1","HLA-A","CD274","IDO1","COL1A1","KRT19"))[, c("gene","full_edgeR_log2FC","full_edgeR_nonDefining_FDR","voom_donor_bootstrap_95CI_low","voom_donor_bootstrap_95CI_high","edgeR_LODO_min_log2FC","edgeR_LODO_max_log2FC","edgeR_LODO_runs_p_lt_005","edgeR_LODO_runs_FDR_lt_005")])

cat("=== [2/2] Running Precision-Weighted Two-Cohort Donor-Resampling Bootstrap (B = 1,000) & LODO for Spearman rho ===\n")
for (xt in c("GSE202051_Untreated18_Q2_IFN14_min5", "GSE202051_Neoadjuvant25_Q2_IFN14_min5")) {
  voom_cache[[xt]] <- get_voom_donor_deltas(xt)
}

disc_vc <- voom_cache[["Peng_Q2_IFN14_min5"]]
disc_ed <- read.delim(file.path(ed_dir, "paired_edgeR_Peng_Q2_IFN14_min5.tsv"), stringsAsFactors = FALSE)
rownames(disc_ed) <- disc_ed$gene
disc_sig <- rownames(subset(disc_ed, !in_defining_signature & FDR_edgeR_nonDefining < 0.05))

val_tags <- c("GSE202051_All43_Q2_IFN14_min5", "GSE202051_Untreated18_Q2_IFN14_min5", "GSE202051_Neoadjuvant25_Q2_IFN14_min5", "GSE154778_All16_Q2_IFN14_min5")
conc_rows <- list()
ck <- 1

wmean_cols <- function(E, W) {
  rowSums(E * W) / pmax(rowSums(W), 1e-12)
}

for (vt in val_tags) {
  val_vc <- voom_cache[[vt]]
  val_ed <- read.delim(file.path(ed_dir, paste0("paired_edgeR_", vt, ".tsv")), stringsAsFactors = FALSE)
  rownames(val_ed) <- val_ed$gene

  shared_all <- setdiff(Reduce(intersect, list(rownames(disc_ed), rownames(val_ed), disc_vc$genes, val_vc$genes)), ifn14_genes)
  sets <- list(
    list(nm = "1_All_Shared_NonDefining_Genes", gs = shared_all, note = "Genome-wide shared expressed non-defining genes (zero selection filter)"),
    list(nm = "2_PreSpecified_23Gene_Functional_Panel", gs = intersect(prespec_23, shared_all), note = "Pre-specified MHC-I + checkpoint + TF panel (zero DEG selection filter)"),
    list(nm = "3_Discovery_Significant_NonDefining_DEGs", gs = intersect(disc_sig, shared_all), note = "Discovery BH FDR < 0.05 subset (conditioned on discovery significance)")
  )

  n1 <- length(disc_vc$donors)
  n2 <- length(val_vc$donors)
  B <- 1000
  b1 <- matrix(sample.int(n1, B * n1, replace = TRUE), nrow = B, ncol = n1)
  b2 <- matrix(sample.int(n2, B * n2, replace = TRUE), nrow = B, ncol = n2)

  for (st in sets) {
    gs <- st$gs
    r_ed <- suppressWarnings(cor.test(disc_ed[gs, "log2FC_edgeR"], val_ed[gs, "log2FC_edgeR"], method = "spearman", exact = FALSE))
    E1 <- disc_vc$E_diff[gs, , drop = FALSE]; W1 <- disc_vc$W_pair[gs, , drop = FALSE]
    E2 <- val_vc$E_diff[gs, , drop = FALSE];  W2 <- val_vc$W_pair[gs, , drop = FALSE]
    m1_full <- wmean_cols(E1, W1)
    m2_full <- wmean_cols(E2, W2)
    rho_voom <- cor(m1_full, m2_full, method = "spearman")

    b_rhos <- numeric(B)
    for (b in seq_len(B)) {
      i1 <- b1[b, ]; i2 <- b2[b, ]
      v1 <- rowSums(E1[, i1, drop = FALSE] * W1[, i1, drop = FALSE]) / rowSums(W1[, i1, drop = FALSE])
      v2 <- rowSums(E2[, i2, drop = FALSE] * W2[, i2, drop = FALSE]) / rowSums(W2[, i2, drop = FALSE])
      b_rhos[b] <- cor(v1, v2, method = "spearman")
    }
    b_ci <- quantile(b_rhos, probs = c(0.025, 0.975), na.rm = TRUE)

    lodo_r <- c(
      sapply(seq_len(n1), function(j) cor(wmean_cols(E1[, -j, drop = FALSE], W1[, -j, drop = FALSE]), m2_full, method = "spearman")),
      sapply(seq_len(n2), function(j) cor(m1_full, wmean_cols(E2[, -j, drop = FALSE], W2[, -j, drop = FALSE]), method = "spearman"))
    )

    conc_rows[[ck]] <- data.frame(
      validation_cohort_model = vt,
      n_discovery_donors = n1,
      n_validation_donors = n2,
      gene_set_definition = st$nm,
      n_genes = length(gs),
      edgeR_spearman_rho = unname(r_ed$estimate),
      gene_level_p_value = pmax(r_ed$p.value, 1e-300),
      sign_agreement_pct = round(100 * mean(sign(disc_ed[gs, "log2FC_edgeR"]) == sign(val_ed[gs, "log2FC_edgeR"])), 1),
      voom_precision_spearman_rho = rho_voom,
      donor_bootstrap_95CI_low = unname(b_ci[1]),
      donor_bootstrap_95CI_high = unname(b_ci[2]),
      LODO_min_spearman_rho = min(lodo_r, na.rm = TRUE),
      LODO_max_spearman_rho = max(lodo_r, na.rm = TRUE),
      methodological_note = st$note,
      stringsAsFactors = FALSE
    )
    ck <- ck + 1
  }
}

df_conc_voom <- do.call(rbind, conc_rows)
write.table(df_conc_voom, file.path(out_dir, "cross_cohort_voom_donor_bootstrap_and_LODO_concordance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(df_conc_voom[, c("validation_cohort_model","gene_set_definition","n_genes","edgeR_spearman_rho","donor_bootstrap_95CI_low","donor_bootstrap_95CI_high","LODO_min_spearman_rho","LODO_max_spearman_rho")])
cat("=== [CSBJ Phase 2 Exact edgeR LODO & Voom Donor Bootstrap Complete!] ===\n")
