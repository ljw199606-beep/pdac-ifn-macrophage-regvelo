#!/usr/bin/env Rscript
# Stage 3 Part B (Gate G3):
# 1. Paired donor-level edgeR TMM QL NB GLM (~ patient + group) + limma-voom 95% CIs across all 12 primary & sensitivity tags
# 2. Ambient extracellular RNA soup cancellation table (unpaired Tumor-vs-Normal vs within-donor paired)
# 3. Complete unselected transcriptome-wide & non-IFN cross-cohort concordance audit (Peng vs GSE202051 Untreated/All/Neoadj & GSE154778)
# 4. Non-circular MSigDB Hallmark v2024.1.Hs fgsea across ALL 50 pathways (after removing IFN-14 defining genes from ranked universe)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
})

out_root <- "/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116"
v1_root  <- "/mnt/d/Research/2026_PDAC_RegVelo"
st3_dir  <- file.path(out_root, "05_results/stage3_macrophage")
pb_dir   <- file.path(st3_dir, "pseudobulk_inputs")
edger_dir <- file.path(st3_dir, "edger_outputs")
dir.create(edger_dir, recursive = TRUE, showWarnings = FALSE)

ifn15_genes <- c("STAT1","IRF1","CXCL9","CXCL10","CXCL11","MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2")
ifn14_genes <- setdiff(ifn15_genes, "STAT1")
ifn10_genes <- c("MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2")

tags <- c(
  "Peng_Q2_IFN14_min5",
  "Peng_PDAConly_Q2_IFN14_min5",
  "Peng_Q2_IFN14_min3",
  "Peng_Q2_IFN14_min10",
  "Peng_Q1_IFN14_min5",
  "Peng_Q2_IFN10_min5",
  "GSE202051_All43_Q2_IFN14_min5",
  "GSE202051_Untreated18_Q2_IFN14_min5",
  "GSE202051_Neoadjuvant25_Q2_IFN14_min5",
  "GSE202051_All43_Q1_IFN14_min5",
  "GSE202051_All43_Q2_IFN10_min5",
  "GSE154778_All16_Q2_IFN14_min5"
)

cat("=== [1/4] Fitting Paired Donor edgeR QL + limma-voom 95% CIs ===\n")
model_summary <- list()
de_tables <- list()

for (tg in tags) {
  cnt_f  <- file.path(pb_dir, paste0("counts_", tg, ".csv"))
  meta_f <- file.path(pb_dir, paste0("meta_", tg, ".csv"))
  if (!file.exists(cnt_f)) next

  cnt  <- as.matrix(read.csv(cnt_f, row.names = 1, check.names = FALSE))
  storage.mode(cnt) <- "double"
  meta <- read.csv(meta_f, stringsAsFactors = FALSE)
  cnt  <- cnt[, meta$sample_id, drop = FALSE]

  meta$patient <- factor(meta$patient)
  meta$group   <- factor(meta$group, levels = c("low", "high"))

  y <- DGEList(counts = cnt)
  # Retain genes with >=1 UMI in at least 3 pseudobulk samples & total >=5 UMIs so sparse immune checkpoints are tested without exclusion
  keep_relaxed <- (rowSums(y$counts >= 1) >= 3) & (rowSums(y$counts) >= 5)
  keep_strict  <- filterByExpr(y, group = meta$group)
  y <- y[keep_relaxed, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ patient + group, data = meta)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "grouphigh")
  et <- topTags(qlf, n = Inf, sort.by = "none")$table

  # Signed quasi-likelihood t-like statistic + approximate 95% CI on log2FC from F-statistic!
  # Since F = t^2 with df2 = qlf$df.total, |SE| = |logFC| / sqrt(max(F, 1e-12))
  df_tot <- qlf$df.total
  t_crit <- qt(0.975, df = pmax(df_tot, 1))
  se_ql  <- ifelse(et$F > 1e-10, abs(et$logFC) / sqrt(et$F), NA_real_)
  ci_lo  <- et$logFC - t_crit * se_ql
  ci_hi  <- et$logFC + t_crit * se_ql
  signed_t <- sign(et$logFC) * sqrt(pmax(et$F, 0))

  # Also fit limma-voom paired model for sensitivity check
  v <- voom(y, design, plot = FALSE)
  lfit <- eBayes(lmFit(v, design), robust = TRUE)
  ltab <- topTable(lfit, coef = "grouphigh", number = Inf, sort.by = "none", confint = TRUE)
  m_idx <- match(rownames(et), rownames(ltab))

  def_set <- if (grepl("IFN10", tg)) ifn10_genes else ifn14_genes
  is_def  <- rownames(et) %in% def_set
  # Compute FDR strictly across non-defining genes as well!
  fdr_non_def <- rep(NA_real_, nrow(et))
  fdr_non_def[!is_def] <- p.adjust(et$PValue[!is_def], method = "BH")

  out_df <- data.frame(
    gene = rownames(et),
    in_defining_signature = is_def,
    passed_strict_filterByExpr = rownames(et) %in% names(keep_strict)[keep_strict],
    log2FC_edgeR = et$logFC,
    ci95_low_edgeR = ci_lo,
    ci95_high_edgeR = ci_hi,
    logCPM = et$logCPM,
    QL_F = et$F,
    signed_QL_t = signed_t,
    PValue_edgeR = pmax(et$PValue, 1e-300),
    FDR_edgeR_all = pmax(et$FDR, 1e-300),
    FDR_edgeR_nonDefining = pmax(fdr_non_def, 1e-300),
    log2FC_voom = ltab$logFC[m_idx],
    ci95_low_voom = ltab$CI.L[m_idx],
    ci95_high_voom = ltab$CI.R[m_idx],
    PValue_voom = pmax(ltab$P.Value[m_idx], 1e-300),
    FDR_voom = pmax(ltab$adj.P.Val[m_idx], 1e-300),
    stringsAsFactors = FALSE
  )
  out_df <- out_df[order(out_df$PValue_edgeR), ]
  write.table(out_df, file.path(edger_dir, paste0("paired_edgeR_", tg, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  de_tables[[tg]] <- out_df

  n_donors <- length(levels(meta$patient))
  n_sig_nondef <- sum(!out_df$in_defining_signature & out_df$FDR_edgeR_nonDefining < 0.05, na.rm = TRUE)
  n_sig_strict <- sum(!out_df$in_defining_signature & out_df$passed_strict_filterByExpr & out_df$FDR_edgeR_nonDefining < 0.05, na.rm = TRUE)
  model_summary[[tg]] <- data.frame(
    tag = tg,
    n_paired_donors = n_donors,
    n_pseudobulks = nrow(meta),
    n_genes_tested_relaxed = nrow(out_df),
    n_genes_strict_filterByExpr = sum(out_df$passed_strict_filterByExpr),
    n_sig_FDR005_nonDefining_relaxed = n_sig_nondef,
    n_sig_FDR005_nonDefining_strict = n_sig_strict,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s: n=%d paired donors | %d genes tested | %d non-defining DEGs (FDR<0.05)\n",
              tg, n_donors, nrow(out_df), n_sig_nondef))
}

df_mod_sum <- do.call(rbind, model_summary)
write.table(df_mod_sum, file.path(st3_dir, "paired_edger_models_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# =========================================================================
# [2/4] Ambient Extracellular RNA Soup Cancellation Diagnostic
# =========================================================================
cat("\n=== [2/4] Ambient RNA Soup Cancellation Proof (Unpaired Tumor-vs-Normal vs Within-Donor Paired) ===\n")
unp_mac <- read.csv(file.path(v1_root, "04_results/pseudobulk/deg_Macrophage_cell_edgeR_voom.csv"), stringsAsFactors = FALSE)
lfc_unp_col <- if ("log2FC" %in% names(unp_mac)) "log2FC" else "logFC_edgeR"
fdr_unp_col <- if ("FDR" %in% names(unp_mac)) "FDR" else "FDR_edgeR"

amb_genes <- c("COL1A1","COL1A2","COL3A1","ACTA2","KRT19","KRT7","PRSS1","CPA1","CELA3A",
               "STAT1","IRF1","LGALS9","CD47","CD80","CD274","B2M","TAP1","HLA-A")
p_main <- de_tables[["Peng_Q2_IFN14_min5"]]
amb_rows <- lapply(amb_genes, function(g) {
  ru <- unp_mac[unp_mac$gene == g, ]
  rp <- p_main[p_main$gene == g, ]
  data.frame(
    gene = g,
    marker_class = ifelse(g %in% c("COL1A1","COL1A2","COL3A1","ACTA2","KRT19","KRT7","PRSS1","CPA1","CELA3A"),
                          "Extracellular_Ambient_Soup_Marker", "True_Macrophage_IFN_Checkpoint_Target"),
    unpaired_TvN_log2FC = if (nrow(ru) > 0) ru[[lfc_unp_col]][1] else NA_real_,
    unpaired_TvN_FDR    = if (nrow(ru) > 0) ru[[fdr_unp_col]][1] else NA_real_,
    paired_withinDonor_log2FC = if (nrow(rp) > 0) rp$log2FC_edgeR[1] else NA_real_,
    paired_withinDonor_FDR    = if (nrow(rp) > 0) rp$FDR_edgeR_all[1] else NA_real_,
    stringsAsFactors = FALSE
  )
})
df_amb <- do.call(rbind, amb_rows)
write.table(df_amb, file.path(st3_dir, "ambient_rna_cancellation_diagnostic.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(df_amb, row.names = FALSE)

# =========================================================================
# [3/4] Complete Cross-Cohort Concordance Audit (Addressing C05/P1!)
# =========================================================================
cat("\n=== [3/4] Cross-Cohort Concordance Audit (All Shared Genes vs Non-IFN Sig vs v0 37-Gene Subset) ===\n")
v0_37_genes <- unique(c(ifn15_genes, "CD274","CD276","CD47","VSIR","IDO1","PDCD1LG2",
                        "CD68","CD163","C1QA","LYZ","AIF1","HLA-DRA","HLA-DRB1","TAP1",
                        "ISG15","CXCL16","CCL2","CCL3","CCL4","IL1B","TNF","NFKB1","RELA"))
prespec_non_ifn <- c("STAT1","LGALS9","CD47","CD80","CD86","NECTIN2","PVR","CD274","PDCD1LG2","VSIR","CD276","IDO1",
                     "HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","HLA-DRA","HLA-DRB1","CD74","TGFB1","ISG15","OAS1","OAS2","OAS3")

disc <- de_tables[["Peng_Q2_IFN14_min5"]]
val_targets <- c(
  "GSE202051_Untreated18_Q2_IFN14_min5",
  "GSE202051_All43_Q2_IFN14_min5",
  "GSE202051_Neoadjuvant25_Q2_IFN14_min5",
  "GSE154778_All16_Q2_IFN14_min5",
  "Peng_PDAConly_Q2_IFN14_min5"
)

conc_rows <- list()
for (vt in val_targets) {
  vdf <- de_tables[[vt]]
  if (is.null(vdf)) next
  m <- merge(disc, vdf, by = "gene", suffixes = c("_disc", "_val"))

  subsets <- list(
    "1_All_Shared_Expressed_Genes_NoIFN15" = !(m$gene %in% ifn15_genes),
    "2_Strict_Expressed_Genes_NoIFN15"     = !(m$gene %in% ifn15_genes) & m$passed_strict_filterByExpr_disc & m$passed_strict_filterByExpr_val,
    "3_Discovery_FDR005_Strictly_Excluding_All_IFN15" = !(m$gene %in% ifn15_genes) & (m$FDR_edgeR_nonDefining_disc < 0.05),
    "4_PreSpecified_26_NonIFN14_Checkpoint_MHC_ISG_Panel" = m$gene %in% prespec_non_ifn,
    "5_Audit_v0_37Genes_With_IFN15_Included" = m$gene %in% v0_37_genes,
    "6_Audit_v0_22Genes_After_Removing_IFN15" = (m$gene %in% v0_37_genes) & !(m$gene %in% ifn15_genes)
  )

  for (sname in names(subsets)) {
    sub <- m[subsets[[sname]], , drop = FALSE]
    if (nrow(sub) < 5) next
    ct <- suppressWarnings(cor.test(sub$log2FC_edgeR_disc, sub$log2FC_edgeR_val, method = "spearman", exact = FALSE))
    # Compute 95% Fisher-z approximate CI for Spearman rho
    r_val <- unname(ct$estimate)
    z_r   <- atanh(pmin(pmax(r_val, -0.9999), 0.9999))
    se_z  <- 1.03 / sqrt(max(nrow(sub) - 3, 1))
    r_lo  <- tanh(z_r - 1.96 * se_z)
    r_hi  <- tanh(z_r + 1.96 * se_z)
    # Exact two-sided p-value from t-distribution to avoid p=0 underflow!
    t_sp  <- abs(r_val) * sqrt((nrow(sub) - 2) / pmax(1 - r_val^2, 1e-15))
    p_safe <- pmax(2 * pt(t_sp, df = nrow(sub) - 2, lower.tail = FALSE), 1e-300)
    dir_conc <- mean(sign(sub$log2FC_edgeR_disc) == sign(sub$log2FC_edgeR_val)) * 100.0

    conc_rows[[paste(vt, sname, sep = "__")]] <- data.frame(
      discovery_model = "Peng_Q2_IFN14_min5",
      validation_model = vt,
      gene_set_definition = sname,
      n_genes = nrow(sub),
      spearman_rho = r_val,
      ci95_low = r_lo,
      ci95_high = r_hi,
      p_value = p_safe,
      pct_direction_concordant = dir_conc,
      stringsAsFactors = FALSE
    )
  }

  # Save merged gene-level discovery vs validation table for GSE202051_All43 & Untreated18
  if (vt %in% c("GSE202051_All43_Q2_IFN14_min5", "GSE202051_Untreated18_Q2_IFN14_min5", "GSE154778_All16_Q2_IFN14_min5")) {
    write.table(m, file.path(st3_dir, paste0("merged_disc_vs_", vt, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

df_conc <- do.call(rbind, conc_rows)
write.table(df_conc, file.path(st3_dir, "cross_cohort_concordance_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(df_conc[, c("validation_model", "gene_set_definition", "n_genes", "spearman_rho", "p_value", "pct_direction_concordant")], row.names = FALSE)

# =========================================================================
# [4/4] Non-Circular MSigDB Hallmark v2024.1.Hs fgsea (All 50 Pathways Reported!)
# =========================================================================
cat("\n=== [4/4] Running Non-Circular MSigDB Hallmark v2024.1.Hs fgsea (Excluding All IFN-14 Defining Genes!) ===\n")
gmt_file <- file.path(out_root, "02_metadata/h.all.v2024.1.Hs.symbols.gmt")
hallmark_pws <- gmtPathways(gmt_file)
# Add 2 explicit non-circular functional modules (disjoint from IFN-14)
hallmark_pws[["CUSTOM_ADAPTIVE_CHECKPOINT_LIGANDS_NO_IFN14"]] <- c("CD274","PDCD1LG2","CD47","LGALS9","NECTIN2","PVR","CD80","CD86","VSIR","CD276","IDO1","TNFRSF14","SIGLEC10","LILRB1","LILRB2")
hallmark_pws[["CUSTOM_ANTIGEN_PRESENTATION_MHC_NO_IFN14"]]    <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","B2M","TAP1","TAP2","TAPBP","CD74","CIITA","PSMB8","PSMB9","PSMB10")

fg_list <- list()
for (tg in c("Peng_Q2_IFN14_min5", "GSE202051_All43_Q2_IFN14_min5", "GSE202051_Untreated18_Q2_IFN14_min5", "GSE154778_All16_Q2_IFN14_min5")) {
  tdf <- de_tables[[tg]]
  if (is.null(tdf)) next
  # STRICT NON-CIRCULARITY: exclude all 14 IFN-14 defining genes from the ranked gene universe!
  sub <- tdf[!tdf$in_defining_signature & is.finite(tdf$signed_QL_t), ]
  stats_vec <- sub$signed_QL_t
  names(stats_vec) <- sub$gene
  stats_vec <- sort(stats_vec, decreasing = TRUE)

  set.seed(42)
  fg <- fgseaMultilevel(pathways = hallmark_pws, stats = stats_vec, minSize = 8, maxSize = 500)
  fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x, 10), collapse = ";"), character(1))
  fg <- as.data.frame(fg)
  fg$pval <- pmax(fg$pval, 1e-300)
  fg$padj <- pmax(fg$padj, 1e-300)
  fg$model <- tg
  fg$gene_set_source <- ifelse(grepl("^HALLMARK_", fg$pathway), "MSigDB_Hallmark_v2024.1.Hs", "Curated_NonCircular_Immunoregulation")
  fg$background_universe_size <- length(stats_vec)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  fg_list[[tg]] <- fg
}
df_fg_all <- do.call(rbind, fg_list)
write.table(df_fg_all, file.path(st3_dir, "noncircular_hallmark_fgsea_all52_pathways.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nTop 10 & Bottom 5 Non-Circular Hallmark Pathways in Peng_Q2_IFN14_min5:\n")
p_fg <- fg_list[["Peng_Q2_IFN14_min5"]]
p_fg_ord <- p_fg[order(-p_fg$NES), c("pathway", "size", "NES", "pval", "padj")]
print(rbind(head(p_fg_ord, 10), tail(p_fg_ord, 5)), row.names = FALSE)
cat("=== [Stage 3 R Part Complete] ===\n")
