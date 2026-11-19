#!/usr/bin/env Rscript
# Stage 6 Part A (Gate G6): Bounded Clinical & Compositional Survival Analysis in TCGA-PAAD
# - Primary cohort: Histologically confirmed Pancreas Adenocarcinoma, Ductal Type (n = 146; 145 OS valid; 140 complete-case multivariable)
# - Sensitivity cohort: All primary TCGA-PAAD tumors (n = 177; 176 OS valid) demonstrating histologic confounding of STAT1
# - Endpoints explicitly defined: OS (Overall Survival), PFI (Progression-Free Interval, cBioPortal PFS_MONTHS), DFI (Disease-Free Interval, DFS_MONTHS)
# - Three pre-specified models:
#     Model 1: Univariate continuous per-1-SD log2(RSEM+1)
#     Model 2 (Primary): Clinicopathologic multivariable (Age + Sex + AJCC Stage + Histologic Grade)
#     Model 3 (Sensitivity): Purity & immune-composition adjusted (Model 2 + Epithelial purity proxy [EPCAM/KRT19/CDH1] + Macrophage abundance [CD68/CD163/CSF1R/C1QA] + PTPRC)
# - Proportional hazards (Schoenfeld cox.zph) + Restricted natural spline nonlinearity LRT (df=3 vs linear)

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

out_root <- "/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116"
v1_root  <- "/mnt/d/Research/2026_PDAC_RegVelo"
tcga_dir <- file.path(v1_root, "02_data/tcga")
st6_dir  <- file.path(out_root, "05_results/stage6_clinical_drug")
dir.create(st6_dir, recursive = TRUE, showWarnings = FALSE)

clin_pat <- read.csv(file.path(tcga_dir, "clinical_patient.csv"), stringsAsFactors = FALSE, check.names = FALSE)
clin_smp <- read.csv(file.path(tcga_dir, "clinical_sample.csv"), stringsAsFactors = FALSE, check.names = FALSE)
expr_raw <- read.csv(file.path(tcga_dir, "mrna_rsem_log_not_applied.csv"), stringsAsFactors = FALSE, check.names = FALSE)

genes <- expr_raw[[1]]
mat <- as.matrix(expr_raw[, -1, drop = FALSE])
storage.mode(mat) <- "double"
if (any(duplicated(genes))) {
  means <- rowMeans(mat, na.rm = TRUE)
  ord <- order(means, decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; genes <- genes[ord]
  keep <- !duplicated(genes)
  mat <- mat[keep, , drop = FALSE]; genes <- genes[keep]
}
rownames(mat) <- genes
expr_log2 <- log2(pmax(mat, 0) + 1)

clin_smp$patientId <- substr(clin_smp$sampleId, 1, 12)
meta <- merge(clin_smp, clin_pat, by = "patientId", all.x = TRUE)
meta <- meta[meta$sampleId %in% colnames(expr_log2) & meta$SAMPLE_TYPE == "Primary", ]
meta <- meta[!duplicated(meta$patientId), ]
rownames(meta) <- meta$sampleId

parse_ev <- function(x) ifelse(is.na(x) | x == "", NA_real_, as.numeric(substr(trimws(x), 1, 1)))
meta$os_months  <- suppressWarnings(as.numeric(meta$OS_MONTHS))
meta$os_event   <- parse_ev(meta$OS_STATUS)
meta$pfi_months <- suppressWarnings(as.numeric(meta$PFS_MONTHS))
meta$pfi_event  <- parse_ev(meta$PFS_STATUS)
meta$dfi_months <- suppressWarnings(as.numeric(meta$DFS_MONTHS))
meta$dfi_event  <- parse_ev(meta$DFS_STATUS)

meta$age <- suppressWarnings(as.numeric(meta$AGE))
meta$sex_male <- ifelse(meta$SEX == "Male", 1, ifelse(meta$SEX == "Female", 0, NA_real_))
st <- toupper(trimws(meta$AJCC_PATHOLOGIC_TUMOR_STAGE))
meta$stage_ord <- ifelse(grepl("^STAGE I($|A|B)", st), 1,
                  ifelse(grepl("^STAGE II($|A|B)", st), 2,
                  ifelse(grepl("^STAGE (III|IV)", st), 3, NA_real_)))
gr <- toupper(trimws(meta$GRADE))
meta$grade_ord <- ifelse(gr == "G1", 1, ifelse(gr == "G2", 2, ifelse(gr %in% c("G3","G4"), 3, NA_real_)))
meta$is_ductal <- meta$TUMOR_TYPE == "Pancreas Adenocarcinoma, Ductal Type"

# Compute composite purity & immune covariates across samples
epi_genes <- intersect(c("EPCAM", "KRT19", "CDH1", "KRT7"), rownames(expr_log2))
mac_genes <- intersect(c("CD68", "CD163", "CSF1R", "C1QA", "AIF1"), rownames(expr_log2))
ifn14_genes <- intersect(c("IRF1","CXCL9","CXCL10","CXCL11","MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2"), rownames(expr_log2))
ifn10_genes <- intersect(c("MX1","MX2","IFIT1","IFIT3","GBP1","GBP5","SOCS1","SOCS3","JAK1","JAK2"), rownames(expr_log2))

z_by_row <- function(m) t(scale(t(m)))

# Official Lambert 2018 verified human sequence-specific TFs ONLY!
lambert_f <- file.path(v1_root, "02_data/reference/prior_networks/lambert2018_human_tfs.txt")
lambert_tfs <- unique(trimws(readLines(lambert_f)))
lambert_tfs <- lambert_tfs[lambert_tfs != "" & lambert_tfs %in% rownames(expr_log2)]
cat(sprintf("Verified Lambert 2018 human TFs expressed in TCGA-PAAD: %d\n", length(lambert_tfs)))

# Reverse-KM median follow-up helper
median_followup <- function(time_v, ev_v) {
  sf <- survfit(Surv(time_v, 1 - ev_v) ~ 1)
  as.numeric(summary(sf)$table["median"])
}

# Missingness & endpoint accounting table
acct_rows <- list()
for (coh_lbl in c("Pure_Ductal_PDAC_n146", "Unfiltered_All_PAAD_n177")) {
  ids <- if (coh_lbl == "Pure_Ductal_PDAC_n146") rownames(meta)[meta$is_ductal] else rownames(meta)
  sm <- meta[ids, ]
  for (ep in c("os", "pfi", "dfi")) {
    tv <- sm[[paste0(ep, "_months")]]
    ev <- sm[[paste0(ep, "_event")]]
    ok_u <- !is.na(tv) & tv > 0 & !is.na(ev)
    ok_m <- ok_u & !is.na(sm$age) & !is.na(sm$sex_male) & !is.na(sm$stage_ord) & !is.na(sm$grade_ord)
    acct_rows[[paste(coh_lbl, ep, sep = "__")]] <- data.frame(
      cohort = coh_lbl,
      endpoint = toupper(ep),
      n_cohort_total = length(ids),
      n_valid_univariate = sum(ok_u),
      n_events_univariate = sum(ev[ok_u]),
      n_censored_univariate = sum(1 - ev[ok_u]),
      median_followup_months_revKM = median_followup(tv[ok_u], ev[ok_u]),
      n_complete_multivariable = sum(ok_m),
      n_events_multivariable = sum(ev[ok_m]),
      n_missing_covariates = sum(ok_u) - sum(ok_m),
      stringsAsFactors = FALSE
    )
  }
}
df_acct <- do.call(rbind, acct_rows)
write.table(df_acct, file.path(st6_dir, "tcga_endpoint_and_missingness_accounting.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(df_acct, row.names = FALSE)

# Function to fit all 3 models + PH + spline nonlinearity for candidate features
key_candidates <- c("STAT1", "IRF1", "IFN14_Signature", "IFN10_Signature", "HMGA2", "FOXM1", "KLF5", "GATA6", "E2F1", "TFAP2C", "SATB1", "SPI1", "IRF8", "NFKB1")

run_full_survival <- function(ids, cohort_tag, feature_vec, compute_all_tfs = FALSE) {
  sm <- meta[ids, , drop = FALSE]
  ex <- expr_log2[, ids, drop = FALSE]

  epi_sc <- colMeans(z_by_row(ex[epi_genes, , drop = FALSE]), na.rm = TRUE)
  mac_sc <- colMeans(z_by_row(ex[mac_genes, , drop = FALSE]), na.rm = TRUE)
  ptprc_sc <- as.numeric(scale(as.numeric(ex["PTPRC", ])))
  ifn14_sc <- colMeans(z_by_row(ex[ifn14_genes, , drop = FALSE]), na.rm = TRUE)
  ifn10_sc <- colMeans(z_by_row(ex[ifn10_genes, , drop = FALSE]), na.rm = TRUE)

  feat_mat <- rbind(
    ex[intersect(feature_vec, rownames(ex)), , drop = FALSE],
    IFN14_Signature = ifn14_sc,
    IFN10_Signature = ifn10_sc
  )
  target_names <- unique(c(intersect(key_candidates, rownames(feat_mat)), if (compute_all_tfs) intersect(lambert_tfs, rownames(feat_mat)) else character(0)))

  out_rows <- list()
  for (ep in c("os", "pfi", "dfi")) {
    if (ep == "dfi" && compute_all_tfs) next  # OS & PFI primary across all TFs; DFI included for key candidates
    tv <- sm[[paste0(ep, "_months")]]
    ev <- sm[[paste0(ep, "_event")]]
    ok_u <- !is.na(tv) & tv > 0 & !is.na(ev)

    ep_rows <- vector("list", length(target_names))
    for (i in seq_along(target_names)) {
      fn <- target_names[i]
      raw_x <- as.numeric(feat_mat[fn, ])
      if (sd(raw_x[ok_u], na.rm = TRUE) < 1e-6) next

      df_u <- data.frame(
        time = tv[ok_u], event = ev[ok_u],
        x_z = as.numeric(scale(raw_x[ok_u])),
        age_z = as.numeric(scale(sm$age[ok_u])),
        sex_male = sm$sex_male[ok_u],
        stage_ord = sm$stage_ord[ok_u],
        grade_ord = sm$grade_ord[ok_u],
        epi_z = as.numeric(scale(epi_sc[ok_u])),
        mac_z = as.numeric(scale(mac_sc[ok_u])),
        ptprc_z = as.numeric(scale(ptprc_sc[ok_u]))
      )

      # Model 1: Univariate
      f1 <- tryCatch(coxph(Surv(time, event) ~ x_z, data = df_u), error = function(e) NULL)
      if (is.null(f1)) next
      s1 <- summary(f1)

      # Model 2: Primary Multivariable (Age + Sex + Stage + Grade)
      df_m <- df_u[complete.cases(df_u[, c("time","event","x_z","age_z","sex_male","stage_ord","grade_ord")]), ]
      df_m$x_z <- as.numeric(scale(df_m$x_z))
      f2 <- tryCatch(coxph(Surv(time, event) ~ x_z + age_z + sex_male + stage_ord + grade_ord, data = df_m), error = function(e) NULL)
      s2 <- if (!is.null(f2)) summary(f2) else NULL
      ph2 <- if (!is.null(f2)) tryCatch(cox.zph(f2)$table["x_z", "p"], error = function(e) NA_real_) else NA_real_

      # Nonlinearity check (natural spline df=3 vs linear in Model 2) for key candidates
      nonlin_p <- NA_real_
      if (!is.null(f2) && (fn %in% key_candidates)) {
        f2_sp <- tryCatch(coxph(Surv(time, event) ~ ns(x_z, df = 3) + age_z + sex_male + stage_ord + grade_ord, data = df_m), error = function(e) NULL)
        if (!is.null(f2_sp)) {
          lrt_sp <- anova(f2, f2_sp)
          nonlin_p <- lrt_sp$`Pr(>|Chi|)`[2]
        }
      }

      # Model 3: Purity + Macrophage + Leukocyte Composition-Adjusted Sensitivity
      f3 <- if (fn %in% key_candidates) {
        tryCatch(coxph(Surv(time, event) ~ x_z + age_z + sex_male + stage_ord + grade_ord + epi_z + mac_z + ptprc_z, data = df_m), error = function(e) NULL)
      } else NULL
      s3 <- if (!is.null(f3)) summary(f3) else NULL

      ep_rows[[i]] <- data.frame(
        cohort = cohort_tag,
        endpoint = toupper(ep),
        feature = fn,
        is_lambert_tf = fn %in% lambert_tfs,
        expression_scale = "per_1_SD_log2_RSEM_plus_1",
        n_uni = nrow(df_u),
        events_uni = sum(df_u$event),
        HR_uni = s1$conf.int["x_z", "exp(coef)"],
        CI95_low_uni = s1$conf.int["x_z", "lower .95"],
        CI95_high_uni = s1$conf.int["x_z", "upper .95"],
        p_uni = pmax(s1$coefficients["x_z", "Pr(>|z|)"], 1e-300),
        n_multi = nrow(df_m),
        events_multi = sum(df_m$event),
        covariates_multi = "Age + Sex + AJCC_Stage + Histologic_Grade",
        HR_multi = if (!is.null(s2)) s2$conf.int["x_z", "exp(coef)"] else NA_real_,
        CI95_low_multi = if (!is.null(s2)) s2$conf.int["x_z", "lower .95"] else NA_real_,
        CI95_high_multi = if (!is.null(s2)) s2$conf.int["x_z", "upper .95"] else NA_real_,
        p_multi = if (!is.null(s2)) pmax(s2$coefficients["x_z", "Pr(>|z|)"], 1e-300) else NA_real_,
        schoenfeld_ph_p_multi = ph2,
        spline_nonlinearity_p = nonlin_p,
        HR_purity_immune_adj = if (!is.null(s3)) s3$conf.int["x_z", "exp(coef)"] else NA_real_,
        CI95_low_purity_immune = if (!is.null(s3)) s3$conf.int["x_z", "lower .95"] else NA_real_,
        CI95_high_purity_immune = if (!is.null(s3)) s3$conf.int["x_z", "upper .95"] else NA_real_,
        p_purity_immune_adj = if (!is.null(s3)) pmax(s3$coefficients["x_z", "Pr(>|z|)"], 1e-300) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    ep_df <- do.call(rbind, ep_rows)
    tf_m <- ep_df$is_lambert_tf
    ep_df$FDR_uni_allLambertTFs <- NA_real_
    ep_df$FDR_multi_allLambertTFs <- NA_real_
    if (sum(tf_m) > 10) {
      ep_df$FDR_uni_allLambertTFs[tf_m] <- p.adjust(ep_df$p_uni[tf_m], method = "BH")
      ep_df$FDR_multi_allLambertTFs[tf_m] <- p.adjust(ep_df$p_multi[tf_m], method = "BH")
    }
    out_rows[[ep]] <- ep_df
  }
  do.call(rbind, out_rows)
}

cat("\nRunning Pure Ductal PDAC (n=146) + Unfiltered All PAAD (n=177) Cox models...\n")
res_ductal <- run_full_survival(rownames(meta)[meta$is_ductal], "Pure_Ductal_PDAC_n146", c(key_candidates, lambert_tfs), compute_all_tfs = TRUE)
res_allpaad <- run_full_survival(rownames(meta), "Unfiltered_All_PAAD_n177", key_candidates, compute_all_tfs = FALSE)

surv_frozen <- rbind(res_ductal, res_allpaad)
write.table(surv_frozen, file.path(out_root, "05_results/survival_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(surv_frozen[surv_frozen$feature %in% key_candidates, ], file.path(st6_dir, "tcga_key_candidates_survival_3models.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Key Candidates in Pure Ductal PDAC (OS) ===\n")
sub_print <- surv_frozen[surv_frozen$cohort == "Pure_Ductal_PDAC_n146" & surv_frozen$endpoint == "OS" & surv_frozen$feature %in% key_candidates,
                         c("feature", "HR_uni", "p_uni", "HR_multi", "p_multi", "schoenfeld_ph_p_multi", "spline_nonlinearity_p", "HR_purity_immune_adj", "p_purity_immune_adj")]
print(sub_print[order(sub_print$p_multi), ], row.names = FALSE)
cat("=== [Stage 6 Part A TCGA Complete] ===\n")
