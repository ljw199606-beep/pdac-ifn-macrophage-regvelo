#!/usr/bin/env python3
"""
Portable, clean-environment reproduction of core CSBJ manuscript figure panels directly
from the archived Figure_Source_Data/ TSV files (zero hardcoded machine paths).

Usage:
  python scripts/reproduce_core_figures_portable.py --source-dir ./Figure_Source_Data --out-dir reproduced_figures
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CB_PAL = {
    "ifn_high": "#D55E00",
    "ifn_mid": "#E69F00",
    "control": "#0072B2",
    "neoadjuvant": "#D55E00",
    "untreated": "#0072B2",
    "purple": "#6A3D9A",
    "slate": "#64748B",
}

def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def reproduce_fig3b(src_dir: Path, out_dir: Path) -> dict:
    df = pd.read_csv(src_dir / "Figure3b_donor_bootstrap_concordance_source.tsv", sep="\t")
    specs = [
        ("GSE202051_All43_Q2_IFN14_min5", "1_", "GSE202051 All (n=34): all shared", CB_PAL["slate"]),
        ("GSE202051_All43_Q2_IFN14_min5", "2_", "GSE202051 All (n=34): panel", CB_PAL["neoadjuvant"]),
        ("GSE202051_All43_Q2_IFN14_min5", "3_", "GSE202051 All (n=34): disc. DEGs", CB_PAL["neoadjuvant"]),
        ("GSE202051_Untreated18_Q2_IFN14_min5", "2_", "Untreated (n=12): panel", CB_PAL["untreated"]),
        ("GSE202051_Untreated18_Q2_IFN14_min5", "3_", "Untreated (n=12): disc. DEGs", CB_PAL["untreated"]),
        ("GSE202051_Neoadjuvant25_Q2_IFN14_min5", "3_", "Neoadjuvant (n=22): disc. DEGs", "#B45309"),
        ("GSE154778_All16_Q2_IFN14_min5", "1_", "GSE154778 (n=15): all shared", CB_PAL["slate"]),
        ("GSE154778_All16_Q2_IFN14_min5", "2_", "GSE154778 (n=15): panel", CB_PAL["purple"]),
        ("GSE154778_All16_Q2_IFN14_min5", "3_", "GSE154778 (n=15): disc. DEGs", CB_PAL["purple"]),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.1), dpi=300)
    y_pos = np.arange(len(specs))[::-1]
    ax.axvline(0.0, color="#666666", ls="--", lw=0.8)
    ylabels, verified_rows = [], []
    for i, (vm, gp, base_lbl, col) in enumerate(specs):
        yv = y_pos[i]
        rc = df[(df["validation_cohort_model"] == vm) & (df["gene_set_definition"].str.startswith(gp))].iloc[0]
        ng = int(rc["n_genes"])
        rho_e, rho_v = float(rc["edgeR_spearman_rho"]), float(rc["voom_precision_spearman_rho"])
        lo_b, hi_b = float(rc["donor_bootstrap_95CI_low"]), float(rc["donor_bootstrap_95CI_high"])
        lo_l, hi_l = float(rc["LODO_min_spearman_rho"]), float(rc["LODO_max_spearman_rho"])
        ylabels.append(f"{base_lbl} (g={ng:,})")
        ax.plot([lo_b, hi_b], [yv, yv], color=col, lw=1.6, zorder=2)
        ax.plot([lo_l, hi_l], [yv, yv], color="#111111", lw=3.0, alpha=0.35, zorder=3)
        ax.scatter([rho_v], [yv], s=26, color=col, edgecolors="#111111", linewidths=0.5, zorder=4,
                   label="voom rho (matched boot 95% CI & LODO)" if i == 0 else None)
        ax.scatter([rho_e], [yv], s=26, marker="D", facecolors="white", edgecolors=col, linewidths=1.1, zorder=5,
                   label="edgeR QL rho" if i == 0 else None)
        ax.text(0.86, yv, f"voom {rho_v:.2f} [{lo_b:.2f},{hi_b:.2f}] | edgeR {rho_e:.2f}",
                va="center", ha="left", fontsize=6.5, color="#111111")
        verified_rows.append({"model": vm, "tier": gp, "n_genes": ng, "voom_rho": round(rho_v, 4),
                              "voom_boot_ci": [round(lo_b, 4), round(hi_b, 4)], "edgeR_rho": round(rho_e, 4)})
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ylabels, fontsize=7.0)
    ax.set_xlim(-0.25, 1.68)
    ax.set_xlabel("Gene-level paired log2FC Spearman rho (thin: voom boot 95% CI; thick: voom LODO)", fontsize=7.2)
    ax.set_title("Reproduced Figure 3b (exact n_genes & matched voom/edgeR)", fontsize=8.3, fontweight="bold", pad=20)
    ax.legend(loc="upper center", bbox_to_anchor=(0.48, 1.11), ncol=2, frameon=False, fontsize=6.6)
    despine(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "Figure3b_reproduced.png", dpi=300)
    fig.savefig(out_dir / "Figure3b_reproduced.pdf")
    plt.close(fig)
    return {"n_rows_plotted": len(verified_rows), "rows": verified_rows}

def reproduce_fig4c(src_dir: Path, out_dir: Path) -> dict:
    df = pd.read_csv(src_dir / "Figure4c_regvelo_1000perm_macrophage_source.tsv", sep="\t")
    tfs = ["IRF4", "REL", "NFKB1", "NR4A3", "CREB1", "HSF1", "RUNX1", "STAT1+IRF1", "STAT1", "IRF1", "STAT2", "IRF2", "GAPDH"]
    sub = df.set_index("perturbation").loc[[t for t in tfs if t in set(df["perturbation"])]].reset_index()
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=300)
    y = np.arange(len(sub))[::-1]
    ax.axvline(0.0, color="#333333", ls="--", lw=0.8)
    for i, r in sub.iterrows():
        yv = y[i]
        dv = float(r["mean_delta_pp"])
        lo, hi = float(r["donor16_bootstrap_95CI_low_pp"]), float(r["donor16_bootstrap_95CI_high_pp"])
        nsd = float(r["null_matched_sd_pp"]) if pd.notna(r["null_matched_sd_pp"]) else 0.0
        k = int(r["n_muted_prior_edges"])
        p = float(r["empirical_perm_p"]) if pd.notna(r["empirical_perm_p"]) else 1.0
        q = float(r["BH_adjusted_q"]) if pd.notna(r["BH_adjusted_q"]) else 1.0
        col = CB_PAL["control"] if p < 0.05 else (CB_PAL["ifn_mid"] if r["perturbation"] in ("STAT1", "IRF1", "STAT2", "IRF2") else CB_PAL["slate"])
        if nsd > 0:
            ax.barh(yv, 3.92 * nsd, left=-1.96 * nsd, height=0.72, color="#E2E8F0", zorder=1)
        ax.barh(yv, dv, height=0.48, color=col, edgecolor="#111111", linewidth=0.45, zorder=2)
        ax.plot([lo, hi], [yv, yv], color="#111111", lw=1.1, zorder=3)
        x_txt = max(hi, dv, 0.0) + 0.14 if dv >= 0 else min(lo, dv, 0.0) - 0.14
        ha = "left" if dv >= 0 else "right"
        lbl = f"{dv:+.2f} pp (k={k}, p={p:.3f}, q={q:.3f})" if k > 0 else "0.00 pp (0 prior edges)"
        ax.text(x_txt, yv, lbl, va="center", ha=ha, fontsize=6.2)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["perturbation"].tolist(), fontsize=7.5, fontweight="bold")
    ax.set_xlim(-6.2, 6.2)
    ax.set_xlabel("Predicted TAM_STAT1_Inflam fate shift (pp; 3 GPU seeds, 16 patients, N_perm=1,000)", fontsize=7.6)
    ax.set_title("Reproduced Figure 4c: 1,000-permutation degree+expression-matched RegVelo screen", fontsize=8.5, fontweight="bold")
    despine(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "Figure4c_reproduced.png", dpi=300)
    fig.savefig(out_dir / "Figure4c_reproduced.pdf")
    plt.close(fig)
    return {
        "n_tfs_plotted": len(sub),
        "NFKB1_delta_pp": float(sub.loc[sub["perturbation"] == "NFKB1", "mean_delta_pp"].iloc[0]),
        "NFKB1_bh_q": float(sub.loc[sub["perturbation"] == "NFKB1", "BH_adjusted_q"].iloc[0]),
        "STAT1_delta_pp": float(sub.loc[sub["perturbation"] == "STAT1", "mean_delta_pp"].iloc[0]),
        "GAPDH_delta_pp": float(sub.loc[sub["perturbation"] == "GAPDH", "mean_delta_pp"].iloc[0]),
    }

def reproduce_fig6a(src_dir: Path, out_dir: Path) -> dict:
    df = pd.read_csv(src_dir / "Figure6a_tcga_cox_with_BH_q_source.tsv", sep="\t")
    feats = [
        ("Pure_Ductal_PDAC_n146", "FOXM1", "FOXM1 (ductal prolif.)"),
        ("Pure_Ductal_PDAC_n146", "HMGA2", "HMGA2 (basal-like)"),
        ("Pure_Ductal_PDAC_n146", "GATA6", "GATA6 (classical)"),
        ("Pure_Ductal_PDAC_n146", "STAT1", "STAT1 (pure ductal, n=144)"),
        ("Pure_Ductal_PDAC_n146", "IFN14_Signature", "IFN-14 sig (pure ductal)"),
        ("Pure_Ductal_PDAC_n146", "IRF1", "IRF1 (pure ductal)"),
        ("Unfiltered_All_PAAD_n177", "STAT1", "STAT1 (unfiltered n=172*)"),
    ]
    fig, ax = plt.subplots(figsize=(6.8, 3.8), dpi=300)
    y = np.arange(len(feats))[::-1]
    ax.axvline(1.0, color="#333333", ls="--", lw=0.8)
    out_summary = {}
    for i, (coh, fn, lbl) in enumerate(feats):
        yv = y[i]
        rf = df[(df["cohort"] == coh) & (df["endpoint"] == "OS") & (df["feature"] == fn)].iloc[0]
        hr2 = float(rf["HR_multi"])
        lo2 = float(rf["CI95_low_multi"] if "CI95_low_multi" in rf else rf["ci95_multi_low"])
        hi2 = float(rf["CI95_high_multi"] if "CI95_high_multi" in rf else rf["ci95_multi_high"])
        p2, q2 = float(rf["p_multi"]), float(rf["BH_q_multi"])
        col = CB_PAL["ifn_high"] if (p2 < 0.05 and hr2 > 1) else (CB_PAL["control"] if (p2 < 0.05 and hr2 < 1) else CB_PAL["slate"])
        if coh == "Unfiltered_All_PAAD_n177":
            col = CB_PAL["ifn_mid"]
        ax.plot([lo2, hi2], [yv, yv], color=col, lw=1.6)
        ax.scatter([hr2], [yv], s=28, color=col, edgecolors="#111111", linewidths=0.5, zorder=4)
        ax.text(3.28, yv, f"HR {hr2:.2f} [{lo2:.2f},{hi2:.2f}] (p={p2:.4f}, q={q2:.3f})", va="center", ha="right", fontsize=6.4, color=col)
        out_summary[f"{coh}:{fn}"] = {"HR": round(hr2, 4), "p": round(p2, 5), "BH_q": round(q2, 4)}
    ax.set_yticks(y)
    ax.set_yticklabels([f[2] for f in feats], fontsize=7.4)
    ax.set_xlim(0.30, 3.35)
    ax.set_xlabel("Multivariable OS Hazard Ratio per 1-SD log2(RSEM+1) (95% CI)", fontsize=7.6)
    ax.set_title("Reproduced Figure 6a: Pure-ductal TCGA-PAAD multivariable Cox OS with BH q", fontsize=8.5, fontweight="bold")
    despine(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "Figure6a_reproduced.png", dpi=300)
    fig.savefig(out_dir / "Figure6a_reproduced.pdf")
    plt.close(fig)
    return out_summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("reproduced_figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "Figure3b": reproduce_fig3b(args.source_dir, args.out_dir),
        "Figure4c": reproduce_fig4c(args.source_dir, args.out_dir),
        "Figure6a": reproduce_fig6a(args.source_dir, args.out_dir),
    }
    (args.out_dir / "reproduction_check_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("Successfully reproduced Figure3b, Figure4c, and Figure6a into:", args.out_dir.resolve())

if __name__ == "__main__":
    main()
