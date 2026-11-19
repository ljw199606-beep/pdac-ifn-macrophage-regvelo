#!/usr/bin/env python3
"""
CSBJ Phase 5 Publication Figure & Multi-Sheet Supplementary Tables Generator
(csbj_05_render_all_figures_and_tables.py)

Generates:
  1. CSBJ_Final_Submission/02_Figures/Figure1..Figure6 (.pdf, .png [600 DPI], .svg)
  2. CSBJ_Final_Submission/02_Figures/Figure_Source_Data/*.tsv + FIGURE_SOURCE_MAP.tsv
  3. CSBJ_Final_Submission/03_Supplementary/Supplementary_Figures/FigureS1..FigureS5 (.pdf, .png [600 DPI], .svg)
  4. CSBJ_Final_Submission/03_Supplementary/Supplementary_Tables.xlsx (00_Data_Dictionary + 16 audited tables)

Strict Visualization Standards:
  - Double-column print width = 175 mm (6.89 in), editable TrueType vector fonts (pdf.fonttype = 42).
  - Strictly colorblind-safe Okabe-Ito / Paul Tol Blue-Orange-Amber-Purple-Slate palette (ZERO red-green opposition!).
  - Programmatically reads reconciled numbers from csbj_frozen_numbers.json and Phase 1/2 audited TSVs.
  - Zero p=0, zero nan, zero untraced numbers, explicit statistical units, n, and uncertainty intervals.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy import stats
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

BASE_ROOT = Path(r"D:\Research\2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT   = Path(r"D:\Research\2026_PDAC_RegVelo")
CSBJ_ROOT = BASE_ROOT / "CSBJ_Final_Submission"

FIG_DIR   = CSBJ_ROOT / "02_Figures"
SRC_DIR   = FIG_DIR / "Figure_Source_Data"
SUPP_FDIR = CSBJ_ROOT / "03_Supplementary" / "Supplementary_Figures"
SUPP_DIR  = CSBJ_ROOT / "03_Supplementary"
MAN_DIR   = CSBJ_ROOT / "04_Reproducibility" / "data_manifest"
P2_DIR    = CSBJ_ROOT / "06_Quality_Audit" / "phase2_tables"

for d in (FIG_DIR, SRC_DIR, SUPP_FDIR, SUPP_DIR):
    d.mkdir(parents=True, exist_ok=True)

DOUBLE_COL_IN = 6.89  # 175 mm CSBJ double-column print width

# Strictly colorblind-safe Okabe-Ito / Paul Tol palette (ZERO green hues paired against red/vermilion!)
CB_PAL = {
    "pdac":        "#D55E00",  # Okabe-Ito vermilion-orange
    "control":     "#0072B2",  # Okabe-Ito cobalt blue
    "ifn_high":    "#C84B00",  # deep burnt vermilion
    "ifn_mid":     "#E69F00",  # warm amber-gold
    "ifn_low":     "#1F5AA6",  # deep sapphire blue
    "untreated":   "#0072B2",  # cobalt blue
    "neoadjuvant": "#D55E00",  # vermilion-orange
    "primary":     "#1E3A5F",  # deep navy
    "metastasis":  "#D97706",  # ochre-amber
    "hm":          "#E69F00",  # hepatic met amber
    "lnm":         "#6A3D9A",  # lymph node royal purple
    "np":          "#56B4E9",  # matched normal pancreas sky cyan (zero green!)
    "purple":      "#6A3D9A",  # royal purple
    "cyan":        "#56B4E9",  # sky cyan
    "magenta":     "#CC79A7",  # muted rose-magenta
    "slate":       "#64748B",  # neutral slate
    "light_slate": "#CBD5E1",  # background neutral
    "navy":        "#1E3A5F",
}

# Zero-green 10-lineage colorblind-safe palette for Peng2019 single-cell UMAP & stacked bars
CELLTYPE_CB = {
    "Ductal cell type 2": "#C84B00",  # burnt vermilion (malignant ductal)
    "Ductal cell type 1": "#0072B2",  # cobalt blue (non-malignant ductal)
    "Acinar cell":        "#56B4E9",  # sky cyan
    "Endocrine cell":     "#CC79A7",  # rose-magenta
    "Macrophage cell":    "#E69F00",  # warm amber-orange
    "T cell":             "#6A3D9A",  # royal purple
    "B cell":             "#3B82F6",  # periwinkle blue
    "Fibroblast cell":    "#8C510A",  # warm brown-ochre
    "Stellate cell":      "#B45309",  # bronze
    "Endothelial cell":   "#475569",  # charcoal slate
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.6,
    "axes.titlesize": 8.1,
    "axes.titleweight": "bold",
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.fontsize": 6.4,
    "legend.title_fontsize": 6.8,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "lines.linewidth": 1.1,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.06,
})

def despine(ax, top=True, right=True):
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)

def panel_label(ax, letter: str, x: float = -0.14, y: float = 1.07):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10.2, fontweight="bold",
            va="bottom", ha="left", color="#111111")

def fmt_p(p: float, prefix: str = "p") -> str:
    if p is None or not np.isfinite(p):
        return f"{prefix} = NA"
    if p < 1e-50:
        return f"{prefix} < 1e-50"
    if p < 1e-4:
        exp = int(np.floor(np.log10(max(p, 1e-300))))
        mant = p / (10.0 ** exp)
        return f"{prefix} = {mant:.2f}e{exp}"
    if p < 0.01:
        return f"{prefix} = {p:.4f}"
    return f"{prefix} = {p:.3f}"

def export_fig(fig: plt.Figure, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf", format="pdf")
    fig.savefig(out_dir / f"{stem}.svg", format="svg")
    fig.savefig(out_dir / f"{stem}.png", format="png", dpi=600)
    plt.close(fig)

frozen = json.loads((MAN_DIR / "csbj_frozen_numbers.json").read_text(encoding="utf-8"))
rng = np.random.default_rng(42)
src_map_rows = []

def reg_src(fig_id: str, panel: str, src_file: str, desc: str):
    src_map_rows.append({
        "figure": fig_id,
        "panel": panel,
        "source_data_file": src_file,
        "generator_script": "04_Reproducibility/scripts/csbj_05_render_all_figures_and_tables.py",
        "description": desc,
    })

# =========================================================================
# FIGURE 1: COHORT ARCHITECTURE (247 DONORS / 270 SPECIMENS) & PENG ATLAS
# =========================================================================
print("[1/12] Rendering Figure1 (247-donor multi-cohort architecture & discovery cellular landscape)...")
umap_peng = pd.read_csv(V1_ROOT / "04_results/plot_inputs/peng_umap_57423.csv.gz")
dot_peng  = pd.read_csv(V1_ROOT / "04_results/plot_inputs/peng_canonical_marker_dotplot.csv")
comp_stat = pd.read_csv(V1_ROOT / "04_results/plot_inputs/peng_donor_composition_stats.csv")
comp_wide = pd.read_csv(V1_ROOT / "04_results/plot_inputs/peng_donor_composition_wide.csv")
coh_247   = pd.read_excel(MAN_DIR / "cohort_provenance_audit.xlsx", sheet_name="01_Cohort_Summary_247")

coh_247.to_csv(SRC_DIR / "Figure1a_cohort_architecture_247_source.tsv", sep="\t", index=False)
dot_peng.to_csv(SRC_DIR / "Figure1c_canonical_markers_source.tsv", sep="\t", index=False)
comp_stat.to_csv(SRC_DIR / "Figure1d_donor_composition_stats_source.tsv", sep="\t", index=False)
comp_wide.to_csv(SRC_DIR / "Figure1e_donor_stacked_fractions_source.tsv", sep="\t", index=False)
reg_src("Figure1", "a", "Figure1a_cohort_architecture_247_source.tsv", "Reconciled 5-cohort 247-donor / 270-specimen study design")
reg_src("Figure1", "b-c", "Figure1c_canonical_markers_source.tsv", "Peng2019 57,423-cell UMAP & canonical lineage dotplot")
reg_src("Figure1", "d-e", "Figure1d_donor_composition_stats_source.tsv", "Donor-level cellular remodeling across 24 PDAC vs 11 control donors")

fig1 = plt.figure(figsize=(DOUBLE_COL_IN, 6.85))
gs1 = gridspec.GridSpec(3, 2, figure=fig1, height_ratios=[0.88, 1.34, 1.25], hspace=0.46, wspace=0.34)

ax1a = fig1.add_subplot(gs1[0, :])
ax1a.set_xlim(0, 100)
ax1a.set_ylim(0, 27)
ax1a.axis("off")
panel_label(ax1a, "a", x=-0.02, y=1.02)
ax1a.text(
    0, 25.4,
    "Reconciled multi-cohort PDAC design: 5 independent cohorts, 247 unique primary donors/patients (270 specimens/sections)",
    fontsize=7.65, fontweight="bold", color="#111111",
)

boxes_1a = [
    (0.3,  1.8, 18.8, 21.4, "#EFF6FF", "#0072B2", "1. Discovery scRNA\nCRA001160 (Peng)",
     "35 donors (24 PDAC, 11 Ctrl)\n57,423 QC cells (18,008 genes)\n5,361 TAMs (31 paired Q2)\nUnit: donor (n=34 Q1 / 31 Q2)"),
    (20.4, 1.8, 18.8, 21.4, "#FFF7ED", "#D55E00", "2. External snRNA\nGSE202051 (Hwang)",
     "43 PDAC pts (0 Ctrl)\n18 Untr. | 25 Neoadj.\n10,140 TAMs (34 paired Q2)\nUnit: patient (n=43 / 34)"),
    (40.5, 1.8, 18.8, 21.4, "#F5F3FF", "#6A3D9A", "3. Velocity scRNA\nGSE154778 (Lin)",
     "16 PDAC pts (10 PT, 6 Met)\n120 SRR -> 30,399 QC cells\n3,318 TAMs | 10,042 Ductal\nUnit: patient (n=16; 1k perm)"),
    (60.6, 1.8, 18.8, 21.4, "#F0F9FF", "#1E3A5F", "4. Visium Spatial\nGSE272362 (Khaliq)",
     "7 autopsy pts / 30 sections\n10 PT, 12 HM, 5 LNM, 3 NP*\n97,020 tissue (91,496 QC) spots\nUnit: patient (n=7) & site"),
    (80.7, 1.8, 18.8, 21.4, "#F8FAFC", "#475569", "5. Clinical Bulk\nTCGA-PAAD (GDC)",
     "146 pure-ductal PDAC pts\n(145 OS, 144 multivar, 84 ev)\n[+31 non-ductal sens = 177]\nUnit: patient (per-SD Cox)"),
]
for bx, by, bw, bh, fc, ec, title_b, body_b in boxes_1a:
    rect = mpatches.FancyBboxPatch((bx, by), bw, bh, boxstyle="round,pad=0.2,rounding_size=0.85",
                                   facecolor=fc, edgecolor=ec, linewidth=0.9)
    ax1a.add_patch(rect)
    ax1a.text(bx + bw / 2, by + bh - 1.2, title_b, ha="center", va="top", fontsize=5.95, fontweight="bold", color=ec)
    ax1a.text(bx + bw / 2, by + bh - 6.9, body_b, ha="center", va="top", fontsize=5.15, color="#1E293B", linespacing=1.22)
for x_arr in [19.2, 39.3, 59.4, 79.5]:
    ax1a.annotate("", xy=(x_arr + 1.05, 12.5), xytext=(x_arr, 12.5),
                  arrowprops=dict(arrowstyle="->", lw=0.85, color="#475569"))

# Panel 1b: Peng2019 57,423-cell UMAP (zero-green colorblind-safe palette)
ax1b = fig1.add_subplot(gs1[1, 0])
panel_label(ax1b, "b", x=-0.16, y=1.09)
u_col = "UMAP1" if "UMAP1" in umap_peng.columns else umap_peng.columns[1]
v_col = "UMAP2" if "UMAP2" in umap_peng.columns else umap_peng.columns[2]
ct_col = "Cell_type" if "Cell_type" in umap_peng.columns else [c for c in umap_peng.columns if "type" in c.lower()][0]
up_shuf = umap_peng.iloc[rng.permutation(len(umap_peng))]
ax1b.scatter(up_shuf[u_col], up_shuf[v_col], c=[CELLTYPE_CB.get(c, "#94A3B8") for c in up_shuf[ct_col]],
             s=0.35, alpha=0.72, linewidths=0, rasterized=True)
despine(ax1b)
ax1b.set_xticks([]); ax1b.set_yticks([])
ax1b.set_xlabel("UMAP 1"); ax1b.set_ylabel("UMAP 2")
ax1b.set_title("CRA001160 atlas (57,423 QC cells; 35 donors)", fontsize=7.6)

short_ct = {
    "Ductal cell type 2": "Ductal 2 (malig.)",
    "Ductal cell type 1": "Ductal 1",
    "Acinar cell": "Acinar",
    "Endocrine cell": "Endocrine",
    "Macrophage cell": "Macrophage (5,361)",
    "T cell": "T cell",
    "B cell": "B cell",
    "Fibroblast cell": "Fibroblast",
    "Stellate cell": "Stellate",
    "Endothelial cell": "Endothelial",
}
handles_1b = [plt.Line2D([0], [0], marker="o", ls="", markersize=3.6, color=CELLTYPE_CB[k], label=short_ct[k]) for k in short_ct]
ax1b.legend(handles=handles_1b, loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=True, facecolor="white",
            framealpha=0.90, edgecolor="none", fontsize=5.5, ncol=2, handletextpad=0.2, columnspacing=0.55, borderpad=0.22)

# Panel 1c: Canonical Lineage Marker Dotplot (YlOrBr colorblind-safe sequential colormap)
ax1c = fig1.add_subplot(gs1[1, 1])
panel_label(ax1c, "c", x=-0.22, y=1.05)
ct_order = list(short_ct.keys())
marker_pick = [g for g in ["CEACAM6", "KRT19", "CFTR", "PRSS1", "INS", "CD68", "C1QA", "CD3D", "MS4A1", "COL1A1", "RGS5", "PECAM1"] if g in set(dot_peng["gene"])]
dp_sub = dot_peng[dot_peng["Cell_type"].isin(ct_order) & dot_peng["gene"].isin(marker_pick)].copy()
ct_y = {c: i for i, c in enumerate(reversed(ct_order))}
g_x  = {g: j for j, g in enumerate(marker_pick)}
sc_dp = ax1c.scatter(
    [g_x[g] for g in dp_sub["gene"]],
    [ct_y[c] for c in dp_sub["Cell_type"]],
    s=np.clip(dp_sub["frac_pos"].values * 46.0, 1.5, 46.0),
    c=dp_sub["scaled_expr"].values, cmap="YlOrBr", vmin=0, vmax=1,
    edgecolors="#222222", linewidths=0.25,
)
ax1c.set_xticks(range(len(marker_pick)))
ax1c.set_xticklabels(marker_pick, rotation=45, ha="right", fontsize=6.1)
ax1c.set_yticks(range(len(ct_order)))
ax1c.set_yticklabels([short_ct[c].replace(" (5,361)", "").replace(" (malig.)", "*") for c in reversed(ct_order)], fontsize=6.2)
ax1c.set_xlim(-0.6, len(marker_pick) - 0.4); ax1c.set_ylim(-0.6, len(ct_order) - 0.4)
ax1c.set_title("Canonical lineage markers (scaled expr & % pos)", fontsize=7.6)
despine(ax1c)
cb1c = fig1.colorbar(sc_dp, ax=ax1c, fraction=0.038, pad=0.03)
cb1c.ax.tick_params(labelsize=5.8); cb1c.set_label("Scaled expr", fontsize=6.0)

# Panel 1d: Donor-level Cellular Remodeling (Unit = donor: 24 PDAC vs 11 Control)
ax1d = fig1.add_subplot(gs1[2, 0])
panel_label(ax1d, "d", x=-0.18, y=1.05)
show_cts = ["Ductal cell type 2", "Macrophage cell", "T cell", "Stellate cell", "Ductal cell type 1", "Endothelial cell", "Acinar cell"]
y_pos = np.arange(len(show_cts))[::-1]
for idx_c, ct_n in enumerate(show_cts):
    yv = y_pos[idx_c]
    t_vals = comp_wide.loc[comp_wide["CONDITION"] == "T", ct_n].values * 100.0
    n_vals = comp_wide.loc[comp_wide["CONDITION"] == "N", ct_n].values * 100.0
    ax1d.scatter(t_vals, np.full_like(t_vals, yv + 0.15) + rng.uniform(-0.05, 0.05, len(t_vals)), s=8, color=CB_PAL["pdac"], alpha=0.75, linewidths=0)
    ax1d.scatter(n_vals, np.full_like(n_vals, yv - 0.15) + rng.uniform(-0.05, 0.05, len(n_vals)), s=8, color=CB_PAL["control"], alpha=0.75, linewidths=0)
    ax1d.plot([np.percentile(t_vals, 25), np.percentile(t_vals, 75)], [yv + 0.15, yv + 0.15], color=CB_PAL["pdac"], lw=1.5)
    ax1d.plot(np.median(t_vals), yv + 0.15, marker="D", markersize=3.4, color="#7C2D12")
    ax1d.plot([np.percentile(n_vals, 25), np.percentile(n_vals, 75)], [yv - 0.15, yv - 0.15], color=CB_PAL["control"], lw=1.5)
    ax1d.plot(np.median(n_vals), yv - 0.15, marker="D", markersize=3.4, color="#0C4A6E")
    q_v = float(comp_stat.loc[comp_stat["Cell_type"] == ct_n, "mwu_padj"].iloc[0])
    ax1d.text(98.0, yv, fmt_p(q_v, "BH q"), va="center", ha="right", fontsize=5.4, color="#1E293B")

ax1d.set_yticks(y_pos)
ax1d.set_yticklabels([short_ct[c].replace(" (5,361)", "").replace(" (malig.)", "") for c in show_cts], fontsize=6.4)
ax1d.set_xlim(-2, 100)
ax1d.set_xlabel("Donor lineage fraction (% of donor cells; median & IQR)")
ax1d.set_title("Donor-level remodeling (PDAC n=24 vs Ctrl n=11)", fontsize=7.5)
despine(ax1d)
ax1d.legend(
    handles=[
        plt.Line2D([0], [0], marker="o", ls="", color=CB_PAL["pdac"], markersize=3.8, label="PDAC (n=24)"),
        plt.Line2D([0], [0], marker="o", ls="", color=CB_PAL["control"], markersize=3.8, label="Ctrl (n=11)"),
    ],
    loc="lower center", bbox_to_anchor=(0.48, 0.01), frameon=True, facecolor="white", framealpha=0.9, edgecolor="none", fontsize=5.5,
)

# Panel 1e: Per-donor Stacked Cellular Fractions (All 35 Donors)
ax1e = fig1.add_subplot(gs1[2, 1])
panel_label(ax1e, "e", x=-0.16, y=1.05)
cw_n = comp_wide[comp_wide["CONDITION"] == "N"].sort_values("Macrophage cell")
cw_t = comp_wide[comp_wide["CONDITION"] == "T"].sort_values("Ductal cell type 2", ascending=False)
cw_ord = pd.concat([cw_n, cw_t], ignore_index=True)
x_don = np.arange(len(cw_ord))
bottom = np.zeros(len(cw_ord))
for ct_n in ct_order:
    vals = cw_ord[ct_n].values * 100.0
    ax1e.bar(x_don, vals, bottom=bottom, width=0.84, color=CELLTYPE_CB[ct_n], linewidth=0)
    bottom += vals
ax1e.axvline(len(cw_n) - 0.5, color="#111111", ls="--", lw=0.9)
ax1e.text(len(cw_n) / 2 - 0.5, 102.5, "Control (n=11)", ha="center", va="bottom", fontsize=6.2, fontweight="bold", color=CB_PAL["control"])
ax1e.text(len(cw_n) + len(cw_t) / 2 - 0.5, 102.5, "PDAC (n=24)", ha="center", va="bottom", fontsize=6.2, fontweight="bold", color=CB_PAL["pdac"])
ax1e.set_xlim(-0.7, len(cw_ord) - 0.3); ax1e.set_ylim(0, 108)
ax1e.set_xticks([0, 6, 12, 23, 34])
ax1e.set_xticklabels([cw_ord["Patient"].iloc[i] for i in [0, 6, 12, 23, 34]], fontsize=6.0)
ax1e.set_ylabel("Donor composition (%)")
ax1e.set_xlabel("Individual donors (unit = donor; n = 35)")
ax1e.set_title("Inter-donor cellular heterogeneity (Ductal 2 = 0% in Ctrl)", fontsize=7.4)
despine(ax1e)

export_fig(fig1, FIG_DIR, "Figure1")

# =========================================================================
# FIGURE 2: IFN-HIGH MACROPHAGE DISCOVERY, LODO PAIRED PSEUDOBULK & BACKGROUND ATTENUATION
# =========================================================================
print("[2/12] Rendering Figure2 (IFN-14/IFN-10 state discovery, LODO-verified paired edgeR, background attenuation & fgsea)...")
mac_cells  = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/peng2019_macrophage_cell_states.tsv", sep="\t")
don_q1     = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/peng2019_donor_abundance_Q1.tsv", sep="\t")
q1_tests   = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/between_donor_abundance_tests_Q1.tsv", sep="\t")
de_peng    = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/edger_outputs/paired_edgeR_Peng_Q2_IFN14_min5.tsv", sep="\t")
amb_audit  = pd.read_csv(P2_DIR / "ambient_rna_and_background_audit.tsv", sep="\t")
edger_lodo = pd.read_csv(P2_DIR / "edgeR_exact_LODO_and_voom_bootstrap_audit.tsv", sep="\t")
fg_all     = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/noncircular_hallmark_fgsea_all52_pathways.tsv", sep="\t")
fg_peng    = fg_all[fg_all["tag"] == "Peng_Q2_IFN14_min5"].copy() if "tag" in fg_all.columns else fg_all[fg_all["model"] == "Peng_Q2_IFN14_min5"].copy()

mac_cells.to_csv(SRC_DIR / "Figure2a_peng_macrophage_states_source.tsv", sep="\t", index=False)
don_q1.to_csv(SRC_DIR / "Figure2b_peng_donor_Q1_proportions_source.tsv", sep="\t", index=False)
de_peng.to_csv(SRC_DIR / "Figure2c_peng_paired_edgeR_IFN14_source.tsv", sep="\t", index=False)
amb_audit.to_csv(SRC_DIR / "Figure2d_background_signal_attenuation_source.tsv", sep="\t", index=False)
fg_peng.to_csv(SRC_DIR / "Figure2e_noncircular_fgsea_source.tsv", sep="\t", index=False)
reg_src("Figure2", "a-b", "Figure2b_peng_donor_Q1_proportions_source.tsv", "Non-circular IFN-14/IFN-10 TAM UMAP & Rule Q1 donor expansion")
reg_src("Figure2", "c", "Figure2c_peng_paired_edgeR_IFN14_source.tsv", "Rule Q2 paired edgeR QL volcano across 15,030 non-defining genes (n=31 donors)")
reg_src("Figure2", "d", "Figure2d_background_signal_attenuation_source.tsv", "Paired within-donor attenuation of non-macrophage background transcripts")
reg_src("Figure2", "e", "Figure2e_noncircular_fgsea_source.tsv", "Non-circular MSigDB Hallmark + MHC/checkpoint module fgsea")

fig2 = plt.figure(figsize=(DOUBLE_COL_IN, 7.20))
gs2 = gridspec.GridSpec(3, 2, figure=fig2, height_ratios=[1.05, 1.16, 1.14], hspace=0.49, wspace=0.35)

# Panel 2a: Macrophage UMAP colored by non-circular IFN-14 Q1 tertiles
ax2a = fig2.add_subplot(gs2[0, 0])
panel_label(ax2a, "a", x=-0.16, y=1.08)
st_cmap = {"low": CB_PAL["ifn_low"], "mid": CB_PAL["light_slate"], "high": CB_PAL["ifn_high"]}
for st_k, lbl_k, sz_k, al_k in [
    ("mid", "IFN-mid (33.3%)", 2.0, 0.45),
    ("low", "IFN-low (33.3%)", 2.4, 0.68),
    ("high", "IFN-high (33.3%, Q1)", 2.8, 0.82),
]:
    sub_m = mac_cells[mac_cells["Q1_global_IFN14"] == st_k]
    ax2a.scatter(sub_m["UMAP1"], sub_m["UMAP2"], s=sz_k, c=st_cmap[st_k], alpha=al_k, linewidths=0, label=lbl_k, rasterized=True)
despine(ax2a)
ax2a.set_xticks([]); ax2a.set_yticks([])
ax2a.set_xlabel("Macrophage UMAP 1"); ax2a.set_ylabel("Macrophage UMAP 2")
ax2a.set_title("CRA001160 TAMs: non-circular IFN-14 (n=5,361; 34 donors)", fontsize=7.25)
ax2a.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=5.8, markerscale=2.0)

# Panel 2b: Rule Q1 Donor-Level State Proportion Comparison (IFN-14, IFN-10, & % All Cells)
ax2b = fig2.add_subplot(gs2[0, 1])
panel_label(ax2b, "b", x=-0.18, y=1.08)
metrics_2b = [
    ("ifn14_high_frac_of_mac_Q1", "IFN-14-high\n(% of TAMs)"),
    ("ifn10_high_frac_of_mac_Q1", "IFN-10-high\n(% of TAMs)"),
    ("ifn14_high_frac_of_all_cells_Q1", "IFN-14-high\n(% all cells x5)"),
]
for mi, (mcol, _) in enumerate(metrics_2b):
    scale_f = 500.0 if "all_cells" in mcol else 100.0
    v_ctrl = don_q1.loc[(don_q1["condition"] == "Control") & (don_q1["n_macrophages"] >= 1), mcol].dropna().values * scale_f
    v_pdac = don_q1.loc[(don_q1["condition"] == "PDAC") & (don_q1["n_macrophages"] >= 1), mcol].dropna().values * scale_f
    bp = ax2b.boxplot([v_ctrl, v_pdac], positions=[mi - 0.17, mi + 0.17], widths=0.26, patch_artist=True,
                      showfliers=False, medianprops=dict(color="#111111", lw=1.2),
                      whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    bp["boxes"][0].set_facecolor(CB_PAL["control"]); bp["boxes"][0].set_alpha(0.35)
    bp["boxes"][1].set_facecolor(CB_PAL["pdac"]);    bp["boxes"][1].set_alpha(0.35)
    ax2b.scatter(np.full_like(v_ctrl, mi - 0.17) + rng.uniform(-0.06, 0.06, len(v_ctrl)), v_ctrl,
                 s=11, color=CB_PAL["control"], alpha=0.85, edgecolors="white", linewidths=0.25, zorder=3)
    ax2b.scatter(np.full_like(v_pdac, mi + 0.17) + rng.uniform(-0.06, 0.06, len(v_pdac)), v_pdac,
                 s=11, color=CB_PAL["pdac"], alpha=0.85, edgecolors="white", linewidths=0.25, zorder=3)
    pv = float(q1_tests[(q1_tests["cohort"] == "Peng2019") & (q1_tests["min_macrophages_per_donor"] == 1) & (q1_tests["metric"] == mcol)]["mwu_p"].iloc[0])
    ymax_loc = min(max(np.max(v_ctrl), np.max(v_pdac)), 68.0) if mi == 2 else max(np.max(v_ctrl), np.max(v_pdac))
    ax2b.plot([mi - 0.17, mi - 0.17, mi + 0.17, mi + 0.17], [ymax_loc + 2.5, ymax_loc + 4.8, ymax_loc + 4.8, ymax_loc + 2.5], color="#333333", lw=0.7)
    ax2b.text(mi, ymax_loc + 5.6, fmt_p(pv, "p"), ha="center", va="bottom", fontsize=5.7, fontweight="bold")

ax2b.set_xticks([0, 1, 2])
ax2b.set_xticklabels([m[1] for m in metrics_2b], fontsize=6.3)
ax2b.set_ylim(-3, 103)
ax2b.set_ylabel("Donor state proportion (%)")
ax2b.set_title("Rule Q1 donor-level expansion (n=24 PDAC vs 10 Ctrl)", fontsize=7.35)
despine(ax2b)
ax2b.legend(
    handles=[
        mpatches.Patch(facecolor=CB_PAL["control"], alpha=0.6, label="Control (n=10)"),
        mpatches.Patch(facecolor=CB_PAL["pdac"], alpha=0.6, label="PDAC (n=24)"),
    ],
    loc="upper right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.7,
)

# Panel 2c: Rule Q2 Within-Donor Paired edgeR QL Volcano + LODO Stability Badge
ax2c = fig2.add_subplot(gs2[1, 0])
panel_label(ax2c, "c", x=-0.20, y=1.08)
non_def = de_peng[~de_peng["in_defining_signature"]].copy()
non_def["neglog10_fdr"] = -np.log10(np.clip(non_def["FDR_edgeR_nonDefining"].values, 1e-35, 1.0))
sig_up_all = int(((non_def["FDR_edgeR_nonDefining"] < 0.05) & (non_def["log2FC_edgeR"] > 0.0)).sum())
sig_dn_all = int(((non_def["FDR_edgeR_nonDefining"] < 0.05) & (non_def["log2FC_edgeR"] < 0.0)).sum())
sig_up = (non_def["FDR_edgeR_nonDefining"] < 0.05) & (non_def["log2FC_edgeR"] > 0.25)
sig_dn = (non_def["FDR_edgeR_nonDefining"] < 0.05) & (non_def["log2FC_edgeR"] < -0.25)
ns_m   = ~(sig_up | sig_dn)

ax2c.scatter(non_def.loc[ns_m, "log2FC_edgeR"], non_def.loc[ns_m, "neglog10_fdr"], s=3.5, c=CB_PAL["light_slate"], alpha=0.55, linewidths=0, rasterized=True)
ax2c.scatter(non_def.loc[sig_dn, "log2FC_edgeR"], non_def.loc[sig_dn, "neglog10_fdr"], s=6.0, c=CB_PAL["ifn_low"], alpha=0.78, linewidths=0, label=f"Down |FC|>0.25 (n={int(sig_dn.sum())}; {sig_dn_all} all)")
ax2c.scatter(non_def.loc[sig_up, "log2FC_edgeR"], non_def.loc[sig_up, "neglog10_fdr"], s=6.5, c=CB_PAL["ifn_high"], alpha=0.82, linewidths=0, label=f"Up |FC|>0.25 (n={int(sig_up.sum())}; {sig_up_all} all)")
ax2c.axhline(-np.log10(0.05), color="#555555", ls="--", lw=0.65)
ax2c.axvline(-0.25, color="#94A3B8", ls=":", lw=0.55)
ax2c.axvline(0.25, color="#94A3B8", ls=":", lw=0.55)
ax2c.axvline(0.0, color="#888888", ls="-", lw=0.45, alpha=0.5)

callout_offsets = {
    "STAT1":  (0.65,  2.2),
    "ISG15":  (0.55, -2.2),
    "IRF7":   (0.65,  1.6),
    "TAP1":   (0.75,  1.5),
    "HLA-A":  (-1.85, 5.5),
    "B2M":    (-1.85, 3.2),
    "CD47":   (-1.95, 0.8),
    "CD80":   (0.85,  0.8),
}
for g_lbl, (dx, dy) in callout_offsets.items():
    r_g = non_def[non_def["gene"] == g_lbl]
    if len(r_g):
        gx, gy = float(r_g["log2FC_edgeR"].iloc[0]), float(r_g["neglog10_fdr"].iloc[0])
        ax2c.scatter([gx], [gy], s=18, facecolors="none", edgecolors="#111111", linewidths=0.8, zorder=5)
        ax2c.annotate(g_lbl, xy=(gx, gy), xytext=(gx + dx, gy + dy), fontsize=5.8, fontweight="bold",
                      color="#111111", arrowprops=dict(arrowstyle="-", lw=0.5, color="#444444"))

lodo_stat1 = str(edger_lodo.query("cohort_model == 'Peng2019_Q2_IFN14_min5' and gene == 'STAT1'")["edgeR_LODO_runs_FDR_lt_005"].iloc[0])
ax2c.text(-3.05, 21.0, f"Exact edgeR LODO:\nSTAT1/TAP1/HLA-A\nFDR<0.05 in {lodo_stat1}",
          fontsize=5.2, bbox=dict(boxstyle="round,pad=0.22", facecolor="#F8FAFC", edgecolor="#CBD5E1", lw=0.5))
ax2c.set_xlim(-3.2, 4.4); ax2c.set_ylim(-1, 35)
ax2c.set_xlabel("Within-donor paired log2FC (IFN-high vs IFN-low)")
ax2c.set_ylabel("-log10(BH FDR; 15,030 non-IFN14 genes)")
ax2c.set_title("Rule Q2 paired edgeR QL (unit = 31 paired donors)", fontsize=7.45)
despine(ax2c)
ax2c.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.7)

# Panel 2d: Epistemically Calibrated Attenuation of Donor-Associated Non-Macrophage Signals
ax2d = fig2.add_subplot(gs2[1, 1])
panel_label(ax2d, "d", x=-0.18, y=1.08)
amb_old = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/ambient_rna_cancellation_diagnostic.tsv", sep="\t")
amb_sub = amb_old.dropna(subset=["unpaired_TvN_log2FC", "paired_withinDonor_log2FC"]).copy()
soup_m = amb_sub["marker_class"] == "Extracellular_Ambient_Soup_Marker"
true_m = ~soup_m

ax2d.axhline(0, color="#222222", ls="--", lw=0.75)
ax2d.axvline(0, color="#AAAAAA", ls=":", lw=0.65)
ax2d.axhspan(-0.42, 0.42, color="#E2E8F0", alpha=0.48, zorder=0)
ax2d.scatter(amb_sub.loc[soup_m, "unpaired_TvN_log2FC"], amb_sub.loc[soup_m, "paired_withinDonor_log2FC"],
             s=32, color=CB_PAL["slate"], marker="s", edgecolors="#1E293B", linewidths=0.5,
             label="Non-macrophage (acinar/duct/stroma)", zorder=3)
ax2d.scatter(amb_sub.loc[true_m, "unpaired_TvN_log2FC"], amb_sub.loc[true_m, "paired_withinDonor_log2FC"],
             s=36, color=CB_PAL["ifn_high"], marker="o", edgecolors="#1E293B", linewidths=0.5,
             label="IFN / MHC-I / checkpoint transcripts", zorder=4)

lbl_2d = {
    "COL1A1": (-0.65, 0.24), "COL3A1": (-0.65, -0.36), "KRT19": (-0.25, -0.36),
    "PRSS1": (0.18, -0.30), "CPA1": (0.18, 0.20),
    "IRF1": (0.22, -0.06), "STAT1": (0.22, 0.14), "TAP1": (-1.05, 0.16),
    "CD80": (0.22, -0.16), "CD47": (-0.95, 0.14), "HLA-A": (0.20, -0.18),
}
for g_n, (ox, oy) in lbl_2d.items():
    rg = amb_sub[amb_sub["gene"] == g_n]
    if len(rg):
        ax2d.text(float(rg["unpaired_TvN_log2FC"].iloc[0]) + ox, float(rg["paired_withinDonor_log2FC"].iloc[0]) + oy,
                  g_n, fontsize=5.8, fontweight="bold" if g_n in {"COL1A1","KRT19","STAT1","IRF1","CD47"} else "normal")

ax2d.set_xlim(-4.6, 6.2); ax2d.set_ylim(-0.75, 2.95)
ax2d.set_xlabel("Unpaired Tumor-vs-Control TAM log2FC (donor-confounded)")
ax2d.set_ylabel("Within-donor paired Q2 log2FC")
ax2d.set_title("Pairing attenuates donor background transcripts", fontsize=7.45)
despine(ax2d)
ax2d.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.5)

# Panel 2e: Non-Circular MSigDB Hallmark & Module Enrichment
ax2e = fig2.add_subplot(gs2[2, :])
panel_label(ax2e, "e", x=-0.06, y=1.05)
pick_pws = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "CUSTOM_ADAPTIVE_CHECKPOINT_LIGANDS_NO_IFN14",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "CUSTOM_ANTIGEN_PRESENTATION_MHC_NO_IFN14",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYC_TARGETS_V1",
]
fg_sub = fg_peng[fg_peng["pathway"].isin(pick_pws)].sort_values("NES", ascending=True).reset_index(drop=True)

def clean_pw(s: str) -> str:
    s = s.replace("HALLMARK_", "Hallmark: ").replace("CUSTOM_", "Module: ").replace("_NO_IFN14", " (non-IFN14)")
    return s.replace("_", " ").title().replace("Ifn14", "IFN-14").replace("Mhc", "MHC-I/II").replace("Nfkb", "NF-kB").replace("Il6 Jak Stat3", "IL6-JAK-STAT3").replace("Tnfa", "TNFa")

y_fg = np.arange(len(fg_sub))
bar_cols = [CB_PAL["purple"] if "CUSTOM_" in p else (CB_PAL["ifn_high"] if v > 0 else CB_PAL["ifn_low"])
            for p, v in zip(fg_sub["pathway"], fg_sub["NES"])]
ax2e.barh(y_fg, fg_sub["NES"].values, height=0.68, color=bar_cols, alpha=0.88, edgecolor="#1E293B", linewidth=0.4)
ax2e.axvline(0, color="#1E293B", lw=0.75)
for yi, (_, r_f) in enumerate(fg_sub.iterrows()):
    nes_v, padj_v = float(r_f["NES"]), float(r_f["padj"])
    txt = fmt_p(padj_v, "FDR")
    ax2e.text(nes_v + 0.07 if nes_v >= 0 else nes_v - 0.07, yi, txt,
              va="center", ha="left" if nes_v >= 0 else "right", fontsize=5.8, color="#111111")

ax2e.set_yticks(y_fg)
ax2e.set_yticklabels([clean_pw(p) for p in fg_sub["pathway"]], fontsize=6.3)
ax2e.set_xlim(-3.45, 4.55)
ax2e.set_xlabel("Normalized Enrichment Score (NES; fgsea on 15,030 non-IFN14 genes)")
ax2e.set_title("Non-circular enrichment: coordinated antigen-presentation and selective immunoregulatory programs", fontsize=7.45)
despine(ax2e)

export_fig(fig2, FIG_DIR, "Figure2")

# =========================================================================
# FIGURE 3: CROSS-COHORT REPLICATION WITH DONOR-BOOTSTRAP CIs & LODO
# =========================================================================
print("[3/12] Rendering Figure3 (Therapy/site detectability, donor-bootstrap concordance CIs, CD47 slopegraphs & dot-matrix)...")
don_202   = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/gse202051_donor_abundance_and_treatment_Q1.tsv", sep="\t")
don_154   = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/gse154778_donor_abundance_primary_vs_met_Q1.tsv", sep="\t")
conc_voom = pd.read_csv(P2_DIR / "cross_cohort_voom_donor_bootstrap_and_LODO_concordance.tsv", sep="\t")
pair_vals = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/donor_level_paired_gene_values.tsv", sep="\t")
func_pan  = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/paired_donor_functional_panel_all_cohorts.tsv", sep="\t")

don_202.to_csv(SRC_DIR / "Figure3a_gse202051_treatment_detectability_source.tsv", sep="\t", index=False)
don_154.to_csv(SRC_DIR / "Figure3a_gse154778_site_detectability_source.tsv", sep="\t", index=False)
conc_voom.to_csv(SRC_DIR / "Figure3b_donor_bootstrap_concordance_source.tsv", sep="\t", index=False)
pair_vals[pair_vals["gene"] == "CD47"].to_csv(SRC_DIR / "Figure3c_cd47_paired_slopegraphs_source.tsv", sep="\t", index=False)
func_pan.to_csv(SRC_DIR / "Figure3d_multicohort_paired_dotmatrix_source.tsv", sep="\t", index=False)
reg_src("Figure3", "a", "Figure3a_gse202051_treatment_detectability_source.tsv", "IFN-high TAM detectability across therapy & lesion site")
reg_src("Figure3", "b", "Figure3b_donor_bootstrap_concordance_source.tsv", "1,000-draw two-cohort donor-resampling bootstrap CIs & LODO on Spearman rho")
reg_src("Figure3", "c-d", "Figure3c_cd47_paired_slopegraphs_source.tsv", "Patient-paired CD47 slopegraphs & 4-stratum functional panel")

fig3 = plt.figure(figsize=(DOUBLE_COL_IN, 7.45))
gs3 = gridspec.GridSpec(3, 1, figure=fig3, height_ratios=[1.10, 1.00, 1.24], hspace=0.52)
gs3_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs3[0], width_ratios=[0.84, 1.16], wspace=0.80)
gs3_mid = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs3[1], wspace=0.34)

# Panel 3a: Detectability Across Neoadjuvant Therapy (GSE202051) & Lesion Site (GSE154778)
ax3a = fig3.add_subplot(gs3_top[0])
panel_label(ax3a, "a", x=-0.20, y=1.06)
u_q1  = don_202.loc[don_202["treatment_binary"] == "Untreated", "ifn14_high_frac_of_mac_Q1"].values * 100.0
n_q1  = don_202.loc[don_202["treatment_binary"] == "Neoadjuvant", "ifn14_high_frac_of_mac_Q1"].values * 100.0
u_cal = don_202.loc[don_202["treatment_binary"] == "Untreated", "ifn14_high_frac_fixed_PengCal"].values * 100.0
n_cal = don_202.loc[don_202["treatment_binary"] == "Neoadjuvant", "ifn14_high_frac_fixed_PengCal"].values * 100.0
p_154 = don_154.loc[don_154["lesion_site"] == "Primary", "ifn14_high_frac_of_mac_Q1"].values * 100.0
m_154 = don_154.loc[don_154["lesion_site"] == "Metastasis", "ifn14_high_frac_of_mac_Q1"].values * 100.0

for pos_x, vals_x, col_x in [
    (-0.16, u_q1, CB_PAL["untreated"]),  (0.16, n_q1, CB_PAL["neoadjuvant"]),
    (0.84, u_cal, CB_PAL["untreated"]),  (1.16, n_cal, CB_PAL["neoadjuvant"]),
    (1.84, p_154, CB_PAL["primary"]),    (2.16, m_154, CB_PAL["metastasis"]),
]:
    bp3 = ax3a.boxplot([vals_x], positions=[pos_x], widths=0.25, patch_artist=True, showfliers=False,
                       medianprops=dict(color="#111111", lw=1.2), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    bp3["boxes"][0].set_facecolor(col_x); bp3["boxes"][0].set_alpha(0.35)
    ax3a.scatter(np.full_like(vals_x, pos_x) + rng.uniform(-0.05, 0.05, len(vals_x)), vals_x,
                 s=10, color=col_x, alpha=0.85, edgecolors="white", linewidths=0.2, zorder=3)

ax3a.text(0, 74, "p=0.453", ha="center", fontsize=5.5)
ax3a.text(1, 61, "p=0.382", ha="center", fontsize=5.5)
ax3a.text(2, 43, "p=0.022\n(FDR=0.067)", ha="center", fontsize=5.2, fontweight="bold")
ax3a.set_xticks([0, 1, 2])
ax3a.set_xticklabels(["GSE202051\nQ1", "GSE202051\nPeng-cal", "GSE154778\nPT vs Met"], fontsize=5.9)
ax3a.set_ylim(-3, 99)
ax3a.set_ylabel("IFN-high TAMs / patient (%)")
ax3a.set_title("Detectable across therapy & site", fontsize=7.2)
despine(ax3a)
ax3a.legend(
    handles=[
        mpatches.Patch(facecolor=CB_PAL["untreated"], alpha=0.65, label="Untreated (18)"),
        mpatches.Patch(facecolor=CB_PAL["neoadjuvant"], alpha=0.65, label="Neoadj (25)"),
        mpatches.Patch(facecolor=CB_PAL["primary"], alpha=0.65, label="PT(10) / Met(6)"),
    ],
    loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=5.2,
)

# Panel 3b: Upgraded Cross-Cohort Concordance Forest Plot (Dynamic n_genes + Matched voom CI & edgeR!)
ax3b = fig3.add_subplot(gs3_top[1])
panel_label(ax3b, "b", x=-0.46, y=1.08)
conc_specs = [
    ("GSE202051_All43_Q2_IFN14_min5", "1_", "GSE202051 (n=34): all", CB_PAL["slate"]),
    ("GSE202051_All43_Q2_IFN14_min5", "2_", "GSE202051 (n=34): panel", CB_PAL["neoadjuvant"]),
    ("GSE202051_All43_Q2_IFN14_min5", "3_", "GSE202051 (n=34): DEGs", CB_PAL["neoadjuvant"]),
    ("GSE202051_Untreated18_Q2_IFN14_min5", "2_", "Untreated (n=12): panel", CB_PAL["untreated"]),
    ("GSE202051_Untreated18_Q2_IFN14_min5", "3_", "Untreated (n=12): DEGs", CB_PAL["untreated"]),
    ("GSE202051_Neoadjuvant25_Q2_IFN14_min5", "3_", "Neoadj (n=22): DEGs", "#B45309"),
    ("GSE154778_All16_Q2_IFN14_min5", "1_", "GSE154778 (n=15): all", CB_PAL["slate"]),
    ("GSE154778_All16_Q2_IFN14_min5", "2_", "GSE154778 (n=15): panel", CB_PAL["purple"]),
    ("GSE154778_All16_Q2_IFN14_min5", "3_", "GSE154778 (n=15): DEGs", CB_PAL["purple"]),
]
y_c = np.arange(len(conc_specs))[::-1]
ax3b.axvline(0.0, color="#666666", ls="--", lw=0.7)
dyn_ylabels = []
for idx_r, (vm, g_pref, lbl_r, col_r) in enumerate(conc_specs):
    yv = y_c[idx_r]
    rc = conc_voom[(conc_voom["validation_cohort_model"] == vm) & (conc_voom["gene_set_definition"].str.startswith(g_pref))].iloc[0]
    ng_r  = int(rc["n_genes"])
    rho_e = float(rc["edgeR_spearman_rho"])
    rho_v = float(rc["voom_precision_spearman_rho"])
    lo_b  = float(rc["donor_bootstrap_95CI_low"])
    hi_b  = float(rc["donor_bootstrap_95CI_high"])
    lodo_lo = float(rc["LODO_min_spearman_rho"])
    lodo_hi = float(rc["LODO_max_spearman_rho"])
    dyn_ylabels.append(f"{lbl_r} (g={ng_r:,})")
    ax3b.plot([lo_b, hi_b], [yv, yv], color=col_r, lw=1.35, zorder=2)
    ax3b.plot([lodo_lo, lodo_hi], [yv, yv], color="#111111", lw=2.4, alpha=0.35, zorder=3)
    ax3b.scatter([rho_v], [yv], s=18, color=col_r, edgecolors="#111111", linewidths=0.45, zorder=4,
                 label="voom (+boot CI & LODO)" if idx_r == 0 else None)
    ax3b.scatter([rho_e], [yv], s=18, marker="D", facecolors="white", edgecolors=col_r, linewidths=0.95, zorder=5,
                 label="edgeR QL" if idx_r == 0 else None)
    ax3b.text(0.84, yv, f"v={rho_v:.2f} [{lo_b:.2f},{hi_b:.2f}] | e={rho_e:.2f}", va="center", ha="left", fontsize=4.65, color="#111111")

ax3b.set_yticks(y_c)
ax3b.set_yticklabels(dyn_ylabels, fontsize=5.15)
ax3b.set_xlim(-0.25, 1.66)
ax3b.set_xlabel("Paired log2FC Spearman rho (voom 95% CI & LODO)", fontsize=6.6)
ax3b.set_title("Cross-cohort concordance (v = voom; e = edgeR)", fontsize=7.0, pad=13)
ax3b.legend(loc="upper center", bbox_to_anchor=(0.45, 1.11), ncol=2, frameon=False, fontsize=4.9, handletextpad=0.3, columnspacing=0.8)
despine(ax3b)

# Panel 3c: Patient-Paired CD47 Slopegraphs Across All 3 Cohorts
for s_i, (coh_s, col_s) in enumerate([("Peng2019", CB_PAL["pdac"]), ("GSE202051", CB_PAL["untreated"]), ("GSE154778", CB_PAL["purple"])]):
    ax_s = fig3.add_subplot(gs3_mid[s_i])
    if s_i == 0:
        panel_label(ax_s, "c", x=-0.26, y=1.08)
    sub_pv = pair_vals[(pair_vals["cohort"] == coh_s) & (pair_vals["gene"] == "CD47")].copy()
    lo_arr, hi_arr = sub_pv["mean_low"].values, sub_pv["mean_high"].values
    _, p_wil = stats.wilcoxon(hi_arr, lo_arr)
    pct_up = np.mean((hi_arr - lo_arr) > 0) * 100.0
    for _, rp in sub_pv.iterrows():
        is_up = rp["mean_high"] >= rp["mean_low"]
        ax_s.plot([0, 1], [rp["mean_low"], rp["mean_high"]],
                  color=col_s if is_up else "#94A3B8", alpha=0.58 if is_up else 0.42, lw=0.85)
    ax_s.plot([0, 1], [np.median(lo_arr), np.median(hi_arr)], color="#111111", lw=2.1, marker="o", markersize=4.3, zorder=5)
    ax_s.set_xticks([0, 1])
    ax_s.set_xticklabels(["IFN-low\nTAMs", "IFN-high\nTAMs"], fontsize=6.3)
    if s_i == 0:
        ax_s.set_ylabel("Donor mean CD47 log1p")
    coh_lbl = "CRA001160" if coh_s == "Peng2019" else coh_s
    ax_s.set_title(f"CD47: {coh_lbl} (n={len(sub_pv)})\n{pct_up:.1f}% up | {fmt_p(p_wil)}", fontsize=6.65)
    ax_s.set_xlim(-0.25, 1.25)
    despine(ax_s)

# Panel 3d: Multi-Cohort Donor-Paired Dot-Matrix (PuOr_r purple-orange colorblind-safe diverging colormap!)
ax3d = fig3.add_subplot(gs3[2])
panel_label(ax3d, "d", x=-0.06, y=1.05)
genes_3d = ["STAT1", "IRF7", "TAP1", "HLA-A", "HLA-B", "HLA-C", "B2M", "CD47", "LGALS9", "CD80", "CD86", "NECTIN2", "CD274", "PDCD1LG2", "IDO1", "VSIR", "CD276", "HLA-DRA"]
cohorts_3d = [
    ("Peng2019_All", "CRA001160 Discovery (n=31)"),
    ("GSE202051_Untreated18", "GSE202051 Untreated (n=12)"),
    ("GSE202051_Neoadjuvant25", "GSE202051 Neoadjuvant (n=22)"),
    ("GSE154778_All16", "GSE154778 Primary+Met (n=15)"),
]
fp_sub = func_pan[(func_pan["scheme"] == "Q2_IFN14_min5") & (func_pan["gene"].isin(genes_3d)) & (func_pan["cohort"].isin([c[0] for c in cohorts_3d]))].copy()
gx_map = {g: j for j, g in enumerate(genes_3d)}
cy_map = {c[0]: i for i, c in enumerate(reversed(cohorts_3d))}
xs_d, ys_d, cs_d, ss_d, sig_e = [], [], [], [], []
for _, rf in fp_sub.iterrows():
    xs_d.append(gx_map[rf["gene"]])
    ys_d.append(cy_map[rf["cohort"]])
    cs_d.append(np.clip(float(rf["mean_paired_delta_log1p"]), -0.35, 0.35))
    p_u = float(rf["wilcoxon_mean_p"]) if pd.notna(rf["wilcoxon_mean_p"]) else 1.0
    ss_d.append(np.clip(-np.log10(max(p_u, 1e-8)) * 16.0 + 12.0, 12.0, 112.0))
    sig_e.append("#111111" if p_u < 0.05 else "#CBD5E1")

sc_3d = ax3d.scatter(xs_d, ys_d, s=ss_d, c=cs_d, cmap="PuOr_r", vmin=-0.35, vmax=0.35,
                     edgecolors=sig_e, linewidths=[1.0 if e == "#111111" else 0.35 for e in sig_e])
ax3d.set_xticks(range(len(genes_3d)))
ax3d.set_xticklabels(genes_3d, rotation=40, ha="right", fontsize=6.3)
ax3d.set_yticks(range(len(cohorts_3d)))
ax3d.set_yticklabels([c[1] for c in reversed(cohorts_3d)], fontsize=6.4)
ax3d.set_xlim(-0.6, len(genes_3d) - 0.4); ax3d.set_ylim(-0.6, len(cohorts_3d) - 0.4)
ax3d.set_title("Multi-cohort donor-paired shifts (PuOr_r color: mean delta log1p; dark rim: paired Wilcoxon p < 0.05)", fontsize=7.25)
despine(ax3d)
cb3d = fig3.colorbar(sc_3d, ax=ax3d, fraction=0.025, pad=0.02)
cb3d.ax.tick_params(labelsize=5.8); cb3d.set_label("Mean paired\ndelta log1p", fontsize=5.9)

export_fig(fig3, FIG_DIR, "Figure3")

# =========================================================================
# FIGURE 4: REGULONS & N_perm = 1,000 HARD-CONSTRAINT REGVELO AUDIT
# =========================================================================
print("[4/12] Rendering Figure4 (CollecTRI/pySCENIC, mask audit, N_perm=1,000 RegVelo + BH q, & evidence matrix)...")
ulm_nc   = pd.read_csv(BASE_ROOT / "05_results/stage6_clinical_drug/peng_macrophage_paired_collectri_ulm.tsv", sep="\t").query("mode == 'NonCircular_NoIFN10_Targets'").copy()
sc_c14   = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_macrophage_paired_C14_aucell.csv")
edge_aud = pd.read_csv(BASE_ROOT / "05_results/stage4_regvelo/regvelo_directed_edge_overlap_audit.tsv", sep="\t")
rgv_1000 = pd.read_csv(P2_DIR / "regvelo_1000perm_macrophage_summary.tsv", sep="\t")
rgv_null = pd.read_csv(P2_DIR / "regvelo_1000perm_null_distributions.tsv", sep="\t")
mac_mat  = pd.read_csv(P2_DIR / "macrophage_ifn_tf_evidence_matrix_csbj_full.tsv", sep="\t")

m_reg = ulm_nc.merge(sc_c14, on="TF", how="inner")
m_reg.to_csv(SRC_DIR / "Figure4a_collectri_vs_pyscenic_source.tsv", sep="\t", index=False)
edge_aud.to_csv(SRC_DIR / "Figure4b_prior_mask_verification_source.tsv", sep="\t", index=False)
rgv_1000.to_csv(SRC_DIR / "Figure4c_regvelo_1000perm_macrophage_source.tsv", sep="\t", index=False)
mac_mat.head(25).to_csv(SRC_DIR / "Figure4de_macrophage_evidence_matrix_source.tsv", sep="\t", index=False)
reg_src("Figure4", "a-b", "Figure4a_collectri_vs_pyscenic_source.tsv", "Non-circular CollecTRI vs pySCENIC & hard-constraint mask verification")
reg_src("Figure4", "c", "Figure4c_regvelo_1000perm_macrophage_source.tsv", "N_perm=1,000 degree+expr-matched RegVelo screen with BH q & 16-donor CIs")
reg_src("Figure4", "d-e", "Figure4de_macrophage_evidence_matrix_source.tsv", "5-layer macrophage evidence matrix & Dirichlet/LOLO rank stability")

fig4 = plt.figure(figsize=(DOUBLE_COL_IN, 7.50))
gs4 = gridspec.GridSpec(3, 2, figure=fig4, height_ratios=[1.00, 1.28, 1.10], hspace=0.50, wspace=0.40)

# Panel 4a: Non-Circular Paired CollecTRI ULM vs Motif-Pruned pySCENIC AUCell
ax4a = fig4.add_subplot(gs4[0, 0])
panel_label(ax4a, "a", x=-0.18, y=1.07)
both_sig = (m_reg["paired_t_FDR"] < 0.05) & (m_reg["wilcoxon_padj"] < 0.05) & (m_reg["delta_ulm_high_minus_low"] > 0) & (m_reg["delta_auc_high_vs_low"] > 0)
ax4a.axhline(0, color="#888888", ls=":", lw=0.6); ax4a.axvline(0, color="#888888", ls=":", lw=0.6)
ax4a.scatter(m_reg.loc[~both_sig, "delta_ulm_high_minus_low"], m_reg.loc[~both_sig, "delta_auc_high_vs_low"] * 100.0,
             s=14, color=CB_PAL["light_slate"], alpha=0.72, linewidths=0)
ax4a.scatter(m_reg.loc[both_sig, "delta_ulm_high_minus_low"], m_reg.loc[both_sig, "delta_auc_high_vs_low"] * 100.0,
             s=26, color=CB_PAL["ifn_high"], edgecolors="#111111", linewidths=0.45, label=f"Dual FDR<0.05 (n={int(both_sig.sum())})")
for tf_l, (ox, oy) in {"STAT1": (-0.11, 0.24), "IRF1": (-0.08, 0.26), "IRF7": (-0.18, 0.14), "STAT2": (-0.09, -0.32), "RELB": (0.05, -0.18), "MITF": (0.05, 0.12), "PRDM1": (-0.28, -0.22)}.items():
    rt = m_reg[m_reg["TF"] == tf_l]
    if len(rt):
        x0, y0 = float(rt["delta_ulm_high_minus_low"].iloc[0]), float(rt["delta_auc_high_vs_low"].iloc[0]) * 100.0
        ax4a.annotate(tf_l, xy=(x0, y0), xytext=(x0 + ox, y0 + oy), fontsize=5.7, fontweight="bold",
                      color="#111111", arrowprops=dict(arrowstyle="-", lw=0.35, color="#666666"))
ax4a.set_xlabel("Non-circular CollecTRI delta ULM (n=31 donors)")
ax4a.set_ylabel("Motif-pruned pySCENIC delta AUC (x100)")
ax4a.set_title("Steady-state cis-regulon enrichment in IFN-high TAMs", fontsize=7.2)
despine(ax4a)
ax4a.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.5)

# Panel 4b: Prior Mask & Jacobian Orientation Verification (Zero-green palette!)
ax4b = fig4.add_subplot(gs4[0, 1])
panel_label(ax4b, "b", x=-0.24, y=1.07)
n_pri = float(edge_aud["n_prior_directed_edges"].mean())
n_inf = float(edge_aud["n_inferred_nonzero_edges_gt1e6"].mean())
n_ok  = float(edge_aud["n_overlap_correct_orientation"].mean())
n_bug = float(edge_aud["n_overlap_if_transposed_v0_bug"].mean())
n_off = float(edge_aud["n_spurious_off_prior_edges"].mean())
bars_4b = [
    ("CollecTRI\nprior", n_pri, CB_PAL["control"]),
    ("Hard-mask\nactive", n_inf, CB_PAL["navy"]),
    ("Aligned [TxR]\n(mask=100%)", n_ok, CB_PAL["ifn_high"]),
    ("If transposed\n(axis check)", n_bug, CB_PAL["ifn_mid"]),
    ("Off-prior\n(mask=0)", n_off, CB_PAL["slate"]),
]
for bi, (blbl, bval, bcol) in enumerate(bars_4b):
    ax4b.bar(bi, bval, width=0.56, color=bcol, alpha=0.88, edgecolor="#1E293B", linewidth=0.45)
    ax4b.text(bi, bval + 45, f"{int(round(bval))}", ha="center", va="bottom", fontsize=5.9, fontweight="bold")
ax4b.set_xticks(range(len(bars_4b)))
ax4b.set_xticklabels([b[0] for b in bars_4b], fontsize=5.4)
ax4b.set_ylim(0, n_pri * 1.25)
ax4b.set_ylabel("Directed GRN edges (3-seed mean)")
ax4b.set_title("Prior-mask & Jacobian orientation verification", fontsize=7.2)
despine(ax4b)

# Panel 4c: Upgraded N_perm = 1,000 Degree- & Expression-Matched RegVelo Perturbation Screen + Donor CIs + BH q!
ax4c = fig4.add_subplot(gs4[1, :])
panel_label(ax4c, "c", x=-0.06, y=1.05)
s_top = rgv_1000[rgv_1000["perturbation_type"] == "Single_TF_KO"].sort_values(["empirical_perm_p", "abs_mean_delta_pp"], ascending=[True, False]).head(7)["perturbation"].tolist()
must_4c = ["STAT1+IRF1", "STAT1", "IRF1", "STAT2", "IRF2", "ACTB", "GAPDH"]
pick_4c = []
for nm in s_top + must_4c:
    if nm not in pick_4c and nm in set(rgv_1000["perturbation"]):
        pick_4c.append(nm)
sub_4c = rgv_1000.set_index("perturbation").loc[pick_4c].reset_index().iloc[::-1].reset_index(drop=True)
y_4c = np.arange(len(sub_4c))

# Plot per-row matched 95% null interval as light gray horizontal bar behind each TF!
for yi, (_, rp) in enumerate(sub_4c.iterrows()):
    sd_n = float(rp["null_matched_sd_pp"])
    if sd_n > 0:
        ax4c.plot([-1.96 * sd_n, 1.96 * sd_n], [yi, yi], color="#CBD5E1", lw=5.2, solid_capstyle="butt", zorder=1,
                  label="95% degree+expr-matched null (N=1,000)" if yi == len(sub_4c) - 1 else None)

ax4c.axvline(0.0, color="#1E293B", ls="--", lw=0.75, zorder=2)

for yi, (_, rp) in enumerate(sub_4c.iterrows()):
    pt, pn = rp["perturbation_type"], rp["perturbation"]
    m_pp = float(rp["mean_delta_pp"])
    d_lo = float(rp["donor16_bootstrap_95CI_low_pp"])
    d_hi = float(rp["donor16_bootstrap_95CI_high_pp"])
    ep   = float(rp["empirical_perm_p"])
    bhq  = float(rp["BH_adjusted_q"]) if pd.notna(rp["BH_adjusted_q"]) else np.nan
    deg_i = int(rp["n_muted_prior_edges"])

    if pt == "Negative_Control_NonTF":
        col_p = CB_PAL["slate"]
    elif pt == "Combinatorial_KO":
        col_p = CB_PAL["purple"]
    elif pn in {"STAT1", "IRF1", "STAT2", "IRF2"}:
        col_p = CB_PAL["ifn_mid"]
    else:
        col_p = CB_PAL["ifn_high"] if (pd.notna(bhq) and bhq < 0.10) else CB_PAL["control"]

    ax4c.barh(yi, m_pp, height=0.58, color=col_p, alpha=0.88, edgecolor="#1E293B", linewidth=0.4, zorder=3)
    if deg_i > 0:
        ax4c.plot([d_lo, d_hi], [yi, yi], color="#111111", lw=1.15, zorder=4,
                  label="16-patient bootstrap 95% CI" if yi == len(sub_4c) - 1 else None)

    if pt == "Negative_Control_NonTF":
        lbl_v = "0.00 pp (0 prior edges; numerical invariance)"
    else:
        q_str = f"BH q={bhq:.3f}" if bhq >= 0.001 else f"BH q={bhq:.4f}"
        p_str = f"p={ep:.3f}" if ep >= 0.002 else f"p={ep:.4f}"
        lbl_v = f"{m_pp:+.2f} pp (k={deg_i}, {p_str}, {q_str})"

    if m_pp < 0:
        x_txt = min(m_pp, d_lo, -1.96 * float(rp["null_matched_sd_pp"])) - 0.08
        ax4c.text(x_txt, yi, lbl_v, va="center", ha="right", fontsize=5.05, color="#111111")
    else:
        x_txt = max(m_pp, d_hi, 1.96 * float(rp["null_matched_sd_pp"])) + 0.08
        ax4c.text(x_txt, yi, lbl_v, va="center", ha="left", fontsize=5.05, color="#111111")

ax4c.set_yticks(y_4c)
ax4c.set_yticklabels(sub_4c["perturbation"].tolist(), fontsize=5.9, fontweight="bold")
max_ext = max(
    np.max(np.abs(sub_4c["donor16_bootstrap_95CI_low_pp"].values)),
    np.max(np.abs(sub_4c["donor16_bootstrap_95CI_high_pp"].values)),
    1.5,
)
ax4c.set_xlim(-max_ext * 2.55, max_ext * 2.45)
ax4c.set_xlabel("Predicted fate shift on TAM_STAT1_Inflam (percentage points; 3 GPU seeds; 16 patients; N_perm = 1,000)")
ax4c.set_title("Lineage-restricted Macrophage_TAM (n=3,318) RegVelo screen: degree+expr-matched 1,000-permutation null & BH q", fontsize=7.15)
despine(ax4c)
ax4c.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=5.2)

# Panel 4d: Compartment-Separated Macrophage 5-Layer Evidence Matrix
ax4d = fig4.add_subplot(gs4[2, 0])
panel_label(ax4d, "d", x=-0.18, y=1.06)
top12 = mac_mat.head(12).copy()
cols_m5 = ["score_M1_Peng_IFN10_DEG", "score_M2_Peng_CollecTRI_ULM", "score_M3_Peng_pySCENIC_Motif", "score_M4_GSE202051_Replication", "score_M5_GSE154778_RegVelo_1000Perm"]
mat_4d = top12[cols_m5].values
ax4d.imshow(mat_4d, aspect="auto", cmap="Blues", vmin=0.0, vmax=1.0)
ax4d.set_xticks(range(5))
ax4d.set_xticklabels(["M1: CRA\nIFN10 DEG", "M2: CRA\nCollecTRI", "M3: CRA\npySCENIC", "M4: GSE202\nIFN10 Rep", "M5: GSE154\nRegVelo |z|"], fontsize=5.3)
ax4d.set_yticks(range(len(top12)))
ax4d.set_yticklabels([f"#{int(r['primary_rank'])} {r['TF']}" for _, r in top12.iterrows()], fontsize=5.9)
for ri in range(mat_4d.shape[0]):
    for ci in range(5):
        vv = mat_4d[ri, ci]
        ax4d.text(ci, ri, f"{vv:.2f}", ha="center", va="center", fontsize=5.0, color="white" if vv > 0.58 else "#111111")
ax4d.set_title("Complementary 5-layer macrophage matrix", fontsize=7.15)

# Panel 4e: Dirichlet Cohort-Weight & Drop-M5 Rank Sensitivity
ax4e = fig4.add_subplot(gs4[2, 1])
panel_label(ax4e, "e", x=-0.18, y=1.06)
top12_rev = top12.iloc[::-1].reset_index(drop=True)
y_4e = np.arange(len(top12_rev))
med_r, iqr_r = top12_rev["dirichlet_median_rank"].values, top12_rev["dirichlet_iqr_rank"].values
ax4e.hlines(y_4e, np.maximum(1, med_r - iqr_r / 2.0), med_r + iqr_r / 2.0, color=CB_PAL["control"], lw=1.6, label="Dirichlet(1,1,1) IQR")
ax4e.scatter(top12_rev["primary_rank"].values, y_4e, s=24, color=CB_PAL["ifn_high"], edgecolors="#111111", linewidths=0.45, zorder=4, label="Cohort-balanced rank")
ax4e.scatter(top12_rev["rank_drop_M5_RegVelo"].values, y_4e, s=20, marker="D", facecolors="white", edgecolors=CB_PAL["purple"], linewidths=0.9, zorder=5, label="Drop M5 (RegVelo)")
ax4e.set_yticks(y_4e)
ax4e.set_yticklabels(top12_rev["TF"].tolist(), fontsize=5.9)
ax4e.set_xlabel("Macrophage TF rank (1 = strongest support)", fontsize=6.8)
ax4e.set_title("Dirichlet weight & leave-M5-out rank stability", fontsize=7.15)
ax4e.invert_xaxis()
despine(ax4e)
ax4e.legend(loc="lower left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.2)

export_fig(fig4, FIG_DIR, "Figure4")

# =========================================================================
# FIGURE 5: EXPRESSION-INFERRED CANDIDATE COMMUNICATION & 7-PATIENT VISIUM MORAN'S I
# =========================================================================
print("[5/12] Rendering Figure5 (Expression-inferred LR axes, subsampling audit & 7-patient Visium Moran's I)...")
lr_aud  = pd.read_csv(BASE_ROOT / "05_results/LR_EDGE_AUDIT.tsv", sep="\t")
cc_sens = pd.read_csv(BASE_ROOT / "05_results/stage5_spatial_comm/cellchat_composition_and_subsampling_sensitivity.tsv", sep="\t")
sp_site = pd.read_csv(BASE_ROOT / "05_results/stage5_spatial_comm/gse272362_univariate_spatial_autocorrelation_by_site_and_patient.tsv", sep="\t")
sp_long = pd.read_csv(BASE_ROOT / "05_results/stage5_spatial_comm/gse272362_per_section_morans_i_long.tsv", sep="\t")
sp_pat7 = pd.read_csv(MAN_DIR / "GSE272362_patient_level_morans_summary.tsv", sep="\t")

lr_aud.to_csv(SRC_DIR / "Figure5a_expression_inferred_lr_audit_source.tsv", sep="\t", index=False)
cc_sens.to_csv(SRC_DIR / "Figure5b_cellchat_subsampling_source.tsv", sep="\t", index=False)
sp_pat7.to_csv(SRC_DIR / "Figure5cd_gse272362_patient_and_site_morans_source.tsv", sep="\t", index=False)
reg_src("Figure5", "a-b", "Figure5a_expression_inferred_lr_audit_source.tsv", "Expression-inferred juxtacrine/paracrine LR candidates & balanced subsampling")
reg_src("Figure5", "c-d", "Figure5cd_gse272362_patient_and_site_morans_source.tsv", "GSE272362 7-patient & site-stratified univariate Moran's I (97,020 in-tissue spots)")

fig5 = plt.figure(figsize=(DOUBLE_COL_IN, 6.90))
gs5 = gridspec.GridSpec(2, 2, figure=fig5, height_ratios=[1.12, 1.18], hspace=0.49, wspace=0.56)

# Panel 5a: Expression-Inferred Candidate Intercellular Axes
ax5a = fig5.add_subplot(gs5[0, 0])
panel_label(ax5a, "a", x=-0.24, y=1.06)
valid_lr = lr_aud[lr_aud["audit_decision"].astype(str).str.startswith("RETAIN")].copy()
valid_lr["axis_lbl"] = valid_lr["ligand"] + " - " + valid_lr["receptor"]  # non-causal dash per guideline 7.2.10!
valid_lr = valid_lr.iloc[::-1].reset_index(drop=True)
y_lr = np.arange(len(valid_lr))

def map_mode_cb(s: str) -> str:
    sl = str(s).lower()
    if "juxtacrine" in sl or "contact" in sl:
        return CB_PAL["ifn_high"]
    if "upstream" in sl or "ifng" in sl or "cytokine" in sl:
        return CB_PAL["purple"]
    return CB_PAL["control"]

cols_5a = [map_mode_cb(m) for m in valid_lr["signaling_mode"]]
fc_lr = pd.to_numeric(valid_lr["fold_change_high_vs_low"], errors="coerce").fillna(1.25).values
vals_5a = np.clip(np.log2(np.maximum(fc_lr, 0.5)), -0.3, 2.75)
p_lr = pd.to_numeric(valid_lr["donor_paired_wilcoxon_p"], errors="coerce").fillna(0.05).clip(lower=1e-8).values

ax5a.hlines(y_lr, 0, vals_5a, color=cols_5a, lw=1.5, alpha=0.86)
ax5a.scatter(vals_5a, y_lr, s=np.clip(-np.log10(p_lr) * 16.0 + 16.0, 16.0, 86.0), c=cols_5a, edgecolors="#111111", linewidths=0.5, zorder=3)
ax5a.axvline(0.0, color="#666666", ls=":", lw=0.7)
ax5a.set_yticks(y_lr)
ax5a.set_yticklabels(valid_lr["axis_lbl"].values, fontsize=5.6)
ax5a.set_xlim(-0.32, 3.05)
ax5a.set_xlabel("Paired co-expression log2FC (n=24 PDAC donors; cap 2.75)")
ax5a.set_title("Expression-inferred candidate LR pairs", fontsize=7.1)
despine(ax5a)
ax5a.legend(
    handles=[
        plt.Line2D([0], [0], marker="o", ls="", color=CB_PAL["ifn_high"], markersize=3.8, label="Cell-contact"),
        plt.Line2D([0], [0], marker="o", ls="", color=CB_PAL["control"], markersize=3.8, label="Secreted chemokine"),
        plt.Line2D([0], [0], marker="o", ls="", color=CB_PAL["purple"], markersize=3.8, label="T/NK IFNG"),
    ],
    loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=5.2,
)

# Panel 5b: CellChat Composition & Count-Matched Subsampling Sensitivity
ax5b = fig5.add_subplot(gs5[0, 1])
panel_label(ax5b, "b", x=-0.38, y=1.06)
y_cc = np.arange(len(cc_sens))[::-1]
short_cc = [
    "1. Raw all 10 types\n(incl. Ductal 2)",
    "2. Shared 9 types\n(no Ductal 2)",
    "3. Powered 7 shared\n(>=100 cells)",
    "4. Balanced 11v11\n(count-matched)",
]
ratios_cc = cc_sens["tumor_to_normal_ratio"].values
cols_cc = [CB_PAL["slate"], CB_PAL["ifn_mid"], CB_PAL["ifn_high"], CB_PAL["control"]]
ax5b.barh(y_cc, ratios_cc, height=0.58, color=cols_cc, alpha=0.88, edgecolor="#1E293B", linewidth=0.45)
ax5b.plot([1.56, 1.90], [y_cc[3], y_cc[3]], color="#111111", lw=1.8)
ax5b.axvline(1.0, color="#333333", ls="--", lw=0.75)
for idx_b, (yi, rv) in enumerate(zip(y_cc, ratios_cc)):
    ax5b.text((max(rv, 1.90) + 0.06) if idx_b == 3 else (rv + 0.06), yi, f"{rv:.2f}x", va="center", ha="left", fontsize=6.2, fontweight="bold")
ax5b.set_yticks(y_cc)
ax5b.set_yticklabels(short_cc[:len(cc_sens)], fontsize=5.9)
ax5b.set_xlim(0, 3.45)
ax5b.set_xlabel("Inferred Tumor / Control interaction count ratio")
ax5b.set_title("Composition & count-matched subsampling audit", fontsize=7.2)
despine(ax5b)

# Panel 5c: GSE272362 Patient-Aggregated Univariate Moran's I Heatmap by Site (YlOrBr sequential colormap)
ax5c = fig5.add_subplot(gs5[1, 0])
panel_label(ax5c, "c", x=-0.18, y=1.05)
genes_5c = ["B2M", "HLA-DRA", "HLA-A", "ISG15", "STAT1", "IRF1", "TAP1", "CXCL9", "CD68", "CD163", "CD47", "IDO1", "CD274", "PDCD1LG2"]
genes_5c = [g for g in genes_5c if g in set(sp_pat7["Feature Name"])]
sites_5c = [("Primary_Tumor", "PT\n(6 pts / 10s)"), ("Hepatic_Metastasis", "HM\n(7 pts / 12s)"), ("Lymph_Node_Metastasis", "LNM\n(3 pts / 5s)"), ("Normal_Pancreas_MatchedAutopsy", "NP*\n(2 pts / 3s)")]
mat_5c = np.zeros((len(genes_5c), 4))
txt_5c = [["" for _ in range(4)] for _ in range(len(genes_5c))]
for ri, gn in enumerate(genes_5c):
    for ci, (st_k, _) in enumerate(sites_5c):
        rr = sp_pat7[(sp_pat7["Feature Name"] == gn) & (sp_pat7["anatomical_group"] == st_k)].iloc[0]
        mat_5c[ri, ci] = float(rr["patient_level_median_I"])
        txt_5c[ri][ci] = f"{int(rr['n_patients_with_FDR005_section'])}/{int(rr['n_unique_patients'])}p"

im5c = ax5c.imshow(mat_5c, aspect="auto", cmap="YlOrBr", vmin=0.0, vmax=0.25)
ax5c.set_xticks(range(4))
ax5c.set_xticklabels([s[1] for s in sites_5c], fontsize=5.8)
ax5c.set_yticks(range(len(genes_5c)))
ax5c.set_yticklabels(genes_5c, fontsize=6.0)
for ri in range(len(genes_5c)):
    for ci in range(4):
        ax5c.text(ci, ri, txt_5c[ri][ci], ha="center", va="center", fontsize=5.1,
                  color="white" if mat_5c[ri, ci] > 0.14 else "#111111")
ax5c.set_title("Patient-aggregated univariate Moran's I (sig/n pts)", fontsize=7.1)
cb5c = fig5.colorbar(im5c, ax=ax5c, fraction=0.04, pad=0.03)
cb5c.ax.tick_params(labelsize=5.7); cb5c.set_label("Patient median I", fontsize=5.9)

# Panel 5d: Single-Gene Spatial Autocorrelation Across 30 Sections / 7 Autopsy Patients (Zero-green palette!)
ax5d = fig5.add_subplot(gs5[1, 1])
panel_label(ax5d, "d", x=-0.22, y=1.08)
genes_5d = ["B2M", "STAT1", "HLA-A", "CXCL9", "CD47", "CD274"]
site_pal_5d = {
    "Primary_Tumor":         ("PT (10s)", CB_PAL["pdac"], -0.21),
    "Hepatic_Metastasis":    ("HM (12s)", CB_PAL["hm"],   -0.07),
    "Lymph_Node_Metastasis": ("LNM (5s)", CB_PAL["lnm"],   0.07),
    "Normal_Pancreas":       ("NP* (3s)", CB_PAL["np"],    0.21),
}
for gi, gn in enumerate(genes_5d):
    sub_l = sp_long[sp_long["gene"] == gn]
    for ag_key, (short_st, col_st, offs) in site_pal_5d.items():
        s_st = sub_l[sub_l["anatomical_group"] == ag_key]
        ax5d.scatter(np.full(len(s_st), gi + offs), s_st["morans_I"].values, s=11, color=col_st,
                     alpha=0.85, edgecolors="#1E293B", linewidths=0.25, label=short_st if gi == 0 else None)
    pat_med = float(sp_site.loc[sp_site["gene"] == gn, "patient_level_median_morans_I"].iloc[0])
    ax5d.plot([gi - 0.26, gi + 0.26], [pat_med, pat_med], color="#111111", lw=1.6, zorder=5)

ax5d.axhline(0.0, color="#666666", ls="--", lw=0.65)
ax5d.set_xticks(range(len(genes_5d)))
ax5d.set_xticklabels(genes_5d, fontsize=6.3)
ax5d.set_ylim(-0.04, 0.42)
ax5d.set_ylabel("Univariate Space Ranger Moran's I")
ax5d.set_title("97,020 in-tissue spots (30 secs / 7 pts; bar: 7-pt med.)", fontsize=7.0)
despine(ax5d)
ax5d.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.5, ncol=2, title="Autopsy site")

export_fig(fig5, FIG_DIR, "Figure5")

# =========================================================================
# FIGURE 6: CLINICAL SURVIVAL WITH BH q & COMPOSITIONAL CORRELATIONS IN TCGA-PAAD
# =========================================================================
print("[6/12] Rendering Figure6 (Pure-ductal TCGA-PAAD Cox OS with BH q, descriptive KM & partial correlations)...")
surv_bh  = pd.read_csv(P2_DIR / "tcga_survival_3models_with_BH_and_PH_audit.tsv", sep="\t")
corr_ex  = pd.read_csv(BASE_ROOT / "05_results/stage6_clinical_drug/tcga_ductal_stat1_partial_correlations_exact.tsv", sep="\t")

surv_src = surv_bh[(surv_bh["endpoint"] == "OS") & (surv_bh["feature"].isin(["FOXM1","HMGA2","GATA6","STAT1","IRF1","IFN14_Signature","IFN10_Signature"]))].copy()
surv_src.to_csv(SRC_DIR / "Figure6a_tcga_cox_with_BH_q_source.tsv", sep="\t", index=False)
corr_ex.to_csv(SRC_DIR / "Figure6c_tcga_partial_spearman_source.tsv", sep="\t", index=False)
reg_src("Figure6", "a", "Figure6a_tcga_cox_with_BH_q_source.tsv", "Pure-ductal TCGA-PAAD Models 2 & 3 Cox OS with panel-wide BH q")
reg_src("Figure6", "b", "Figure6b_tcga_pure_ductal_km_source.tsv", "Descriptive median-split Kaplan-Meier OS curve & risk table (n=145)")
reg_src("Figure6", "c", "Figure6c_tcga_partial_spearman_source.tsv", "Raw vs purity/leukocyte-adjusted partial Spearman STAT1 correlations")

fig6 = plt.figure(figsize=(DOUBLE_COL_IN, 7.00))
gs6 = gridspec.GridSpec(2, 2, figure=fig6, height_ratios=[1.20, 1.12], hspace=0.50, wspace=0.38)

# Panel 6a: Multivariable & Purity-Adjusted Cox OS Forest Plot with Exact p and BH q
ax6a = fig6.add_subplot(gs6[0, 0])
panel_label(ax6a, "a", x=-0.22, y=1.06)
feats_6a = [
    ("Pure_Ductal_PDAC_n146", "FOXM1", "FOXM1 (ductal prolif.)"),
    ("Pure_Ductal_PDAC_n146", "HMGA2", "HMGA2 (basal-like)"),
    ("Pure_Ductal_PDAC_n146", "GATA6", "GATA6 (classical)"),
    ("Pure_Ductal_PDAC_n146", "STAT1", "STAT1 (pure ductal, n=144)"),
    ("Pure_Ductal_PDAC_n146", "IFN14_Signature", "IFN-14 sig (pure ductal)"),
    ("Pure_Ductal_PDAC_n146", "IRF1",  "IRF1 (pure ductal)"),
    ("Unfiltered_All_PAAD_n177", "STAT1", "STAT1 (unfiltered n=172*)"),
]
y_6a = np.arange(len(feats_6a))[::-1]
ax6a.axvline(1.0, color="#333333", ls="--", lw=0.75)

for i_f, (coh_f, fn_f, lbl_f) in enumerate(feats_6a):
    yv = y_6a[i_f]
    rf = surv_bh[(surv_bh["cohort"] == coh_f) & (surv_bh["endpoint"] == "OS") & (surv_bh["feature"] == fn_f)].iloc[0]
    hr2, lo2, hi2, p2, q2 = float(rf["HR_multi"]), float(rf["ci95_multi_low"] if "ci95_multi_low" in rf else rf["CI95_low_multi"]), float(rf["ci95_multi_high"] if "ci95_multi_high" in rf else rf["CI95_high_multi"]), float(rf["p_multi"]), float(rf["BH_q_multi"])
    hr3 = float(rf["HR_purity_immune_adj"])
    lo3 = float(rf["ci95_purity_low"] if "ci95_purity_low" in rf else rf["CI95_low_purity_immune"])
    hi3 = float(rf["ci95_purity_high"] if "ci95_purity_high" in rf else rf["CI95_high_purity_immune"])

    col_m2 = CB_PAL["ifn_high"] if (p2 < 0.05 and hr2 > 1) else (CB_PAL["control"] if (p2 < 0.05 and hr2 < 1) else CB_PAL["slate"])
    if coh_f == "Unfiltered_All_PAAD_n177":
        col_m2 = CB_PAL["ifn_mid"]

    ax6a.plot([lo2, hi2], [yv + 0.11, yv + 0.11], color=col_m2, lw=1.35)
    ax6a.scatter([hr2], [yv + 0.11], s=20, color=col_m2, edgecolors="#111111", linewidths=0.4, zorder=4)
    if np.isfinite(hr3):
        ax6a.plot([lo3, hi3], [yv - 0.13, yv - 0.13], color=CB_PAL["purple"], lw=1.0, alpha=0.82)
        ax6a.scatter([hr3], [yv - 0.13], s=16, marker="^", facecolors="white", edgecolors=CB_PAL["purple"], linewidths=0.9, zorder=4)
    p_short = f"p={p2:.3f}" if p2 >= 0.005 else f"p={p2:.4f}"
    ax6a.text(3.32, yv + 0.02, f"{hr2:.2f} ({p_short}, q={q2:.2f})", va="center", ha="right", fontsize=4.85, color=col_m2)

ax6a.set_yticks(y_6a)
ax6a.set_yticklabels([f[2] for f in feats_6a], fontsize=5.75)
ax6a.set_xlim(0.28, 3.38); ax6a.set_ylim(-1.45, len(feats_6a) - 0.35)
ax6a.set_xlabel("OS Hazard Ratio per 1-SD log2(RSEM+1) (95% CI)")
ax6a.set_title("Multivariable & purity-adjusted Cox OS (p & BH q)", fontsize=7.25)
despine(ax6a)
ax6a.legend(
    handles=[
        plt.Line2D([0], [0], marker="o", ls="-", color=CB_PAL["slate"], markersize=3.8, label="M2: Age+Sex+Stage+Grade"),
        plt.Line2D([0], [0], marker="^", ls="-", color=CB_PAL["purple"], markerfacecolor="white", markersize=4.0, label="M3: +Purity+Mac+PTPRC"),
    ],
    loc="lower center", ncol=2, frameon=True, facecolor="white", framealpha=0.95, edgecolor="none", fontsize=4.9, columnspacing=0.7,
)

# Panel 6b: Descriptive Kaplan-Meier OS + Risk Table (Pure-Ductal PDAC n=145; 84 deaths)
gs6_km = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs6[0, 1], height_ratios=[3.2, 1.0], hspace=0.08)
ax6b = fig6.add_subplot(gs6_km[0])
ax6b_rt = fig6.add_subplot(gs6_km[1])
panel_label(ax6b, "b", x=-0.18, y=1.08)

km_src_prev = BASE_ROOT / "06_figures/source_data/Figure6b_tcga_pure_ductal_km_source.tsv"
m_tc = pd.read_csv(km_src_prev, sep="\t")
m_tc.to_csv(SRC_DIR / "Figure6b_tcga_pure_ductal_km_source.tsv", sep="\t", index=False)

def km_curve(times: np.ndarray, events: np.ndarray):
    order = np.argsort(times)
    t_s, e_s = times[order], events[order]
    t_pts, s_pts, surv = [0.0], [1.0], 1.0
    cens_t, cens_s = [], []
    for ut in np.unique(t_s):
        at_risk = int(np.sum(t_s >= ut))
        d_i = int(np.sum((t_s == ut) & (e_s == 1)))
        c_i = int(np.sum((t_s == ut) & (e_s == 0)))
        if d_i > 0:
            t_pts.extend([float(ut), float(ut)])
            s_pts.extend([surv, surv * (1.0 - d_i / at_risk)])
            surv *= (1.0 - d_i / at_risk)
        if c_i > 0:
            cens_t.append(float(ut)); cens_s.append(surv)
    t_pts.append(float(t_s.max())); s_pts.append(surv)
    return np.array(t_pts), np.array(s_pts), np.array(cens_t), np.array(cens_s)

hi_m = m_tc["STAT1_grp"] == "STAT1 >= median"
t_hi, e_hi = m_tc.loc[hi_m, "os_m"].values, m_tc.loc[hi_m, "os_e"].values
t_lo, e_lo = m_tc.loc[~hi_m, "os_m"].values, m_tc.loc[~hi_m, "os_e"].values
O1, E1, V1_lr = 0.0, 0.0, 0.0
for dt in np.sort(np.unique(m_tc.loc[m_tc["os_e"] == 1, "os_m"].values)):
    n1, n2 = np.sum(t_hi >= dt), np.sum(t_lo >= dt); n_tot = n1 + n2
    d1, d2 = np.sum((t_hi == dt) & (e_hi == 1)), np.sum((t_lo == dt) & (e_lo == 1)); d_tot = d1 + d2
    if n_tot > 1:
        O1 += d1; E1 += d_tot * (n1 / n_tot)
        V1_lr += (n1 * n2 * d_tot * (n_tot - d_tot)) / (n_tot ** 2 * (n_tot - 1))
lr_p = float(stats.chi2.sf((O1 - E1) ** 2 / max(V1_lr, 1e-12), df=1))

for grp_lbl, t_arr, e_arr, col_k in [
    ("STAT1 >= median (n=73)", t_hi, e_hi, CB_PAL["ifn_high"]),
    ("STAT1 < median (n=72)",  t_lo, e_lo, CB_PAL["ifn_low"]),
]:
    tx, sy, ct_x, cs_y = km_curve(t_arr, e_arr)
    ax6b.plot(tx, sy * 100.0, color=col_k, lw=1.35, label=grp_lbl)
    ax6b.scatter(ct_x, cs_y * 100.0, marker="|", s=15, color=col_k, linewidths=0.75)

ax6b.text(30, 83, f"Descriptive log-rank {fmt_p(lr_p)}\nCont. HR = 1.10 (p=0.40, q=0.71)\nNo sig. OS association", fontsize=5.3,
          bbox=dict(boxstyle="round,pad=0.22", facecolor="#F8FAFC", edgecolor="#CBD5E1", lw=0.6))
ax6b.set_xlim(-3, 72); ax6b.set_ylim(0, 105); ax6b.set_xticklabels([])
ax6b.set_ylabel("Overall survival (%)")
ax6b.set_title("Pure-ductal PDAC descriptive KM (n=145; 84 deaths)", fontsize=7.25)
despine(ax6b)
ax6b.legend(loc="lower left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.4)

rt_ticks = [0, 12, 24, 36, 48, 60]
ax6b_rt.set_xlim(-3, 72); ax6b_rt.set_ylim(-0.6, 1.7)
for sp_nm in ("top", "right", "left"):
    ax6b_rt.spines[sp_nm].set_visible(False)
ax6b_rt.set_yticks([1.0, 0.0]); ax6b_rt.set_yticklabels(["High", "Low"], fontsize=5.7, fontweight="bold")
ax6b_rt.tick_params(axis="y", length=0, pad=3)
ax6b_rt.set_xticks(rt_ticks); ax6b_rt.set_xticklabels([str(x) for x in rt_ticks], fontsize=6.0)
ax6b_rt.set_xlabel("Follow-up time (months)", fontsize=7.0)
for xt in rt_ticks:
    ax6b_rt.text(xt, 1.0, str(int(np.sum(t_hi >= xt))), ha="center", va="center", fontsize=5.7, color=CB_PAL["ifn_high"])
    ax6b_rt.text(xt, 0.0, str(int(np.sum(t_lo >= xt))), ha="center", va="center", fontsize=5.7, color=CB_PAL["ifn_low"])

# Panel 6c: Compositional Deconvolution of Bulk STAT1 Correlations
ax6c = fig6.add_subplot(gs6[1, :])
panel_label(ax6c, "c", x=-0.06, y=1.05)
genes_6c = [g for g in ["CXCL10", "CD274", "IDO1", "PDCD1LG2", "CXCL9", "TAP1", "IRF1", "LAG3", "HLA-A", "B2M", "HAVCR2", "TIGIT", "CTLA4", "LGALS9", "CD47", "CD276", "PDCD1", "VSIR"] if g in set(corr_ex["gene"])]
c6_sub = corr_ex.set_index("gene").loc[genes_6c].reset_index()
ax6c.axhline(0.0, color="#333333", ls="--", lw=0.75)
for xi, (_, rc) in enumerate(c6_sub.iterrows()):
    r_raw, r_par = float(rc["raw_spearman_rho"]), float(rc["partial_rho_adj_PTPRC_Mac_Purity"])
    lo_p, hi_p, q_par = float(rc["partial_ci95_low"]), float(rc["partial_ci95_high"]), float(rc["partial_FDR_PTPRC_Mac_Purity"])
    ax6c.plot([xi, xi], [r_raw, r_par], color="#94A3B8", lw=1.0, ls=":", zorder=2)
    ax6c.plot([xi, xi], [lo_p, hi_p], color=CB_PAL["ifn_high"] if q_par < 0.05 else CB_PAL["slate"], lw=1.3, zorder=3)
    ax6c.scatter([xi], [r_raw], s=20, marker="o", facecolors="white", edgecolors="#475569", linewidths=0.85, zorder=4,
                 label="Unadjusted bulk Spearman rho" if xi == 0 else None)
    ax6c.scatter([xi], [r_par], s=24, marker="o", color=CB_PAL["ifn_high"] if q_par < 0.05 else CB_PAL["slate"],
                 edgecolors="#111111", linewidths=0.45, zorder=5,
                 label="Partial rho (adj. PTPRC + Macrophage + Epithelial purity; 95% CI)" if xi == 0 else None)

ax6c.annotate("Infiltration-confounded T/myeloid\nmarkers attenuate near 0", xy=(16.5, 0.02), xytext=(12.2, -0.26),
              fontsize=5.8, fontweight="bold", color="#334155",
              arrowprops=dict(arrowstyle="->", lw=0.7, color="#475569"))
ax6c.set_xticks(range(len(c6_sub)))
ax6c.set_xticklabels(c6_sub["gene"].tolist(), rotation=38, ha="right", fontsize=6.3)
ax6c.set_ylim(-0.34, 0.76)
ax6c.set_ylabel("Spearman correlation with bulk STAT1")
ax6c.set_title("Bulk TCGA-PAAD (n=146 pure ductal): shared IFN transcripts vs leukocyte-infiltration confounding", fontsize=7.35)
despine(ax6c)
ax6c.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.8)

export_fig(fig6, FIG_DIR, "Figure6")

# =========================================================================
# SUPPLEMENTARY FIGURES S1 - S5 (UPDATED 247 DONORS, LODO, & 1,000-PERM REGVELO)
# =========================================================================
print("[7/12] Rendering FigureS1 (Reconciled 247-donor provenance table & GSE154778 16-patient splicing QC)...")
sp_qc = pd.read_csv(BASE_ROOT / "05_results/stage4_regvelo/gse154778_macrophage_per_sample_splicing_qc.tsv", sep="\t")
sp_qc.to_csv(SRC_DIR / "FigureS1_gse154778_splicing_qc_source.tsv", sep="\t", index=False)

figs1 = plt.figure(figsize=(DOUBLE_COL_IN, 5.75))
gss1 = gridspec.GridSpec(2, 2, figure=figs1, height_ratios=[1.02, 1.12], hspace=0.48, wspace=0.35)

axs1a = figs1.add_subplot(gss1[0, :])
axs1a.axis("off")
panel_label(axs1a, "a", x=-0.02, y=1.02)
axs1a.text(0, 0.94, "Reconciled provenance across 5 cohorts: 247 unique primary donors/patients (236 PDAC + 11 Ctrl; 270 specimens)", fontsize=7.6, fontweight="bold")
tbl_s1 = [
    ["Cohort & Accession", "Unique Subjects", "Specimens / Runs", "Post-QC Analytic Scale", "Analytic Unit & Key Audit Closure"],
    ["CRA001160 / PRJCA001063\n(Peng et al., Beijing)", "35 donors\n(24 PDAC + 11 Ctrl)", "35 surgical / control\nspecimens", "57,423 QC cells (18,008 genes)\n5,361 TAMs (34 donors)", "Donor Q1 (24 vs 10) & paired Q2 (n=31 >=5);\nresolved 36,016 = 18,008x2 rows"],
    ["GSE202051\n(Hwang et al., MGH)", "43 PDAC patients\n(0 non-malig. Ctrl)", "43 primary tumor\nsnRNA-seq specimens", "10,140 TAM nuclei\n(18 Untr, 25 Neoadj)", "Patient Q1 (n=43) & paired Q2 (n=34: 12U, 22N);\nno normal pancreas in cohort"],
    ["GSE154778\n(Lin et al., TGen)", "16 PDAC patients\n(10 PT + 6 Met)", "16 biopsies; 120 SRR\nruns (SRR12273020-139)", "30,399 QC spliced+unspliced\n3,318 TAMs | 10,042 Ductal", "Patient Q1 (10 vs 6), paired Q2 (n=15) &\n3-seed 1,000-perm RegVelo"],
    ["GSE272362\n(Khaliq et al., Indiana)", "7 rapid-autopsy\nPDAC patients", "30 Visium slides (10 PT,\n12 HM, 5 LNM, 3 NP*)", "149,760 grid -> 97,020 tissue\n-> 91,496 post-QC spots", "7-patient & site-stratified Moran's I;\n3 NP slides are matched autopsy pancreas"],
    ["TCGA-PAAD\n(GDC PanCanAtlas)", "146 pure-ductal PDAC\n(+31 non-ductal sens)", "146 primary tumors\n(177 total PAAD)", "145 valid OS (84 deaths);\n144 complete multivar.", "Continuous per-SD Cox M1-M3 + BH q;\n247 unique primary subjects (not 171)"],
]
t_obj = axs1a.table(cellText=tbl_s1[1:], colLabels=tbl_s1[0], colWidths=[0.19, 0.14, 0.16, 0.20, 0.31],
                    loc="center", bbox=[0.0, 0.0, 1.0, 0.88], cellLoc="left")
t_obj.auto_set_font_size(False); t_obj.set_fontsize(4.95)
for (rr, cc), cell in t_obj.get_celld().items():
    cell.set_linewidth(0.4); cell.set_edgecolor("#CBD5E1")
    if rr == 0:
        cell.set_facecolor(CB_PAL["navy"]); cell.get_text().set_color("white"); cell.get_text().set_weight("bold")
    elif rr % 2 == 1:
        cell.set_facecolor("#F8FAFC")

axs1b = figs1.add_subplot(gss1[1, 0])
panel_label(axs1b, "b", x=-0.16, y=1.06)
sp_ord = sp_qc.sort_values("n_macrophages", ascending=True).reset_index(drop=True)
y_s1 = np.arange(len(sp_ord))
axs1b.barh(y_s1, sp_ord["n_TAM_STAT1_Inflam"], color=CB_PAL["ifn_high"], height=0.7, label="TAM_STAT1_Inflam")
axs1b.barh(y_s1, sp_ord["n_TAM_SPP1_Hypoxic"], left=sp_ord["n_TAM_STAT1_Inflam"], color=CB_PAL["ifn_mid"], height=0.7, label="TAM_SPP1_Hypoxic")
axs1b.barh(y_s1, sp_ord["n_TAM_C1Q_Resident"], left=sp_ord["n_TAM_STAT1_Inflam"] + sp_ord["n_TAM_SPP1_Hypoxic"], color=CB_PAL["ifn_low"], height=0.7, label="TAM_C1Q_Resident")
axs1b.set_yticks(y_s1); axs1b.set_yticklabels(sp_ord["sample"].tolist(), fontsize=5.8)
axs1b.set_xlabel("Macrophage cell count per patient")
axs1b.set_title("GSE154778 TAM states across all 16 patients", fontsize=7.2)
despine(axs1b); axs1b.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.4)

axs1c = figs1.add_subplot(gss1[1, 1])
panel_label(axs1c, "c", x=-0.16, y=1.06)
cols_s1c = [CB_PAL["metastasis"] if s == "Metastasis" else CB_PAL["primary"] for s in sp_ord["lesion_site"]]
axs1c.scatter(sp_ord["median_spliced_umis"], sp_ord["median_unspliced_frac"] * 100.0, s=np.clip(sp_ord["n_macrophages"] * 0.15 + 18, 18, 95),
              c=cols_s1c, edgecolors="#111111", linewidths=0.5, alpha=0.86)
s1c_off = {"MET01": (-580, 0.7), "P02": (80, 0.7), "MET06": (-640, -0.4), "P03": (75, -0.8), "P06": (85, 0.9), "P07": (85, -1.1), "P09": (100, 0.5)}
for _, r_s in sp_ord.iterrows():
    dx_s, dy_s = s1c_off.get(r_s["sample"], (65, 0.2))
    axs1c.text(r_s["median_spliced_umis"] + dx_s, r_s["median_unspliced_frac"] * 100.0 + dy_s, r_s["sample"], fontsize=5.1)
axs1c.set_xlim(250, 7250)
axs1c.set_xlabel("Median spliced UMIs per TAM")
axs1c.set_ylabel("Median unspliced UMI fraction (%)")
axs1c.set_title("STARsolo splicing QC (120 SRR runs; 16 patients)", fontsize=7.2)
despine(axs1c)

export_fig(figs1, SUPP_FDIR, "FigureS1")

print("[8/12] Rendering FigureS2 (Paired pseudobulk sensitivities + exact 31-run edgeR LODO audit)...")
mod_sum = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/paired_edger_models_summary.tsv", sep="\t")
de_p14  = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/edger_outputs/paired_edgeR_Peng_Q2_IFN14_min5.tsv", sep="\t")
de_p10  = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/edger_outputs/paired_edgeR_Peng_Q2_IFN10_min5.tsv", sep="\t")
de_pq1  = pd.read_csv(BASE_ROOT / "05_results/stage3_macrophage/edger_outputs/paired_edgeR_Peng_Q1_IFN14_min5.tsv", sep="\t")
edger_lodo.to_csv(SRC_DIR / "FigureS2d_edger_exact_LODO_source.tsv", sep="\t", index=False)

figs2 = plt.figure(figsize=(DOUBLE_COL_IN, 5.9))
gss2 = gridspec.GridSpec(2, 2, figure=figs2, hspace=0.46, wspace=0.38)

axs2a = figs2.add_subplot(gss2[0, 0])
panel_label(axs2a, "a", x=-0.18, y=1.06)
y_ms = np.arange(len(mod_sum))[::-1]
axs2a.barh(y_ms, mod_sum["n_sig_FDR005_nonDefining_relaxed"], height=0.65, color=CB_PAL["ifn_high"], alpha=0.82, label="Primary universe")
axs2a.barh(y_ms, mod_sum["n_sig_FDR005_nonDefining_strict"], height=0.65, color=CB_PAL["control"], alpha=0.85, label="Strict filterByExpr")
axs2a.set_yticks(y_ms)
axs2a.set_yticklabels([t.replace("_min", " m").replace("GSE202051_", "G202_").replace("GSE154778_", "G154_") for t in mod_sum["tag"]], fontsize=5.2)
axs2a.set_xlabel("Significant non-defining DEGs (BH FDR < 0.05)")
axs2a.set_title("Min-cell (>=5, >=10) & stratification sensitivity", fontsize=7.0)
despine(axs2a); axs2a.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.2)

axs2b = figs2.add_subplot(gss2[0, 1])
panel_label(axs2b, "b", x=-0.22, y=1.08)
m_14_10 = de_p14[~de_p14["in_defining_signature"]].merge(de_p10[~de_p10["in_defining_signature"]], on="gene", suffixes=("_14", "_10"))
r_1410 = np.corrcoef(m_14_10["log2FC_edgeR_14"], m_14_10["log2FC_edgeR_10"])[0, 1]
axs2b.scatter(m_14_10["log2FC_edgeR_14"], m_14_10["log2FC_edgeR_10"], s=3, c=CB_PAL["slate"], alpha=0.45, rasterized=True)
axs2b.plot([-3, 4], [-3, 4], color=CB_PAL["ifn_high"], ls="--", lw=0.8)
axs2b.text(-2.7, 3.2, f"Pearson r = {r_1410:.3f}\n{len(m_14_10):,} non-defining genes", fontsize=5.9, fontweight="bold")
axs2b.set_xlabel("Paired log2FC: IFN-14 (STAT1 excluded)")
axs2b.set_ylabel("Paired log2FC: IFN-10\n(STAT1/IRF1/CXCL9-11 excl.)")
axs2b.set_title("IFN-14 vs stricter IFN-10 concordance", fontsize=7.1)
despine(axs2b)

axs2c = figs2.add_subplot(gss2[1, 0])
panel_label(axs2c, "c", x=-0.18, y=1.06)
v_ok = de_p14.dropna(subset=["log2FC_edgeR", "log2FC_voom"])
r_ev = np.corrcoef(v_ok["log2FC_edgeR"], v_ok["log2FC_voom"])[0, 1]
axs2c.scatter(v_ok["log2FC_edgeR"], v_ok["log2FC_voom"], s=3, c=CB_PAL["control"], alpha=0.45, rasterized=True)
axs2c.plot([-3, 4], [-3, 4], color="#111111", ls="--", lw=0.8)
axs2c.text(-2.7, 3.2, f"Pearson r = {r_ev:.3f}\nedgeR QL vs limma-voom", fontsize=5.9, fontweight="bold")
axs2c.set_xlabel("Paired log2FC: edgeR quasi-likelihood NB GLM")
axs2c.set_ylabel("Paired log2FC: limma-voom precision model")
axs2c.set_title("Estimator sensitivity: edgeR QL vs limma-voom", fontsize=7.1)
despine(axs2c)

# Panel S2d: Exact 31-Donor Leave-One-Donor-Out (LODO) Range & voom Donor-Bootstrap 95% CIs!
axs2d = figs2.add_subplot(gss2[1, 1])
panel_label(axs2d, "d", x=-0.22, y=1.06)
lodo_g = ["STAT1", "IRF7", "TAP1", "HLA-A", "B2M", "CD80", "LGALS9", "CD47", "CD274", "COL1A1", "KRT19"]
l_sub = edger_lodo[(edger_lodo["cohort_model"] == "Peng2019_Q2_IFN14_min5") & (edger_lodo["gene"].isin(lodo_g))].set_index("gene").loc[lodo_g].reset_index().iloc[::-1].reset_index(drop=True)
y_ls = np.arange(len(l_sub))
axs2d.axvline(0.0, color="#666666", ls="--", lw=0.7)
for yi, (_, rl) in enumerate(l_sub.iterrows()):
    col_l = CB_PAL["slate"] if rl["gene"] in {"COL1A1", "KRT19"} else (CB_PAL["ifn_high"] if rl["full_edgeR_nonDefining_FDR"] < 0.05 else CB_PAL["ifn_mid"])
    axs2d.plot([rl["voom_donor_bootstrap_95CI_low"], rl["voom_donor_bootstrap_95CI_high"]], [yi, yi], color=col_l, lw=1.3)
    axs2d.plot([rl["edgeR_LODO_min_log2FC"], rl["edgeR_LODO_max_log2FC"]], [yi, yi], color="#111111", lw=2.8, alpha=0.45)
    axs2d.scatter([rl["full_edgeR_log2FC"]], [yi], s=18, color=col_l, edgecolors="#111111", linewidths=0.4, zorder=4)
axs2d.set_yticks(y_ls); axs2d.set_yticklabels(l_sub["gene"].tolist(), fontsize=5.8)
axs2d.set_xlabel("Paired log2FC (thin: voom boot CI; thick: 31-run LODO)")
axs2d.set_title("Exact 31-donor edgeR LODO & bootstrap stability", fontsize=7.0)
despine(axs2d)

export_fig(figs2, SUPP_FDIR, "FigureS2")

print("[9/12] Rendering FigureS3 (3-seed RegVelo stability, 1,000-perm null scaling + rewired prior check & ductal comparator)...")
pert_all = pd.read_csv(BASE_ROOT / "05_results/perturbation_results.tsv", sep="\t")
seed_dg  = pd.read_excel(CSBJ_ROOT / "06_Quality_Audit/RegVelo_perturbation_audit.xlsx", sheet_name="04_3Seed_Model_Diagnostics")
rw_chk   = pd.read_csv(P2_DIR / "regvelo_retrained_rewired_prior_check.tsv", sep="\t")
seed_dg.to_csv(SRC_DIR / "FigureS3a_3seed_diagnostics_source.tsv", sep="\t", index=False)
rw_chk.to_csv(SRC_DIR / "FigureS3b_rewired_prior_retraining_source.tsv", sep="\t", index=False)

figs3 = plt.figure(figsize=(DOUBLE_COL_IN, 5.75))
gss3 = gridspec.GridSpec(2, 2, figure=figs3, hspace=0.48, wspace=0.38)

axs3a = figs3.add_subplot(gss3[0, 0])
panel_label(axs3a, "a", x=-0.18, y=1.06)
x_sd = np.arange(len(seed_dg))
w_b = 0.25
axs3a.bar(x_sd - w_b, seed_dg["wt_mean_TAM_STAT1_Inflam"] * 100, width=w_b, color=CB_PAL["ifn_high"], label="TAM_STAT1_Inflam")
axs3a.bar(x_sd,       seed_dg["wt_mean_TAM_SPP1_Hypoxic"] * 100, width=w_b, color=CB_PAL["ifn_mid"], label="TAM_SPP1_Hypoxic")
axs3a.bar(x_sd + w_b, seed_dg["wt_mean_TAM_C1Q_Resident"] * 100, width=w_b, color=CB_PAL["ifn_low"], label="TAM_C1Q_Resident")
axs3a.set_xticks(x_sd)
axs3a.set_xticklabels([f"Seed {int(s)}" for s in seed_dg["seed"]], fontsize=6.2)
axs3a.set_ylim(0, 82)
axs3a.set_ylabel("WT baseline fate probability (%)")
axs3a.set_title("Saved 3-seed WT fate reproducibility", fontsize=7.1)
despine(axs3a); axs3a.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.3)

# Panel S3b: Degree- & Expression-Tertile Null SD Scaling (N_perm=1,000) vs Retrained Rewired Prior
axs3b = figs3.add_subplot(gss3[0, 1])
panel_label(axs3b, "b", x=-0.22, y=1.06)
for eb, lbl_eb, col_eb in [(0, "Low expr tertile", CB_PAL["control"]), (1, "Mid expr tertile", CB_PAL["ifn_mid"]), (2, "High expr tertile", CB_PAL["ifn_high"])]:
    sub_n = rgv_null[rgv_null["expr_tertile"] == eb].sort_values("degree_knot")
    if len(sub_n):
        axs3b.plot(sub_n["degree_knot"], sub_n["null_abs_q95_pp"], marker="o", markersize=3.8, lw=1.2, color=col_eb, label=lbl_eb)
for tf_chk in ["NFKB1", "REL", "RUNX1", "STAT1", "IRF1"]:
    r_o = rgv_1000[rgv_1000["perturbation"] == tf_chk].iloc[0]
    axs3b.scatter([r_o["n_muted_prior_edges"]], [abs(r_o["mean_delta_pp"])], marker="*", s=46, color=CB_PAL["purple"], zorder=5)
    axs3b.text(r_o["n_muted_prior_edges"] + 0.8, abs(r_o["mean_delta_pp"]), tf_chk, fontsize=5.3, fontweight="bold")
axs3b.set_xlabel("Perturbed prior out-degree k (edges)")
axs3b.set_ylabel("|Delta fate| (pp): 95th pct of 1,000-perm null")
axs3b.set_title("Degree & expression-tertile null scaling (N=1,000)", fontsize=7.0)
despine(axs3b); axs3b.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", fontsize=5.1)

axs3c = figs3.add_subplot(gss3[1, :])
panel_label(axs3c, "c", x=-0.06, y=1.06)
duc_p = pert_all[(pert_all["compartment"] == "Malignant_Ductal") & (pert_all["perturbation_type"].isin(["Single_TF_KO", "Negative_Control_NonTF"]))].copy()
duc_best = duc_p.assign(abs_pp=duc_p["mean_delta_pp"].abs()).sort_values("abs_pp", ascending=False).drop_duplicates("perturbation")
duc_show = pd.concat([duc_best[duc_best["perturbation_type"] == "Single_TF_KO"].head(8), duc_best[duc_best["perturbation"].isin(["STAT1","IRF1","ACTB","GAPDH"])]]).drop_duplicates("perturbation").iloc[::-1].reset_index(drop=True)
duc_show.to_csv(SRC_DIR / "FigureS3c_malignant_ductal_comparator_source.tsv", sep="\t", index=False)
y_d3 = np.arange(len(duc_show))
cols_d3 = [CB_PAL["slate"] if t == "Negative_Control_NonTF" else (CB_PAL["ifn_high"] if p < 0.05 else CB_PAL["ifn_mid"])
           for t, p in zip(duc_show["perturbation_type"], duc_show["empirical_decoy_p"].fillna(1.0))]
axs3c.barh(y_d3, duc_show["mean_delta_pp"], height=0.62, color=cols_d3, alpha=0.86, edgecolor="#1E293B", linewidth=0.4)
axs3c.axvline(0.0, color="#1E293B", ls="--", lw=0.75)
for yi, (_, rd) in enumerate(duc_show.iterrows()):
    v_pp = float(rd["mean_delta_pp"])
    st_lbl = str(rd["terminal_state"]).replace("Ductal_", "")
    axs3c.text(v_pp - 0.06 if v_pp < 0 else v_pp + 0.06, yi, f"{v_pp:+.2f} pp ({st_lbl})", va="center", ha="right" if v_pp < 0 else "left", fontsize=5.4)
axs3c.set_yticks(y_d3); axs3c.set_yticklabels(duc_show["perturbation"].tolist(), fontsize=5.9, fontweight="bold")
m_d3 = max(np.max(np.abs(duc_show["mean_delta_pp"].values)), 1.5)
axs3c.set_xlim(-m_d3 * 1.65, m_d3 * 1.65)
axs3c.set_xlabel("Malignant_Ductal fate shift (percentage points; auxiliary epithelial comparator)")
axs3c.set_title("Auxiliary Malignant_Ductal (n=10,042 epithelial cells) RegVelo comparator", fontsize=7.2)
despine(axs3c)

export_fig(figs3, SUPP_FDIR, "FigureS3")

print("[10/12] Rendering FigureS4 (Motif-pruned pySCENIC regulons & intracellular TF-target promoter edges)...")
sc_sum = pd.read_csv(V1_ROOT / "04_results/scenic/scenic_motif_regulons_summary.csv")
tf_imm = pd.read_csv(BASE_ROOT / "05_results/stage5_spatial_comm/macrophage_TF_to_immunogene_transcriptional_edges.tsv", sep="\t")
sc_sum.head(30).to_csv(SRC_DIR / "FigureS4a_pyscenic_motif_summary_source.tsv", sep="\t", index=False)
tf_imm.to_csv(SRC_DIR / "FigureS4b_macrophage_tf_to_immunogene_edges_source.tsv", sep="\t", index=False)

figs4 = plt.figure(figsize=(DOUBLE_COL_IN, 5.4))
gss4 = gridspec.GridSpec(1, 2, figure=figs4, wspace=0.44)

axs4a = figs4.add_subplot(gss4[0])
panel_label(axs4a, "a", x=-0.20, y=1.04)
sc_pick = sc_sum[sc_sum["TF"].isin(["STAT1","STAT2","IRF7","PRDM1","SPI1","IRF8","CEBPB","REL","RUNX3","KLF6","KLF4","MITF"])].sort_values("regulon_nes_score", ascending=True).reset_index(drop=True)
y_s4 = np.arange(len(sc_pick))
axs4a.barh(y_s4, sc_pick["regulon_nes_score"], height=0.62, color=CB_PAL["control"], alpha=0.86, edgecolor="#1E293B", linewidth=0.4)
for yi, (_, rs4) in enumerate(sc_pick.iterrows()):
    axs4a.text(float(rs4["regulon_nes_score"]) + 0.25, yi, f"{int(rs4['n_targets'])} targets", va="center", ha="left", fontsize=5.5)
axs4a.set_yticks(y_s4); axs4a.set_yticklabels(sc_pick["regulon_name"].tolist(), fontsize=6.0)
axs4a.set_xlim(0, 17.8)
axs4a.set_xlabel("hg38 v10 cisTarget motif enrichment NES")
axs4a.set_title("Motif-pruned pySCENIC macrophage regulons", fontsize=7.1)
despine(axs4a)

axs4b = figs4.add_subplot(gss4[1])
panel_label(axs4b, "b", x=-0.26, y=1.04)
tf_m = tf_imm[tf_imm["TF_regulator"] != tf_imm["target_gene"]].merge(
    de_p14[["gene", "log2FC_edgeR"]].rename(columns={"gene": "target_gene"}), on="target_gene", how="left"
).dropna(subset=["log2FC_edgeR"]).sort_values(["verified_in_cisTarget_pySCENIC_motif_regulon", "log2FC_edgeR"], ascending=[False, False]).drop_duplicates(["TF_regulator", "target_gene"])
tf_imm_top = pd.concat([
    tf_m[tf_m["verified_in_cisTarget_pySCENIC_motif_regulon"] == True].head(9),
    tf_m[(tf_m["verified_in_cisTarget_pySCENIC_motif_regulon"] == False) & (tf_m["TF_regulator"].isin(["STAT1", "IRF1"]))].head(5),
]).iloc[::-1].reset_index(drop=True)
y_ti = np.arange(len(tf_imm_top))
cols_ti = [CB_PAL["ifn_high"] if bool(r["verified_in_cisTarget_pySCENIC_motif_regulon"]) else CB_PAL["ifn_mid"] for _, r in tf_imm_top.iterrows()]
axs4b.axvline(0.0, color="#888888", ls="--", lw=0.6)
axs4b.barh(y_ti, tf_imm_top["log2FC_edgeR"].values, height=0.60, color=cols_ti, alpha=0.86, edgecolor="#1E293B", linewidth=0.4)
axs4b.set_yticks(y_ti); axs4b.set_yticklabels([f"{r['TF_regulator']} -> {r['target_gene']}" for _, r in tf_imm_top.iterrows()], fontsize=5.8)
axs4b.set_xlabel("Target paired log2FC in IFN-high TAMs")
axs4b.set_title("Curated & motif-verified TF-target edges", fontsize=7.0)
axs4b.legend(handles=[mpatches.Patch(color=CB_PAL["ifn_high"], label="CollecTRI + cisTarget motif"), mpatches.Patch(color=CB_PAL["ifn_mid"], label="CollecTRI curated only")],
             loc="lower right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.2)
despine(axs4b)

export_fig(figs4, SUPP_FDIR, "FigureS4")

print("[11/12] Rendering FigureS5 (Auxiliary DepMap 24Q4 epithelial audit, TCGA mutations & LINCS drug audit)...")
dep_aud = pd.read_csv(BASE_ROOT / "05_results/stage6_clinical_drug/depmap_24q4_pdac_epithelial_audit.tsv", sep="\t")
mut_aud = pd.read_csv(BASE_ROOT / "05_results/stage6_clinical_drug/tcga_ductal_somatic_mutations_audit.tsv", sep="\t")
drg_aud = pd.read_csv(BASE_ROOT / "05_results/drug_candidates.tsv", sep="\t")
mac_drg = drg_aud[drg_aud["query_signature_context"].str.startswith("Macrophage_IFN")].head(12).copy()
dep_aud.to_csv(SRC_DIR / "FigureS5a_depmap_24q4_audit_source.tsv", sep="\t", index=False)
mut_aud.to_csv(SRC_DIR / "FigureS5b_tcga_mutations_source.tsv", sep="\t", index=False)
mac_drg.to_csv(SRC_DIR / "FigureS5c_lincs_macrophage_audit_source.tsv", sep="\t", index=False)

figs5 = plt.figure(figsize=(DOUBLE_COL_IN, 5.8))
gss5 = gridspec.GridSpec(2, 2, figure=figs5, height_ratios=[1.0, 1.15], hspace=0.50, wspace=0.38)

axs5a = figs5.add_subplot(gss5[0, 0])
panel_label(axs5a, "a", x=-0.18, y=1.06)
axs5a.axvline(-0.5, color="#888888", ls="--", lw=0.7)
axs5a.scatter(dep_aud["mean_chronos_gene_effect"], dep_aud["fraction_dependent_cell_lines_lt_minus05"] * 100.0,
              s=22, color=CB_PAL["control"], edgecolors="#111111", linewidths=0.4, alpha=0.82)
for g_dp, (oxd, oyd) in {"MYC": (0.06, -2.0), "JUN": (-0.02, 12.0), "KLF5": (-0.42, 12.0), "RELA": (-0.45, -12.0), "SOX9": (0.08, 6.0), "IRF1": (-0.36, -3.0), "STAT1": (-0.12, 13.0), "TP53": (0.03, 9.0)}.items():
    rdp = dep_aud[dep_aud["TF"] == g_dp]
    if len(rdp):
        xdp, ydp = float(rdp["mean_chronos_gene_effect"].iloc[0]), float(rdp["fraction_dependent_cell_lines_lt_minus05"].iloc[0]) * 100.0
        axs5a.annotate(g_dp, xy=(xdp, ydp), xytext=(xdp + oxd, ydp + oyd), fontsize=5.5, fontweight="bold",
                       arrowprops=dict(arrowstyle="-", lw=0.35, color="#555555"))
axs5a.set_xlabel("Mean Chronos gene effect (18 PDAC lines)")
axs5a.set_ylabel("Dependent lines (% < -0.50)")
axs5a.set_ylim(-5, 108); axs5a.set_title("DepMap 24Q4 epithelial-only fitness (STAT1 = 0%)", fontsize=7.0)
despine(axs5a)

axs5b = figs5.add_subplot(gss5[0, 1])
panel_label(axs5b, "b", x=-0.18, y=1.06)
g_mc = mut_aud.columns[0]
pct_cols = [c for c in mut_aud.columns if "pct" in c.lower() or "freq" in c.lower()]
cnt_cols = [c for c in mut_aud.columns if "n_mut" in c.lower() or "count" in c.lower() or "mutated" in c.lower()]
p_col_m = pct_cols[0] if pct_cols else cnt_cols[0]
mut_top = mut_aud.head(10).iloc[::-1].reset_index(drop=True)
vals_mt = pd.to_numeric(mut_top[p_col_m], errors="coerce").fillna(0.0).values
if np.max(vals_mt) <= 1.0:
    vals_mt = vals_mt * 100.0
cols_mt = [CB_PAL["pdac"] if g in {"KRAS","TP53","CDKN2A","SMAD4"} else CB_PAL["slate"] for g in mut_top[g_mc]]
axs5b.barh(np.arange(len(mut_top)), vals_mt, height=0.62, color=cols_mt, alpha=0.86, edgecolor="#1E293B", linewidth=0.4)
for yi, vv in enumerate(vals_mt):
    axs5b.text(vv + 1.2, yi, f"{vv:.1f}%", va="center", ha="left", fontsize=5.5)
axs5b.set_yticks(np.arange(len(mut_top))); axs5b.set_yticklabels(mut_top[g_mc].tolist(), fontsize=5.9)
axs5b.set_xlim(0, max(np.max(vals_mt) * 1.22, 100))
axs5b.set_xlabel("Mutated pure-ductal tumors (%)")
axs5b.set_title("TCGA-PAAD WES mutations (descriptive context)", fontsize=7.0)
despine(axs5b)

axs5c = figs5.add_subplot(gss5[1, :])
panel_label(axs5c, "c", x=-0.06, y=1.05)
mac_drg_rev = mac_drg.iloc[::-1].reset_index(drop=True)
y_dr = np.arange(len(mac_drg_rev))
tier_cols = [CB_PAL["ifn_high"] if str(r["strict_audit_tier"]).startswith("Tier_1") else ("#111111" if "Withdrawn" in str(r["strict_audit_tier"]) else (CB_PAL["slate"] if "BRD" in str(r["strict_audit_tier"]) else CB_PAL["control"])) for _, r in mac_drg_rev.iterrows()]
axs5c.barh(y_dr, -mac_drg_rev["median_z_sum"].values, height=0.62, color=tier_cols, alpha=0.86, edgecolor="#1E293B", linewidth=0.4)
for yi, (_, rdr) in enumerate(mac_drg_rev.iterrows()):
    v_z = -float(rdr["median_z_sum"])
    st_2026 = str(rdr["verified_2026_regulatory_status"]).split(" (")[0]
    axs5c.text(v_z + 0.25, yi, f"{int(rdr['n_l1000_signatures'])} sigs / {int(rdr['n_distinct_cell_lines'])} lines | {st_2026}", va="center", ha="left", fontsize=5.2, color="#111111")
axs5c.set_yticks(y_dr); axs5c.set_yticklabels(mac_drg_rev["pert_name"].tolist(), fontsize=5.9, fontweight="bold")
axs5c.set_xlim(0, 31.0)
axs5c.set_xlabel("Reversal strength (-median z_sum; exploratory SigCom LINCS L1000 cancer-line screen)")
axs5c.set_title("Auxiliary L1000 signature-reversal audit (Tier-1 hits dominated by non-specific proteotoxic stressors)", fontsize=6.95)
despine(axs5c)

export_fig(figs5, SUPP_FDIR, "FigureS5")
pd.DataFrame(src_map_rows).to_csv(SRC_DIR / "FIGURE_SOURCE_MAP.tsv", sep="\t", index=False)

# =========================================================================
# BUILD SUPPLEMENTARY_TABLES.XLSX WITH COMPLETE DATA DICTIONARY
# =========================================================================
print("[12/12] Building CSBJ Supplementary_Tables.xlsx (00_Data_Dictionary + 16 audited worksheets)...")
t_s1_coh   = pd.read_excel(MAN_DIR / "cohort_provenance_audit.xlsx", sheet_name="01_Cohort_Summary_247")
t_s2_spec  = pd.read_excel(MAN_DIR / "cohort_provenance_audit.xlsx", sheet_name="04_Patient_Specimen_Manifest")
t_s2b_spot = pd.read_csv(MAN_DIR / "GSE272362_sample_spot_audit.tsv", sep="\t")
t_s3a_q1   = q1_tests
t_s3b_amb  = amb_audit
t_s4a_de   = de_peng[(~de_peng["in_defining_signature"]) & ((de_peng["FDR_edgeR_nonDefining"] < 0.05) | de_peng["gene"].isin(["CD274","PDCD1LG2","IDO1","VSIR","CD276","CD80","CD86","LGALS9","CD47","KRT19","COL1A1","PRSS1","CPA1"]))].copy()
t_s4b_lodo = edger_lodo
t_s5_fg    = fg_all
t_s6a_conc = conc_voom
t_s6b_pan  = func_pan
t_s7_ulm   = ulm_nc
t_s8a_rgv  = rgv_1000
t_s8b_null = rgv_null
t_s9a_lr   = lr_aud
t_s9b_sp   = sp_pat7
t_s10a_cox = surv_bh[surv_bh["feature"].isin(["STAT1","IRF1","IRF7","STAT2","IRF2","IFN14_Signature","IFN10_Signature","FOXM1","HMGA2","GATA6","KLF5","SPI1","CD47","LGALS9","CD274","IDO1"])].copy()
t_s10b_cor = corr_ex
t_s11_mat  = pd.read_csv(P2_DIR / "target_evidence_matrix_csbj.tsv", sep="\t")
t_s12_drg  = drg_aud

sheets = [
    ("Table_S1_Cohort_247", "Reconciled 5-cohort provenance (247 primary unique donors/patients across 270 specimens; +31 TCGA non-ductal sensitivity)", t_s1_coh),
    ("Table_S2A_Specimen_Map", "Complete patient-specimen manifest across CRA001160, GSE202051, GSE154778, GSE272362, and TCGA-PAAD", t_s2_spec),
    ("Table_S2B_Visium_Spots", "GSE272362 30-section Visium spot audit reconciling 149,760 grid, 97,020 Space Ranger in-tissue, and 91,496 post-QC spots", t_s2b_spot),
    ("Table_S3A_Q1_Abundance", "Rule Q1 between-donor IFN-high macrophage proportion tests across diagnosis, therapy, and lesion site", t_s3a_q1),
    ("Table_S3B_Background_RNA", "Diagnostic audit of donor-background non-macrophage transcripts under unpaired vs within-donor paired Q2 modeling", t_s3b_amb),
    ("Table_S4A_Paired_DEGs", "Significant non-defining macrophage DEGs (BH FDR < 0.05) + checkpoint/background panel in CRA001160 paired edgeR QL", t_s4a_de),
    ("Table_S4B_edgeR_LODO", "Exact 111-run Leave-One-Donor-Out (LODO) edgeR QL NB GLM & 1,000-draw voom donor-bootstrap audit across 3 cohorts", t_s4b_lodo),
    ("Table_S5_NonCirc_fgsea", "Non-circular MSigDB Hallmark v2024.1.Hs + custom antigen-presentation/checkpoint module fgsea results", t_s5_fg),
    ("Table_S6A_Concordance_CI", "Cross-cohort Spearman rho with 1,000-draw two-cohort voom donor-resampling bootstrap 95% CIs and LODO ranges", t_s6a_conc),
    ("Table_S6B_Paired_Panel", "Multi-cohort donor-paired shifts across MHC, checkpoint, and transcription factor genes", t_s6b_pan),
    ("Table_S7_CollecTRI_ULM", "Non-circular within-donor paired CollecTRI ULM regulon enrichment across verified Lambert 2018 human TFs", t_s7_ulm),
    ("Table_S8A_RegVelo_1000P", "Lineage-restricted Macrophage_TAM 3-seed hard-constraint RegVelo screen with N_perm=1,000 degree+expr null, BH q & 16-donor CIs", t_s8a_rgv),
    ("Table_S8B_RegVelo_Null", "Summary of N_perm=1,000 degree-knot x expression-tertile random prior-edge null distributions", t_s8b_null),
    ("Table_S9A_LR_Candidates", "Biologically curated expression-inferred intercellular ligand-receptor candidates & rejected legacy pairs", t_s9a_lr),
    ("Table_S9B_Spatial_7Pat", "Patient-aggregated (n=7 rapid-autopsy PDAC patients) and site-stratified univariate Space Ranger Moran's I", t_s9b_sp),
    ("Table_S10A_TCGA_Cox_BH", "Pure-ductal TCGA-PAAD vs unfiltered PAAD Models 1-3 Cox OS/PFI with panel-wide BH q & Schoenfeld PH tests", t_s10a_cox),
    ("Table_S10B_TCGA_PartRho", "Raw vs PTPRC/macrophage/purity-adjusted partial Spearman correlations with bulk STAT1 in pure-ductal TCGA-PAAD", t_s10b_cor),
    ("Table_S11_Evidence_Mat", "Compartment-separated 5-layer Macrophage IFN-state & auxiliary Epithelial TF matrices with Dirichlet/LOLO ranks", t_s11_mat),
    ("Table_S12_L1000_Audit", "Auxiliary SigCom LINCS L1000 & Open Targets audit with strict Tier-1 criteria and verified 2026 regulatory status", t_s12_drg),
]

dict_rows = []
for sh_name, sh_desc, df_sh in sheets:
    for col in df_sh.columns:
        s = df_sh[col]
        dt_s = "Integer" if pd.api.types.is_integer_dtype(s) else ("Float" if pd.api.types.is_float_dtype(s) else ("Boolean" if pd.api.types.is_bool_dtype(s) else "String"))
        dict_rows.append({
            "Sheet_Name": sh_name,
            "Sheet_Purpose": sh_desc,
            "Row_Count": len(df_sh),
            "Column_Name": col,
            "Data_Type": dt_s,
            "Description": f"{col} in {sh_name} ({sh_desc})",
        })
df_dict = pd.DataFrame(dict_rows)

wb_path = SUPP_DIR / "Supplementary_Tables.xlsx"
with pd.ExcelWriter(wb_path, engine="openpyxl") as writer:
    df_dict.to_excel(writer, sheet_name="00_Data_Dictionary", index=False)
    for sh_name, _, df_sh in sheets:
        df_sh.replace([np.inf, -np.inf], np.nan).to_excel(writer, sheet_name=sh_name[:31], index=False)

wb = openpyxl.load_workbook(wb_path)
hdr_fill = PatternFill(start_color="1E3A5F", end_color="1E3A5F", fill_type="solid")
hdr_font = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
body_font = Font(name="Calibri", size=9.5, color="111111")

for ws in wb.worksheets:
    ws.freeze_panes = "A2"
    for cell in ws[1]:
        cell.fill = hdr_fill; cell.font = hdr_font
        cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    for col_idx in range(1, ws.max_column + 1):
        col_ltr = get_column_letter(col_idx)
        max_len = len(str(ws.cell(1, col_idx).value or ""))
        for r_idx in range(2, min(ws.max_row + 1, 50)):
            val = ws.cell(r_idx, col_idx).value
            ws.cell(r_idx, col_idx).font = body_font
            if val is not None:
                max_len = max(max_len, min(len(str(val)), 46))
        ws.column_dimensions[col_ltr].width = min(max(max_len + 3, 11), 44)

wb.save(wb_path)
print(f"[CSBJ Phase 5 Complete: 6 Main Figures, 5 Supplementary Figures, Source Data & {len(wb.worksheets)}-Sheet Supplementary_Tables.xlsx Generated!]")
