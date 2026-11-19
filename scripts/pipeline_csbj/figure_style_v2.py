#!/usr/bin/env python3
"""
Unified Springer Nature / BMC Genomics Publication Figure Style Module (v2, 2026-11-16).
- Enforces double-column journal print width (175 mm = 6.89 in), colorblind-safe Okabe-Ito/Tol palettes,
  editable TrueType PDF/SVG vector fonts, 600-DPI PNG + 300-DPI TIFF exports, and zero p=0 / nan formatting.
"""
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

DOUBLE_COL_IN = 6.89   # 175 mm standard Springer Nature / BMC Genomics double-column width
SINGLE_COL_IN = 3.35   # 85 mm single-column width

# Okabe-Ito & Paul Tol colorblind-safe palettes (never red-green opposition!)
PALETTE = {
    "tumor": "#C44E52",         # muted brick-vermilion
    "normal": "#2B6CB0",        # deep cobalt blue
    "ifn_high": "#B2182B",      # deep crimson
    "ifn_mid": "#E69F00",       # warm amber
    "ifn_low": "#2166AC",       # sapphire blue
    "untreated": "#0072B2",     # ocean blue
    "neoadjuvant": "#D55E00",   # burnt orange
    "primary": "#1F4E78",       # navy
    "metastasis": "#C55A11",    # copper
    "hm": "#D95F02",            # hepatic met orange
    "lnm": "#7570B3",           # lymph node purple
    "np": "#1B9E77",            # normal pancreas teal
    "neutral": "#7F8C8D",       # slate gray
    "accent_purple": "#6A3D9A",
    "accent_teal": "#009E73",
    "grid": "#E5E8EC",
}

CELLTYPE_COLORS = {
    "Ductal cell type 2": "#B2182B",
    "Ductal cell type 1": "#4393C3",
    "Acinar cell": "#009E73",
    "Endocrine cell": "#CC79A7",
    "Macrophage cell": "#D55E00",
    "T cell": "#6A3D9A",
    "B cell": "#56B4E9",
    "Fibroblast cell": "#E69F00",
    "Stellate cell": "#8C510A",
    "Endothelial cell": "#35978F",
}

def apply_pub_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7.8,
        "axes.titlesize": 8.4,
        "axes.titleweight": "bold",
        "axes.labelsize": 7.8,
        "axes.labelweight": "normal",
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.8,
        "legend.title_fontsize": 7.2,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.8,
        "ytick.major.size": 2.8,
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
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#111111",
    )

def fmt_p(p: float, prefix: str = "p") -> str:
    """Format p/FDR strictly without p=0 or nan."""
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

def save_pub_figure(fig: plt.Figure, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_p = out_dir / f"{stem}.pdf"
    svg_p = out_dir / f"{stem}.svg"
    png_p = out_dir / f"{stem}.png"
    tif_p = out_dir / f"{stem}.tiff"
    fig.savefig(pdf_p, format="pdf")
    fig.savefig(svg_p, format="svg")
    fig.savefig(png_p, format="png", dpi=600)
    try:
        fig.savefig(tif_p, format="tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    except Exception:
        fig.savefig(tif_p, format="tiff", dpi=300)
    plt.close(fig)
    return {"pdf": str(pdf_p), "svg": str(svg_p), "png": str(png_p), "tiff": str(tif_p)}
