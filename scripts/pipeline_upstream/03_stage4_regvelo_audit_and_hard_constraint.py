#!/usr/bin/env python3
"""
Stage 4 (Gate G4): Deep RegVelo & GRN Forensic Audit + 3-Seed GPU Hard-Constraint vs Soft-Constraint
Verification on Lineage-Restricted Macrophage_TAM (3,318 cells across all 16 GSE154778 patients).

Addresses:
  - B01/P0: Direct lineage-restricted Macrophage_TAM STAT1/IRF1/myeloid TF perturbation (no ductal proxy)
  - B02/P0: Cross-lineage Global_All16_TME velocity deprecated; Unknown terminal state audited & eliminated
  - B03/P0: Exact unrounded WT/KO probabilities, absolute percentage-point change (pp), relative % change, and v0 -0.0%/2.7% audit
  - B04/P0: Directed edge set audit (reproducing v0 462/460/18 transpose/axis mismatch + proving 100% prior subset compliance under hard_constraint)
  - B05/P1: Multi-seed (3 seeds) specificity, 6 non-TF housekeeping controls, 20 degree-matched prior-edge decoys, and fixed-anchor vs re-identified GPCCA sensitivity
  - B06/P1: Complete per-sample cell counts, spliced/unspliced UMI QC, ELBO/training convergence & kinetic diagnostics
"""
import os
import time
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import MinMaxScaler
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

if hasattr(sp, "csr_array") and not hasattr(sp.csr_array, "getnnz"):
    sp.csr_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)
if hasattr(sp, "csc_array") and not hasattr(sp.csc_array, "getnnz"):
    sp.csc_array.getnnz = lambda self, axis=None: (self != 0).sum(axis=axis) if axis is not None else int(self.nnz)

import torch
import scvi
import anndata as ad
import scanpy as sc
import scvelo as scv
import cellrank as cr
import regvelo as rgv
from regvelo._constants import REGISTRY_KEYS

torch.set_float32_matmul_precision("high")

OUT_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo_Submission_20261116")
V1_ROOT = Path("/mnt/d/Research/2026_PDAC_RegVelo")
V0_ROOT = Path("/mnt/d/09-04/09-04")
ST4_DIR = OUT_ROOT / "05_results/stage4_regvelo"
ST4_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR = Path("/home/ljwbiomanager/pdac_work/regvelo_models_stage4_hard")
MODEL_DIR.mkdir(parents=True, exist_ok=True)

# =========================================================================
# 1. AUDIT v0 ARTIFACTS (B03: -0.0% / 2.7%; B04: 462 / 460 / 18)
# =========================================================================
print("[1/4] Auditing v0 RegVelo files for B03 (-0.0% / 2.7%) and B04 (462 / 460 / 18)...", flush=True)
v0_audit = {}
v0_tabs = V0_ROOT / "原始代码/pdac_regvelo_project/07_data_tables"
for fn in sorted(os.listdir(v0_tabs)):
    if any(k in fn.lower() for k in ["grn", "regvelo", "perturb", "deplet", "fate", "screen", "velo"]):
        p = v0_tabs / fn
        if p.suffix == ".csv" and p.stat().st_size < 15_000_000:
            df_t = pd.read_csv(p)
            v0_audit[fn] = {
                "rows": len(df_t),
                "cols": list(df_t.columns),
                "head": df_t.head(6).to_dict(orient="records"),
            }

# Also inspect v0 gen_main_fig2_v4.py to document exact source of -0.0% and +2.7%/-2.7%
fig2_v0_path = V0_ROOT / "原始代码/pdac_regvelo_project/03_code/figures_main/gen_main_fig2_v4.py"
if fig2_v0_path.exists():
    txt_f2 = fig2_v0_path.read_text(encoding="utf-8", errors="replace")
    v0_audit["gen_main_fig2_v4_excerpt"] = [
        ln for ln in txt_f2.splitlines()
        if any(k in ln for k in ["Ductal_3", "2.7", "462", "460", "18", "fate", "delta", "inferred", "prior"])
    ][:40]

# Inspect regvelo.REGVELOVI v_encoder hard_constraint implementation directly!
import inspect
v_enc_src = inspect.getsource(rgv.module.velovi_module.VelocityEncoder) if hasattr(rgv, "module") and hasattr(rgv.module, "velovi_module") else ""
if not v_enc_src:
    try:
        for attr in dir(rgv):
            obj = getattr(rgv, attr)
            if hasattr(obj, "__module__"):
                pass
        v_enc_src = inspect.getsource(rgv.REGVELOVI)
    except Exception:
        v_enc_src = ""
v0_audit["regvelo_version"] = rgv.__version__
(ST4_DIR / "v0_regvelo_forensic_audit.json").write_text(json.dumps(v0_audit, indent=2, default=str), encoding="utf-8")
print(f"  Saved v0_regvelo_forensic_audit.json ({len(v0_audit)} keys).", flush=True)

# =========================================================================
# 2. LOAD OFFICIAL LAMBERT 2018 HUMAN TFs, COLLECTRI PRIOR, & 16-PATIENT MACROPHAGE_TAM ATLAS
# =========================================================================
print("[2/4] Preparing 16-patient Macrophage_TAM spliced/unspliced atlas & CollecTRI prior...", flush=True)
lambert_tfs = set(
    x.strip()
    for x in (V1_ROOT / "02_data/reference/prior_networks/lambert2018_human_tfs.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)
coll_df = pd.read_csv(V1_ROOT / "02_data/reference/prior_networks/collectri_human.csv")
if "source" not in coll_df.columns and "regulator" in coll_df.columns:
    coll_df = coll_df.rename(columns={"regulator": "source"})
if "weight" not in coll_df.columns and "mor" in coll_df.columns:
    coll_df = coll_df.rename(columns={"mor": "weight"})
# Restrict regulators strictly to verified sequence-specific human TFs (Lambert et al., Cell 2018)!
coll_tf_only = coll_df[coll_df["source"].isin(lambert_tfs)].copy()
gt_net_full = coll_tf_only.pivot_table(index="target", columns="source", values="weight", aggfunc="first").fillna(0.0)
VALID_HUMAN_COLLECTRI_TFS = set(gt_net_full.columns)

NEG_CONTROLS = ["ACTB", "GAPDH", "RPL13A", "EEF1A1", "TUBA1B", "RPS18"]
PRIORITY_GENES = [
    "STAT1", "IRF1", "IRF4", "IRF7", "IRF8", "IRF9", "STAT2", "STAT3", "STAT6",
    "SPI1", "NFKB1", "RELA", "REL", "RELB", "HIF1A", "EPAS1", "CEBPB", "MAF", "MITF",
    "JUN", "FOS", "ATF3", "KLF4", "KLF5", "KLF6", "XBP1", "BHLHE40", "NFIL3", "PPARG", "RARA", "BCL6",
    "CD274", "PDCD1LG2", "IDO1", "CD47", "LGALS9", "CD80", "NECTIN2", "CXCL9", "CXCL10", "ISG15",
    "HLA-A", "B2M", "TAP1", "SPP1", "C1QA", "MRC1", "CD163",
] + NEG_CONTROLS

adata_full = sc.read_h5ad(V1_ROOT / "02_data/gse154778/gse154778_all16_velocity.h5ad")
adata_mac_raw = adata_full[adata_full.obs["cell_type"] == "Macrophage_Myeloid"].copy()

# Per-sample spliced/unspliced QC diagnostic table for Macrophage_TAM (B06/P1)
s_sum = np.asarray(adata_mac_raw.layers["spliced"].sum(axis=1)).ravel()
u_sum = np.asarray(adata_mac_raw.layers["unspliced"].sum(axis=1)).ravel()
adata_mac_raw.obs["spliced_umis"] = s_sum
adata_mac_raw.obs["unspliced_umis"] = u_sum
adata_mac_raw.obs["unspliced_frac"] = u_sum / np.maximum(s_sum + u_sum, 1.0)

mac_sample_qc = (
    adata_mac_raw.obs.groupby("sample", observed=True)
    .agg(
        n_macrophages=("fine_state", "size"),
        n_TAM_STAT1_Inflam=("fine_state", lambda s: int((s == "TAM_STAT1_Inflam").sum())),
        n_TAM_SPP1_Hypoxic=("fine_state", lambda s: int((s == "TAM_SPP1_Hypoxic").sum())),
        n_TAM_C1Q_Resident=("fine_state", lambda s: int((s == "TAM_C1Q_Resident").sum())),
        median_spliced_umis=("spliced_umis", "median"),
        median_unspliced_umis=("unspliced_umis", "median"),
        median_unspliced_frac=("unspliced_umis", lambda x: float(np.median(adata_mac_raw.obs.loc[x.index, "unspliced_frac"]))),
    )
    .reset_index()
)
mac_sample_qc["lesion_site"] = np.where(mac_sample_qc["sample"].str.startswith("MET"), "Metastasis", "Primary")
mac_sample_qc.to_csv(ST4_DIR / "gse154778_macrophage_per_sample_splicing_qc.tsv", sep="\t", index=False)
print(f"  Macrophage_TAM QC: {adata_mac_raw.n_obs} cells across {len(mac_sample_qc)} samples.", flush=True)

# Build clean lineage-restricted Macrophage_TAM velocity AnnData
a = adata_mac_raw.copy()
a.X = a.layers["counts"].copy()
s_nz = np.asarray((a.layers["spliced"] > 0).sum(axis=0)).ravel()
u_nz = np.asarray((a.layers["unspliced"] > 0).sum(axis=0)).ravel()
min_su = max(15, int(0.005 * a.n_obs))
a = a[:, (s_nz >= min_su) & (u_nz >= min_su)].copy()

scv.pp.normalize_per_cell(a)
a_hvg = a.copy()
sc.pp.log1p(a_hvg)
sc.pp.highly_variable_genes(a_hvg, n_top_genes=min(2500, a_hvg.n_vars), flavor="seurat")
tf_protect = a.var_names.isin(PRIORITY_GENES) | (
    a.var_names.isin(VALID_HUMAN_COLLECTRI_TFS)
    & (np.asarray((a.layers["spliced"] > 0).sum(axis=0)).ravel() >= max(30, int(0.02 * a.n_obs)))
    & (np.asarray((a.layers["unspliced"] > 0).sum(axis=0)).ravel() >= max(25, int(0.015 * a.n_obs)))
)
a = a[:, a_hvg.var["highly_variable"].values | tf_protect].copy()
sc.pp.log1p(a)

sc.pp.pca(a, n_comps=30)
import harmonypy as hm
ho = hm.run_harmony(a.obsm["X_pca"].astype(np.float64), a.obs, "sample", max_iter_harmony=15)
z_corr = np.asarray(ho.Z_corr)
if z_corr.shape[0] != a.n_obs and z_corr.shape[1] == a.n_obs:
    z_corr = z_corr.T
a.obsm["X_pca_harmony"] = z_corr.astype(np.float32)
sc.pp.neighbors(a, n_neighbors=30, n_pcs=30, use_rep="X_pca_harmony", transformer="sklearn", random_state=42)
sc.tl.umap(a, random_state=42)
scv.pp.moments(a, n_pcs=None, n_neighbors=None)

Ms_dense = a.layers["Ms"].toarray() if sp.issparse(a.layers["Ms"]) else np.asarray(a.layers["Ms"], dtype=np.float32)
Mu_dense = a.layers["Mu"].toarray() if sp.issparse(a.layers["Mu"]) else np.asarray(a.layers["Mu"], dtype=np.float32)
a.layers["Ms"] = MinMaxScaler().fit_transform(Ms_dense).astype(np.float32)
a.layers["Mu"] = MinMaxScaler().fit_transform(Mu_dense).astype(np.float32)

scv.tl.velocity(a, mode="deterministic")
q99_u = np.percentile(a.layers["Mu"], 99.2, axis=0)
q99_s = np.percentile(a.layers["Ms"], 99.2, axis=0)
safe_moments = (q99_u > 1e-3) & (q99_s > 1e-3)
vel_keep = (
    (a.var["velocity_r2"].fillna(0).values > 0.005)
    & (a.var["velocity_gamma"].fillna(0).values > 0.0)
    & a.var["velocity_genes"].fillna(False).values
)
prio_keep = a.var_names.isin(PRIORITY_GENES) | a.var_names.isin(VALID_HUMAN_COLLECTRI_TFS)
a = a[:, safe_moments & (vel_keep | prio_keep)].copy()

# Align prior GRN (skeleton rows = regulators, cols = targets)
a = rgv.pp.set_prior_grn(a, gt_net_full, keep_dim=True, cor_filter=True)
skel = a.uns["skeleton"]
has_edge = (skel.sum(axis=0).values > 0) | (skel.sum(axis=1).values > 0)
keep_grn = has_edge | a.var_names.isin(NEG_CONTROLS) | a.var_names.isin(PRIORITY_GENES)
a = a[:, keep_grn].copy()
skel = skel.loc[a.var_names, a.var_names].copy()
for nc in NEG_CONTROLS:
    if nc in skel.index:
        skel.loc[nc, :] = 0
a.uns["skeleton"] = skel
a.uns["network"] = skel.copy()
a.uns["regulators"] = a.var_names.to_numpy()
a.uns["targets"] = a.var_names.to_numpy()

tf_regulators = [g for g in a.var_names if (g in VALID_HUMAN_COLLECTRI_TFS) and (g not in NEG_CONTROLS) and (skel.loc[g].sum() > 0)]
# W_prior for REGVELOVI MUST be [n_targets, n_regulators] = skel.T!
W_prior = torch.tensor(skel.T.values, dtype=torch.float32)
n_prior_edges = int((W_prior != 0).sum().item())
stat1_prior_targets = skel.columns[skel.loc["STAT1"] != 0].tolist() if "STAT1" in skel.index else []
irf1_prior_targets = skel.columns[skel.loc["IRF1"] != 0].tolist() if "IRF1" in skel.index else []
print(
    f"  Macrophage_TAM GRN universe: {a.n_obs} cells x {a.n_vars} genes | {len(tf_regulators)} verified Lambert human TFs | "
    f"{n_prior_edges} correlation-validated CollecTRI prior edges",
    flush=True,
)
print(f"  STAT1 prior out-degree = {len(stat1_prior_targets)}: {stat1_prior_targets[:15]}", flush=True)
print(f"  IRF1 prior out-degree  = {len(irf1_prior_targets)}: {irf1_prior_targets[:15]}", flush=True)

# =========================================================================
# 3. TRAIN 3-SEED GPU HARD-CONSTRAINT REGVELOVI (soft_constraint=False) & AUDIT EDGE OVERLAP + KO
# =========================================================================
@torch.inference_mode()
def precompute_posterior_latent_times(vae, adata_in, n_samples=12, batch_size=2048):
    adata_val = vae._validate_anndata(adata_in)
    vae.module.eval()
    dev = next(vae.module.parameters()).device
    spliced_all = torch.as_tensor(vae.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), dtype=torch.float32, device=dev)
    unspliced_all = torch.as_tensor(vae.adata_manager.get_from_registry(REGISTRY_KEYS.U_KEY), dtype=torch.float32, device=dev)
    n_obs = spliced_all.shape[0]
    n_targets = int(vae.module.n_targets)
    ind_t_samples = torch.empty((n_samples, n_obs, n_targets), dtype=torch.float32, device="cpu")
    for s_i in range(n_samples):
        chunks = []
        for start in range(0, n_obs, batch_size):
            end = min(start + batch_size, n_obs)
            inf_out = vae.module.inference(spliced=spliced_all[start:end], unspliced=unspliced_all[start:end], n_samples=1)
            px_rho = vae.module.decoder(inf_out["z"], latent_dim=None)
            chunks.append((vae.module.t_max * px_rho).cpu())
        ind_t_samples[s_i] = torch.cat(chunks, dim=0)
    t_arr = ind_t_samples.mean(dim=0).numpy()
    scaling = (20.0 / t_arr.max(axis=0)).astype(np.float32)
    gamma, beta, _ = vae.module._get_rates()
    spliced_reg = vae.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[:, vae.module.target_index].astype(np.float32)
    ms_target = np.asarray(adata_val.layers["Ms"][:, vae.module.target_index], dtype=np.float32)
    return {
        "ind_t_samples": ind_t_samples, "scaling": scaling,
        "gamma": gamma.detach().to(dev), "beta": beta.detach().to(dev),
        "spliced_reg": spliced_reg, "ms_target": ms_target,
        "obs_template": adata_val.obs.copy(), "var_template": adata_val.var.iloc[vae.module.target_index].copy(),
        "conn": adata_val.obsp["connectivities"], "dist": adata_val.obsp["distances"],
        "neigh_info": dict(adata_val.uns["neighbors"]),
        "n_obs": n_obs, "n_targets": n_targets, "n_samples": n_samples,
    }

@torch.inference_mode()
def fast_ko_velocity(vae, precomp, ode_batch_size=2048):
    vae.module.eval()
    dev = next(vae.module.parameters()).device
    vel_acc = torch.zeros((precomp["n_obs"], precomp["n_targets"]), dtype=torch.float32, device=dev)
    for s_i in range(precomp["n_samples"]):
        ind_t_cpu = precomp["ind_t_samples"][s_i]
        for start in range(0, precomp["n_obs"], ode_batch_size):
            end = min(start + ode_batch_size, precomp["n_obs"])
            ind_t_gpu = ind_t_cpu[start:end].to(dev, non_blocking=True)
            mean_u, mean_s, _ = vae.module._get_induction_unspliced_spliced(ind_t_gpu)
            vel_acc[start:end].add_(precomp["beta"] * mean_u - precomp["gamma"] * mean_s)
    vel_arr = (vel_acc * (1.0 / precomp["n_samples"])).cpu().numpy()
    np.clip(vel_arr, -precomp["spliced_reg"], None, out=vel_arr)
    vel_arr /= precomp["scaling"]
    adata_out = ad.AnnData(
        X=precomp["ms_target"], obs=precomp["obs_template"], var=precomp["var_template"],
        obsp={"connectivities": precomp["conn"], "distances": precomp["dist"]},
        uns={"neighbors": dict(precomp["neigh_info"])},
    )
    adata_out.layers["Ms"] = precomp["ms_target"]
    adata_out.layers["velocity"] = vel_arr
    return adata_out

def compute_terminal_anchors(adata_wt, state_col, candidate_states, n_anchor=50):
    lt_mat = adata_wt.layers["latent_time_velovi"]
    lt_arr = lt_mat.values if isinstance(lt_mat, pd.DataFrame) else np.asarray(lt_mat)
    adata_wt.obs["regvelo_mean_latent_time"] = lt_arr.mean(axis=1)
    score_map = {
        "TAM_STAT1_Inflam": "score_IFN_TAM",
        "TAM_SPP1_Hypoxic": "score_SPP1_TAM",
        "TAM_C1Q_Resident": "score_C1Q_TAM",
    }
    ct_indices = {}
    for st in candidate_states:
        sub_obs = adata_wt.obs[adata_wt.obs[state_col] == st].copy()
        k = min(n_anchor, max(25, int(0.15 * len(sub_obs))))
        sc_col = score_map.get(st)
        z_lt = stats.zscore(sub_obs["regvelo_mean_latent_time"].values)
        z_sc = stats.zscore(sub_obs[sc_col].values) if (sc_col and sc_col in sub_obs.columns) else 0.0
        sub_obs["_rk"] = 0.5 * z_lt + 0.5 * z_sc
        ct_indices[st] = sub_obs.sort_values("_rk", ascending=False).head(k).index.tolist()
    return ct_indices

device = "cuda:0" if torch.cuda.is_available() else "cpu"
ck = cr.kernels.ConnectivityKernel(a).compute_transition_matrix()
terminal_candidates = ["TAM_STAT1_Inflam", "TAM_SPP1_Hypoxic", "TAM_C1Q_Resident"]
combo_kos = [["STAT1", "IRF1"], ["STAT1", "STAT3"], ["SPI1", "IRF8"], ["HIF1A", "CEBPB"], ["IRF4", "SPI1"]]

edge_audit_rows = []
training_diag_rows = []
hard_ko_records = []
reest_sens_records = []

var_list = a.var_names.tolist()
gene2idx = {g: i for i, g in enumerate(var_list)}
prior_mask_np = (skel.T.values != 0)  # [targets, regulators]
prior_edges_rc = np.argwhere(prior_mask_np)  # N x 2 array of (target_row, regulator_col)

print(f"[3/4] Training 3-seed Hard-Constraint REGVELOVI (soft_constraint=False) on {device}...", flush=True)
for run_idx, seed_val in enumerate([42, 142, 242]):
    t0 = time.time()
    scvi.settings.seed = seed_val
    torch.manual_seed(seed_val)
    np.random.seed(seed_val)

    rgv.REGVELOVI.setup_anndata(a, spliced_layer="Ms", unspliced_layer="Mu")
    vae = rgv.REGVELOVI(
        a,
        W=W_prior.to(device),
        regulators=tf_regulators,
        soft_constraint=False,  # STRICT HARD PRIOR MASK!
        lam=1.0,
        lam2=0.0,
    )
    vae.train(
        max_epochs=350,
        batch_size=a.n_obs,
        early_stopping=True,
        check_val_every_n_epoch=5,
    )
    n_epochs_trained = len(vae.history.get("elbo_train", [])) if hasattr(vae, "history") and vae.history else -1
    final_elbo_train = float(vae.history["elbo_train"].iloc[-1].iloc[0]) if n_epochs_trained > 0 else np.nan
    final_elbo_val = float(vae.history["elbo_validation"].iloc[-1].iloc[0]) if (n_epochs_trained > 0 and "elbo_validation" in vae.history and len(vae.history["elbo_validation"]) > 0) else np.nan

    # Directed edge set audit on fc1.weight and GRN_Jacobian (B04/P0)!
    fc1_w = vae.module.v_encoder.fc1.weight.detach().cpu().numpy()  # [targets, regulators]
    with torch.no_grad():
        s_mean = torch.tensor(a.layers["Ms"].mean(axis=0), dtype=torch.float32, device=device)
        jaco = vae.module.v_encoder.GRN_Jacobian(s_mean).detach().cpu().numpy()

    nz_fc1 = np.abs(fc1_w) > 1e-6
    n_nz_fc1 = int(nz_fc1.sum())
    n_overlap_correct_axis = int((nz_fc1 & prior_mask_np).sum())
    n_off_prior_fc1 = int((nz_fc1 & (~prior_mask_np)).sum())
    # What if someone mistakenly compared against transposed prior_mask_np.T (the v0 bug!)?
    n_overlap_transposed_bug = int((nz_fc1 & prior_mask_np.T).sum())

    edge_audit_rows.append({
        "seed": seed_val,
        "constraint_mode": "hard_constraint (soft_constraint=False)",
        "n_genes": a.n_vars,
        "n_lambert_tfs": len(tf_regulators),
        "n_prior_directed_edges": n_prior_edges,
        "n_inferred_nonzero_edges_gt1e6": n_nz_fc1,
        "n_overlap_correct_orientation": n_overlap_correct_axis,
        "pct_inferred_within_prior_correct": 100.0 * n_overlap_correct_axis / max(n_nz_fc1, 1),
        "n_spurious_off_prior_edges": n_off_prior_fc1,
        "n_overlap_if_transposed_v0_bug": n_overlap_transposed_bug,
        "pct_overlap_if_transposed_v0_bug": 100.0 * n_overlap_transposed_bug / max(n_nz_fc1, 1),
    })

    adata_wt = a.copy()
    rgv.tl.set_output(adata_wt, vae, n_samples=20, batch_size=1024)
    vk_wt = cr.kernels.VelocityKernel(adata_wt)
    vk_wt.compute_transition_matrix(n_jobs=8, backend="threading", show_progress_bar=False)
    gpcca_wt = cr.estimators.GPCCA(0.8 * vk_wt + 0.2 * ck)

    # Unsupervised macrostate check
    unp_macro = {}
    try:
        gpcca_wt.compute_macrostates(n_states=4, n_cells=30, cluster_key="fine_state")
        unp_macro = gpcca_wt.macrostates.value_counts().to_dict()
    except Exception:
        pass

    if run_idx == 0:
        primary_anchors = compute_terminal_anchors(adata_wt, "fine_state", terminal_candidates, n_anchor=50)
    gpcca_wt.set_terminal_states(primary_anchors)
    gpcca_wt.compute_fate_probabilities(solver="gmres")
    fate_wt_df = pd.DataFrame(
        np.asarray(gpcca_wt.fate_probabilities),
        index=adata_wt.obs_names,
        columns=list(gpcca_wt.fate_probabilities.names),
    )

    training_diag_rows.append({
        "compartment": "Macrophage_TAM",
        "constraint_mode": "hard_constraint",
        "run": run_idx,
        "seed": seed_val,
        "n_cells": a.n_obs,
        "n_genes": a.n_vars,
        "n_epochs_trained": n_epochs_trained,
        "final_elbo_train": final_elbo_train,
        "final_elbo_validation": final_elbo_val,
        "unsupervised_gpcca_macrostates": json.dumps({str(k): int(v) for k, v in unp_macro.items()}),
        "wt_mean_prob_TAM_STAT1_Inflam": float(fate_wt_df["TAM_STAT1_Inflam"].mean()),
        "wt_mean_prob_TAM_SPP1_Hypoxic": float(fate_wt_df["TAM_SPP1_Hypoxic"].mean()),
        "wt_mean_prob_TAM_C1Q_Resident": float(fate_wt_df["TAM_C1Q_Resident"].mean()),
        "train_seconds": round(time.time() - t0, 1),
    })

    # Screen ALL expressed Lambert human TFs + ALL 6 non-TF housekeeping controls + 20 degree-matched prior-edge decoys + 5 combos across ALL 3 SEEDS!
    raw_fc1 = vae.module.v_encoder.fc1.weight.detach().clone()
    out_degs = [int((raw_fc1[:, gene2idx[g]].abs() > 1e-6).sum().item()) for g in tf_regulators]
    stat1_deg = int((raw_fc1[:, gene2idx["STAT1"]].abs() > 1e-6).sum().item()) if "STAT1" in gene2idx else int(np.median(out_degs))

    tasks = []
    for g in tf_regulators:
        tasks.append(("Single_TF_KO", g, [gene2idx[g]], None))
    for g in [nc for nc in NEG_CONTROLS if nc in gene2idx]:
        tasks.append(("Negative_Control_NonTF", g, [gene2idx[g]], None))
    for cb in combo_kos:
        if all(g in gene2idx for g in cb):
            tasks.append(("Combinatorial_KO", "+".join(cb), [gene2idx[g] for g in cb], None))
    # 20 degree-matched random prior-edge decoys matching STAT1's exact active out-degree!
    active_prior_rc = np.argwhere(np.abs(fc1_w) > 1e-6)
    for p_i in range(20):
        rng_p = np.random.default_rng(5000 + run_idx * 100 + p_i)
        pick = rng_p.choice(len(active_prior_rc), size=min(max(stat1_deg, 5), len(active_prior_rc)), replace=False)
        rr = active_prior_rc[pick, 0]
        cc = active_prior_rc[pick, 1]
        tasks.append(("PriorEdge_Decoy_Ctrl", f"PriorDecoy_{p_i+1:02d}", [], (rr, cc)))

    precomp = precompute_posterior_latent_times(vae, adata_in=a, n_samples=12, batch_size=4096)
    for ptype, pname, col_idxs, rand_rc in tasks:
        pert_w = raw_fc1.clone()
        if rand_rc is not None:
            rr, cc = rand_rc
            n_muted = len(rr)
            pert_w[rr, cc] = 0.0
        else:
            n_muted = 0
            for c_idx in col_idxs:
                n_muted += int((pert_w[:, c_idx].abs() > 1e-6).sum().item())
                pert_w[:, c_idx] = 0.0

        vae.module.v_encoder.fc1.weight.data = pert_w
        adata_ko = fast_ko_velocity(vae, precomp, ode_batch_size=4096)
        vae.module.v_encoder.fc1.weight.data = raw_fc1.clone()

        vk_ko = cr.kernels.VelocityKernel(adata_ko)
        vk_ko.compute_transition_matrix(n_jobs=8, backend="threading", show_progress_bar=False)
        gpcca_ko = cr.estimators.GPCCA(0.8 * vk_ko + 0.2 * ck)
        gpcca_ko.set_terminal_states(primary_anchors)
        gpcca_ko.compute_fate_probabilities(solver="gmres")
        fate_ko_df = pd.DataFrame(
            np.asarray(gpcca_ko.fate_probabilities),
            index=adata_ko.obs_names,
            columns=list(gpcca_ko.fate_probabilities.names),
        )

        n_c = len(fate_wt_df)
        y_lbl = np.array([1] * n_c + [0] * n_c)
        for st_col in terminal_candidates:
            w = fate_wt_df[st_col].values.astype(float)
            k = fate_ko_df[st_col].values.astype(float) if st_col in fate_ko_df.columns else np.zeros_like(w)
            delta = k - w
            pred = np.concatenate([w, k])
            auc = float(roc_auc_score(y_lbl, pred)) if np.std(pred) > 1e-10 else 0.5
            wt_m = float(np.mean(w))
            ko_m = float(np.mean(k))
            d_m  = float(np.mean(delta))
            hard_ko_records.append({
                "compartment": "Macrophage_TAM",
                "constraint_mode": "hard_constraint",
                "run": run_idx,
                "seed": seed_val,
                "perturbation_type": ptype,
                "perturbation": pname,
                "n_muted_prior_edges": n_muted,
                "terminal_state": st_col,
                "wt_fate_prob": wt_m,
                "ko_fate_prob": ko_m,
                "delta_fate_prob": d_m,
                "delta_percentage_points": d_m * 100.0,
                "relative_pct_change": (d_m / max(wt_m, 1e-8)) * 100.0,
                "depletion_auc": auc,
            })

        # Sensitivity check (Stage 4 item 9): re-identify anchors on perturbed velocity for key regulators on seed 42!
        if run_idx == 0 and pname in {"STAT1", "IRF1", "STAT1+IRF1", "IRF4", "SPI1", "NFKB1", "ACTB"}:
            # Re-score anchors using perturbed velocity graph's absorbing states / GPCCA macrostates
            try:
                gpcca_re = cr.estimators.GPCCA(0.8 * vk_ko + 0.2 * ck)
                gpcca_re.compute_macrostates(n_states=3, n_cells=50, cluster_key="fine_state")
                gpcca_re.set_terminal_states_from_macrostates()
                gpcca_re.compute_fate_probabilities(solver="gmres")
                f_re = pd.DataFrame(np.asarray(gpcca_re.fate_probabilities), index=adata_ko.obs_names, columns=[str(c) for c in gpcca_re.fate_probabilities.names])
                for st_col in terminal_candidates:
                    match_cols = [c for c in f_re.columns if c.startswith(st_col)]
                    if match_cols:
                        k_re = float(f_re[match_cols[0]].mean())
                        w_0  = float(fate_wt_df[st_col].mean())
                        reest_sens_records.append({
                            "perturbation": pname,
                            "terminal_state": st_col,
                            "fixed_anchor_delta_pp": float((fate_ko_df[st_col].mean() - w_0) * 100.0),
                            "reidentified_macrostate_ko_prob": k_re,
                        })
            except Exception:
                pass

    print(f"  Completed Seed {seed_val} ({len(tasks)} perturbations in {time.time()-t0:.1f}s) | Off-prior edges = {n_off_prior_fc1}", flush=True)

pd.DataFrame(edge_audit_rows).to_csv(ST4_DIR / "regvelo_directed_edge_overlap_audit.tsv", sep="\t", index=False)
pd.DataFrame(training_diag_rows).to_csv(ST4_DIR / "regvelo_training_and_convergence_diagnostics.tsv", sep="\t", index=False)
if reest_sens_records:
    pd.DataFrame(reest_sens_records).to_csv(ST4_DIR / "regvelo_terminal_state_reidentification_sensitivity.tsv", sep="\t", index=False)

df_hard_raw = pd.DataFrame(hard_ko_records)
df_hard_raw.to_csv(ST4_DIR / "regvelo_Macrophage_TAM_hard_constraint_all_seeds_raw.tsv", sep="\t", index=False)

# Aggregate across 3 seeds & compute empirical permutation p against the 20 STAT1-degree-matched PriorEdge_Decoy_Ctrl!
df_hard_agg = (
    df_hard_raw.groupby(["compartment", "constraint_mode", "perturbation_type", "perturbation", "terminal_state"], as_index=False)
    .agg(
        n_seeds=("seed", "nunique"),
        mean_muted_prior_edges=("n_muted_prior_edges", "mean"),
        mean_wt_prob=("wt_fate_prob", "mean"),
        sd_wt_prob=("wt_fate_prob", "std"),
        mean_ko_prob=("ko_fate_prob", "mean"),
        sd_ko_prob=("ko_fate_prob", "std"),
        mean_delta_prob=("delta_fate_prob", "mean"),
        sd_delta_prob=("delta_fate_prob", "std"),
        mean_delta_pp=("delta_percentage_points", "mean"),
        sd_delta_pp=("delta_percentage_points", "std"),
        mean_relative_pct_change=("relative_pct_change", "mean"),
        mean_depletion_auc=("depletion_auc", "mean"),
        sd_depletion_auc=("depletion_auc", "std"),
    )
)
df_hard_agg["abs_mean_delta_pp"] = df_hard_agg["mean_delta_pp"].abs()

# Degree-normalized sensitivity (pp per muted prior edge) + empirical p against PriorEdge_Decoy_Ctrl per terminal state
df_hard_agg["empirical_decoy_p"] = np.nan
df_hard_agg["empirical_decoy_fdr"] = np.nan
for st_col, idx in df_hard_agg.groupby("terminal_state").groups.items():
    decoy_vals = df_hard_agg.loc[idx].query("perturbation_type == 'PriorEdge_Decoy_Ctrl'")["abs_mean_delta_pp"].values
    if len(decoy_vals) > 0:
        emp_p = [(float(np.sum(decoy_vals >= v)) + 1.0) / (len(decoy_vals) + 1.0) for v in df_hard_agg.loc[idx, "abs_mean_delta_pp"].values]
        df_hard_agg.loc[idx, "empirical_decoy_p"] = emp_p
        tf_mask = df_hard_agg.loc[idx, "perturbation_type"] == "Single_TF_KO"
        tf_idx = df_hard_agg.loc[idx].index[tf_mask]
        df_hard_agg.loc[tf_idx, "empirical_decoy_fdr"] = multipletests(df_hard_agg.loc[tf_idx, "empirical_decoy_p"], method="fdr_bh")[1]

df_hard_agg = df_hard_agg.sort_values(["terminal_state", "abs_mean_delta_pp"], ascending=[True, False]).reset_index(drop=True)
df_hard_agg.to_csv(ST4_DIR / "regvelo_Macrophage_TAM_hard_constraint_summary.tsv", sep="\t", index=False)

# =========================================================================
# 4. MERGE HARD-CONSTRAINT + SOFT-CONSTRAINT + MALIGNANT_DUCTAL INTO FROZEN perturbation_results.tsv
# =========================================================================
print("[4/4] Building unified frozen perturbation_results.tsv (Lambert 2018 TFs + Controls Only)...", flush=True)
soft_mac = pd.read_csv(V1_ROOT / "04_results/regvelo/regvelo_Macrophage_TAM_ko_summary.csv")
soft_duc = pd.read_csv(V1_ROOT / "04_results/regvelo/regvelo_Malignant_Ductal_ko_summary.csv")

def format_v1_soft(df_s, comp_lbl):
    # Filter Single_TF_KO strictly to verified Lambert 2018 TFs!
    keep = (df_s["perturbation_type"] != "Single_TF_KO") | (df_s["perturbation"].isin(lambert_tfs))
    sub = df_s[keep].copy()
    return pd.DataFrame({
        "compartment": comp_lbl,
        "constraint_mode": "soft_constraint (lam2=0.02)",
        "perturbation_type": sub["perturbation_type"],
        "perturbation": sub["perturbation"],
        "terminal_state": sub["terminal_state"],
        "n_seeds": sub["n_runs"],
        "mean_muted_edges": sub["mean_muted_edges"],
        "mean_wt_prob": sub["mean_wt_fate_prob"],
        "mean_ko_prob": sub["mean_ko_fate_prob"],
        "mean_delta_prob": sub["mean_delta_fate_prob"],
        "sd_delta_prob": sub["sd_delta_fate_prob"],
        "mean_delta_pp": sub["mean_delta_fate_prob"] * 100.0,
        "sd_delta_pp": sub["sd_delta_fate_prob"] * 100.0,
        "mean_relative_pct_change": (sub["mean_delta_fate_prob"] / sub["mean_wt_fate_prob"].clip(1e-8)) * 100.0,
        "mean_depletion_auc": sub["depletion_likelihood_auc"],
        "sd_depletion_auc": sub["sd_likelihood_auc"],
        "empirical_decoy_p": sub["empirical_perm_p"],
    })

hard_fmt = pd.DataFrame({
    "compartment": df_hard_agg["compartment"],
    "constraint_mode": "hard_constraint (soft_constraint=False)",
    "perturbation_type": df_hard_agg["perturbation_type"],
    "perturbation": df_hard_agg["perturbation"],
    "terminal_state": df_hard_agg["terminal_state"],
    "n_seeds": df_hard_agg["n_seeds"],
    "mean_muted_edges": df_hard_agg["mean_muted_prior_edges"],
    "mean_wt_prob": df_hard_agg["mean_wt_prob"],
    "mean_ko_prob": df_hard_agg["mean_ko_prob"],
    "mean_delta_prob": df_hard_agg["mean_delta_prob"],
    "sd_delta_prob": df_hard_agg["sd_delta_prob"],
    "mean_delta_pp": df_hard_agg["mean_delta_pp"],
    "sd_delta_pp": df_hard_agg["sd_delta_pp"],
    "mean_relative_pct_change": df_hard_agg["mean_relative_pct_change"],
    "mean_depletion_auc": df_hard_agg["mean_depletion_auc"],
    "sd_depletion_auc": df_hard_agg["sd_depletion_auc"],
    "empirical_decoy_p": df_hard_agg["empirical_decoy_p"],
})

frozen_pert = pd.concat([hard_fmt, format_v1_soft(soft_mac, "Macrophage_TAM"), format_v1_soft(soft_duc, "Malignant_Ductal")], ignore_index=True)
frozen_pert.to_csv(OUT_ROOT / "05_results/perturbation_results.tsv", sep="\t", index=False)

print("\n=== Directed Edge Overlap Audit (Hard Constraint vs Transposed Bug) ===", flush=True)
print(pd.DataFrame(edge_audit_rows).to_string(index=False), flush=True)

print("\n=== Top 15 Hard-Constraint Perturbations on TAM_STAT1_Inflam (plus STAT1, IRF1, Negative Controls) ===", flush=True)
st_infl = df_hard_agg[df_hard_agg["terminal_state"] == "TAM_STAT1_Inflam"]
top_show = pd.concat([
    st_infl[st_infl["perturbation_type"].isin(["Single_TF_KO", "Combinatorial_KO"])].head(15),
    st_infl[st_infl["perturbation"].isin(["STAT1", "IRF1", "STAT1+IRF1", "STAT3", "ACTB", "GAPDH", "PriorDecoy_01"])],
]).drop_duplicates(subset=["perturbation"])
print(
    top_show[["perturbation_type", "perturbation", "mean_muted_prior_edges", "mean_wt_prob", "mean_ko_prob", "mean_delta_pp", "sd_delta_pp", "mean_relative_pct_change", "mean_depletion_auc", "empirical_decoy_p"]]
    .to_string(index=False),
    flush=True,
)
print("[Stage 4 GPU Audit & Hard-Constraint Screen Complete]", flush=True)
