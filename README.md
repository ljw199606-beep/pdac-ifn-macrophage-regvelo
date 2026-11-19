# PDAC Interferon-Associated Macrophage States — Reproducible Code & Source Data Archive

[![Journal](https://img.shields.io/badge/Journal-CSBJ%20(Science%20Partner%20Journals)-1E3A5F)]()
[![Code License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Derived%20Data%20License-CC%20BY%204.0-lightgrey.svg)](LICENSE_DATA.txt)

Official code repository accompanying the *Computational and Structural Biotechnology Journal* (`CSBJ`) Research Article:

> **Patient-reproducible interferon-associated macrophage states in pancreatic ductal adenocarcinoma: donor-aware integration of cis-regulatory networks, RNA velocity, and spatial transcriptomics**  
> Jingwei Li, Yuchan Chen, Jiamiao Cheng, Guanwei Wang, Zhongyang Song\*, and Zhiming Zhang\*  
> *Gansu University of Chinese Medicine; Affiliated Hospital of Gansu University of Traditional Chinese Medicine; Gansu Institute of Cardiovascular Diseases; Gansu Provincial Hospital of Traditional Chinese Medicine*

---

## 1. Archive Architecture & Dual-Licensing
- **GitHub Repository (`v1.0.0`)**: Custom Python and R analysis scripts, portable figure-reproduction scripts, environment specifications, cohort provenance manifests, and documentation — licensed under the **MIT License** (`LICENSE`).
- **Zenodo Companion Archive (`v1.0.0`)**: Synchronized code snapshot plus all redistributable derived outputs (`Figure_Source_Data/`, `Supplementary_Tables.xlsx`, `RegVelo_perturbation_audit.xlsx`, canonical `phase2_tables/`, and `SHA256SUMS.txt`) — original derived data licensed under **CC BY 4.0** (`LICENSE_DATA.txt`).
- **Third-Party Raw Datasets**: Excluded from redistribution and obtained directly from their official portals under source repository terms (see [DATA_REDISTRIBUTION_AUDIT.md](DATA_REDISTRIBUTION_AUDIT.md)).

---

## 2. Obtaining Official Public Inputs
1. **CRA001160 / PRJCA001063** (Peng et al., 2019; 35 donors): https://ngdc.cncb.ac.cn/gsa/browse/CRA001160
2. **GSE202051** (Hwang et al., 2022; 43 PDAC snRNA-seq patients): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051 (`GSE202051_totaldata-final-toshare.h5ad.gz`)
3. **GSE154778** (Lin et al., 2020; 16 PDAC patients; 120 SRA runs `SRR12273020`–`SRR12273139`): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778 (`PRJNA647635` / `SRP272677`; quantify via `STARsolo v2.7.11b --soloFeatures Gene Velocyto` on GRCh38 GENCODE v44)
4. **GSE272362** (Khaliq et al., 2024; 7 rapid-autopsy patients, 30 Visium sections `GSM8399443`–`GSM8399472`): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272362 (`GSE272362_RAW.tar`)
5. **TCGA-PAAD** (NCI GDC PanCanAtlas 2018; 146 pure-ductal + 31 non-ductal sensitivity): https://gdc.cancer.gov/about-data/publications/pancanatlas
6. **DepMap Public 24Q4**: https://depmap.org/portal/
7. **SigCom LINCS L1000** & **Open Targets Platform (CC0 1.0)**: https://maayanlab.cloud/sigcom-lincs/ and https://platform.opentargets.org

---

## 3. Clean-Environment Quick Reproduction (~10 s, CPU Only)
After downloading `Figure_Source_Data/` from the companion Zenodo archive:
```bash
python -m venv .venv
source .venv/bin/activate   # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r environment/requirements_figures.txt
python scripts/reproduce_core_figures_portable.py --source-dir ./Figure_Source_Data --out-dir reproduced_figures
```
Full end-to-end execution order (requires downloaded raw datasets, R `4.4.3`, and CUDA GPU for `REGVELOVI`):
1. `scripts/pipeline_upstream/00_*` – `05c_*` (QC, non-circular IFN-14/IFN-10 scoring, paired pseudobulk, CollecTRI/pySCENIC, CellChat, TCGA base models)
2. `scripts/pipeline_csbj/csbj_01_provenance_and_spatial_spot_audit.py`
3. `scripts/pipeline_csbj/csbj_02_pseudobulk_lodo_ambient_and_bootstrap_audit.py`
4. `Rscript scripts/pipeline_csbj/csbj_02b_edger_exact_lodo_and_voom_bootstrap.R`
5. `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 python scripts/pipeline_csbj/csbj_03_regvelo_1000perm_and_donor_sensitivity.py`
6. `scripts/pipeline_csbj/csbj_04_freeze_unified_results.py`
7. `scripts/pipeline_csbj/csbj_05_render_all_figures_and_tables.py`
8. `scripts/pipeline_csbj/csbj_06_build_manuscript_and_submission_docs.py` & `csbj_07_peer_review_and_automated_audit.py`

---

## 4. Figure & Supplementary Table Traceability Map
| Deliverable | Script(s) | Source Data / Output File(s) |
|---|---|---|
| **Figure 1 (a–e)** | `csbj_01`, `02_stage3`, `csbj_05` | `Figure1a`, `1c`, `1d`, `1e` TSVs |
| **Figure 2 (a–e)** | `02_stage3`, `02b_stage3`, `csbj_02`, `csbj_02b`, `csbj_05` | `Figure2a`–`2e` TSVs |
| **Figure 3 (a–d)** | `02_stage3`, `csbj_02b`, `csbj_05`, `reproduce_core_figures_portable.py` | `Figure3a`–`3d` TSVs (`Figure3b` exact dynamic `n_genes` & matched `voom`/`edgeR` intervals) |
| **Figure 4 (a–e)** | `05b_part1`, `03_stage4`, `csbj_03`, `csbj_04`, `csbj_05`, `reproduce_core_figures_portable.py` | `Figure4a`–`4de` TSVs |
| **Figure 5 (a–d)** | `04_stage5`, `csbj_01`, `csbj_05` | `Figure5a`, `5b`, `5cd` TSVs |
| **Figure 6 (a–c)** | `05a_stage6`, `csbj_02`, `csbj_05`, `reproduce_core_figures_portable.py` | `Figure6a`–`6c` TSVs |
| **Supplementary Figures S1–S5** | `csbj_01`, `csbj_02b`, `csbj_03`, `05b_part1`, `csbj_05` | `FigureS1`–`S5` TSVs |
| **Supplementary Tables S1–S12** | `csbj_01`–`csbj_05` | `derived_results/Supplementary_Tables.xlsx` & `RegVelo_perturbation_audit.xlsx` |

---

## 5. Funding & Citation
Supported by the National Natural Science Foundation of China (Grant Nos. 82560862 and 81660730), the Gansu Provincial Science and Technology Program (Grant No. 24YFFA067), the Longyuan Youth Innovation and Entrepreneurship Talent Program (Grant No. 2024QNGR53), the Project of the Gansu Provincial Administration of Traditional Chinese Medicine (Grant No. GZKG-2024-60), Gansu Province Fourth Batch of Longyuan Youth Talents Funding Support (Grant No. 001174001), and the Gansu Provincial Yinling Expert Studio Program (Document No. 4 [2025]). See [CITATION.cff](CITATION.cff) for citation metadata.
