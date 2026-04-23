# DS596 Final Project: CellMentor Projection Diagnostics Extension

**Team:** Katherine Beaty, Alivia Pavel, Jessica Cannon, Mia Vargas
**Course:** DS596 Learning from Large-Scale Biological Data, Spring 2026
**Paper:** Hevdeli O, Petrenko E, Aran D. *CellMentor: cell-type aware dimensionality reduction for single-cell RNA-sequencing data.* Nature Communications 17, 396 (2026). https://doi.org/10.1038/s41467-025-67088-7

---

## Overview

This repository contains:

1. **A replication** of the CellMentor Baron → Muraro pancreas projection
   (the one non-deprecated dataset available in the current package), including
   a publication-style UMAP and ARI/NMI scores at the paper's default
   clustering resolution for direct comparison with the paper's Table 1.

2. **A projection-phase extension** that adds three per-cell diagnostics
   currently absent from CellMentor's output:
   - Normalized reconstruction error (`‖Wh − x‖² / ‖x‖²`)
   - Centroid-based confidence scores and class margin
   - A novel cell type flag calibrated against the reference error distribution

   These diagnostics address a limitation the authors explicitly acknowledge:
   cells belonging to populations absent from the reference are
   projected onto the nearest existing class with no uncertainty signal.

---

## Repository structure

```
DS596-Final-Project/
├── README.md                          ← this file
├── extension/
│   ├── diagnostics.R                  ← core extension functions
│   ├── plotting.R                     ← visualization helpers
│   ├── run_extension_demo.R           ← full end-to-end demo (Baron → Muraro)
│   ├── submit_extension.sh            ← SCC batch submission script
│   ├── quick_test.sh                  ← fast subset test (~5 min)
│   └── ADDING_TO_YOUR_PIPELINE.md    ← guide for applying to other datasets
└── replication/
    ├── replicate_figures.R            ← ARI/NMI + publication-style UMAP
    └── submit_replication.sh          ← SCC batch script
```

---

## Replication instructions

> **These instructions assume access to the BU Shared Computing Cluster (SCC).**
> All required R packages are pre-installed in the project library.
> No external data downloads are required - both datasets are bundled inside
> the CellMentor R package.

### Step 1 - Clone the repository into a working directory:

```bash
git clone https://github.com/Kbeaty26/DS596-Final-Project.git
cd DS596-Final-Project
```

### Step 2 - Install R dependencies (one time, ~20 minutes)
 
R packages need to be installed into your own personal library. Run the setup
script, which installs everything automatically and skips packages that are
already present:
 
```bash
bash setup.sh
```
 
You should see `All packages installed successfully.` at the end. If any
package shows `MISSING`, re-run `setup.sh` - intermittent network issues on
the SCC sometimes cause single packages to fail.
 
### Step 3 - Configure scripts for your username
 
All shell scripts have paths hardcoded to the submitting user's SCC directory.
Run this one-liner to replace `aliviap` with your own username throughout:
 
```bash
USERNAME=$(whoami)
sed -i "s/aliviap/${USERNAME}/g" \
  extension/submit_extension.sh \
  extension/quick_test.sh \
  replication/submit_replication.sh
```
 
Also update the R library path in the two R scripts to point to your library:
 
```bash
USERNAME=$(whoami)
sed -i "s|R_libs_4.5.2|R_libs_4.5.2|g" \
  extension/run_extension_demo.R \
  replication/replicate_figures.R
 
sed -i "s|/projectnb/ds596/students/aliviap/|/projectnb/ds596/students/${USERNAME}/|g" \
  extension/run_extension_demo.R \
  replication/replicate_figures.R
```
 
### Step 4 - Create output directories
 
```bash
USERNAME=$(whoami)
mkdir -p /projectnb/ds596/students/${USERNAME}/extension/logs
mkdir -p /projectnb/ds596/students/${USERNAME}/extension/figures
mkdir -p /projectnb/ds596/students/${USERNAME}/extension/figures/replication
```
 
### Step 5 - Run the extension demo
 
This loads the Baron (reference) and Muraro (query) pancreas datasets from the
CellMentor package, fits the model, projects the query, and runs the full
projection diagnostics extension.
 
```bash
qsub extension/submit_extension.sh
```
 
Monitor progress:
```bash
qstat -u $(whoami)
tail -f /projectnb/ds596/students/$(whoami)/extension/logs/extension.log
```
 
Outputs written to `/projectnb/ds596/students/<your_username>/extension/figures/`:
 
| File | Description |
|---|---|
| `baron_muraro_diagnostics.csv` | Per-cell reconstruction error, confidence, novelty flag |
| `01_error_distributions.png` | Reference vs. query error density with novelty threshold |
| `02_confidence_distribution.png` | Confidence histogram colored by novelty flag |
| `03_umap_all_panels.png` | 5-panel UMAP: true label, predicted, confidence, error, novelty |
| `baron_muraro_results.rds` | Full results bundle used by the replication script |
| `session_info.txt` | Full R session info for reproducibility |
 
### Step 6 - Run the replication figures
 
This loads the saved results from Step 5 (no re-training) and produces
publication-adjacent figures for comparison with the paper.
 
```bash
qsub replication/submit_replication.sh
```
 
Outputs written to `.../extension/figures/replication/`:
 
| File | Description |
|---|---|
| `figA_umap_projection.png` | Publication-style UMAP, mirrors Fig. A (p.7) of paper |
| `figB_ari_nmi_resolution.png` | ARI/NMI vs. clustering resolution sweep |
| `figB_resolution_sweep.csv` | Numeric scores at each resolution |
| `figC_ari_nmi_table.csv` | ARI/NMI at default resolution — compare with paper Table 1 |
 
### Optional: fast sanity check (~5 min, no job submission needed)
 
Runs on a 30-cells-per-type subsample to verify the pipeline end-to-end before
committing cluster time:
 
```bash
bash extension/quick_test.sh
```
 
---
 
## Pre-generated figures
 
The `figures/` directory in this repository contains output figures from our
full run so you can inspect results without re-running everything:
 
- `01_error_distributions.png`
- `02_confidence_distribution.png`
- `03_umap_all_panels.png`