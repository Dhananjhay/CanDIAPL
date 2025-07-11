# CanDIAPL

**CanDIAPL** is a Snakemake-based workflow that queries the Gaia astronomical catalog using ADQL and generates an interactive HTML viewer to inspect the results.

This pipeline is ideal for quick exploration and visualization of nearby stars based on a central coordinate in RA/Dec, using the Gaia DR3 database.


## 📦 Dependencies

The pipeline assumes the following Python packages are available (can be installed via Conda or pip):

```bash
conda create -n candiapl-env snakemake
conda activate candiapl-env
```

## 📁 Project Structure

```
CanDIAPL/
├── README.md ← This file
└── workflow/
    ├── Snakefile ← Snakemake workflow definition
    ├── scripts/ ← Python scripts for querying and visualization
    │ ├── fetch_data.py ← Queries Gaia TAP service using ADQL
    │ └── view.py ← Converts query results into an HTML table
    ├── envs/ ← Conda environment YAMLs
    └── profiles/ ← Snakemake profiles for config
```


## ⚙️ Workflow Overview

The pipeline consists of the following rules:

1. **`fetch_data`**  
   - Executes an ADQL query to the Gaia TAP service using `pyvo`
   - Retrieves the top 100 stars within 0.1° of RA = 36.033303 and Dec = +18.25795
   - Saves the output as `data/gaia_data.csv`

2. **`view_data`**  
   - Reads the CSV data and renders it as an HTML table using `csv`
   - Saves the output as `data/viewer.html`

3. **`all`**  
   - Final rule that depends on `data/viewer.html`

---

## 🚀 Running the Pipeline

From the root of the project:

```bash
cd workflow
snakemake --cores all
```
