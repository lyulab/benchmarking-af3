### 🧰 Requirements

Before running the analyses, install the following software and dependencies:

| Software | Description | Link |
|-----------|--------------|------|
| **AlphaFold3** | Protein structure prediction | [github.com/google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3) |
| **APoc** | Alignment of protein pockets | [sites.gatech.edu/cssb/apoc](https://sites.gatech.edu/cssb/apoc) |
| **DockRMSD** | Ligand RMSD computation | [aideepmed.com/DockRMSD](https://aideepmed.com/DockRMSD) |

You will also need to download the **BioLip2** database: 🔗 [https://zhanggroup.org/BioLiP](https://zhanggroup.org/BioLiP)


### ⚙️ Environment Setup

Install [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) before proceeding.

This work uses the same environment from [AI-driven Structure-enabled Antiviral Platform (ASAP)](https://asapdiscovery.readthedocs.io/), 
and leverages [Biopython](https://biopython.org) (Python bioinformatics toolkit) for reproducible benchmarking.

Run the following commands to create and configure the environment:
```bash
mamba create -n af3-benchmarking python=3.10
mamba activate af3-benchmarking
mamba install -c conda-forge asapdiscovery
mamba install -c openeye openeye-toolkits
```


### 📂 Directory Structure

```
out_of_sample/
├── 1_preparation/
│   ├── 0_DatasetConstruction.ipynb
│   └── 1_create_json_inputs.py
└── 2_analysis/
    ├── 2_copy_biolip_files.py
    ├── 3_find_pocket_residues.py
    ├── 4_add_pocket.sh
    ├── 5_run_apoc.sh
    ├── 6_coord_transform.sh
    ├── 7_save_mol2_array.sh
    ├── 8_run_dockrmsd.sh
    └── 9_save_metrics.py
```


### 📝 Workflow and Script Descriptions

#### 1️⃣ Preparation

##### `0_DatasetConstruction.ipynb`
Processes data downloaded from [RCSB](https://www.rcsb.org/) and cross-references with the BioLip2 database.  
Saves a CSV file containing information for each protein–ligand complex, including receptor sequence and ligand SMILES strings.

##### `1_create_json_inputs.py`
Generates AlphaFold3 JSON input files from a CSV. Outputs are saved in the `af_input` folder.

**Usage:**
```bash
python 1_create_json_inputs.py -c af3_inputs.csv
```


#### 2️⃣ Analysis

> All scripts in this folder should be placed in the same directory as your AF3 output folder, renamed to `finished_outputs`.

##### `2_copy_biolip_files.py`
Copies the corresponding experimental protein and ligand structures from BioLip into the `finished_outputs` directory.

- Update the hardcoded paths to your BioLip directories before running.

**Usage:**
```bash
python 2_copy_biolip_files.py
```

##### `3_find_pocket_residues.py`
Identifies binding pocket residues within 5 Å of the experimental ligand and maps them to the AF3 predicted structure.

**Usage:**
```bash
python 3_find_pocket_residues.py
```

##### `4_add_pocket.sh`
Appends user-defined pocket selections to create a combined PDB file (used as APoc input).

**Usage:**
```bash
./4_add_pocket.sh
```

##### `5_run_apoc.sh`
Runs the **APoc** executable for all AF3 outputs and generates alignment reports.

**Usage:**
```bash
./5_run_apoc.sh
```

##### `6_coord_transform.sh`
Applies transformation matrices from APoc to align predicted and experimental structures (pocket-aligned). References `align_structures.py`.

**Usage:**
```bash
./6_coord_transform.sh
```

##### `7_save_mol2_array.sh`
Extracts and saves ligand `.mol2` files for RMSD computation. References `extract_ligand.py`.

**Usage:**
```bash
./7_save_mol2_array.sh
```

##### `8_run_dockrmsd.sh`
Runs **DockRMSD** to calculate ligand RMSD between predicted and reference structures.

**Usage:**
```bash
./8_run_dockrmsd.sh
```

##### `9_save_metrics.py`
Aggregates all computed metrics into a summary CSV file.

**Usage:**
```bash
python 9_save_metrics.py
```

### 📊 Output

The final results include:
- Pocket-aligned protein-ligand structures 
- A summary CSV of all performance metrics
