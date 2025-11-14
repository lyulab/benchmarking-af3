# Using AF3 to predict on out-of-sample complexes

## 🧰 Requirements

Before running the analyses, install the following software and dependencies:

| Software | Description | Link |
|-----------|--------------|------|
| **AlphaFold3** | Structure prediction | [github.com/google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3) |
| **APoc** | Alignment of protein pockets | [sites.gatech.edu/cssb/apoc](https://sites.gatech.edu/cssb/apoc) |
| **DockRMSD** | Ligand RMSD computation | [aideepmed.com/DockRMSD](https://aideepmed.com/DockRMSD) |

You will also need to download the **BioLiP2** database: 🔗 [https://zhanggroup.org/BioLiP](https://zhanggroup.org/BioLiP). 

Instructions for downloading this database using the provided Perl script are available [here](https://aideepmed.com/BioLiP/weekly.html). Note that the script retrieves both the _redundant_ and _non-redundant_ datasets, where the non-redundant set is “a subset of the redundant dataset by protein sequence clustering at 90% sequence identity” ([BioLiP website](https://aideepmed.com/BioLiP/weekly.html)). For our dataset curation, we use the full redundant dataset. Because the Perl script downloads both versions, we remove the non-redundant protein and ligand directories after download, then rename the remaining directories to `protein` and `ligand` within the automatically generated `BioLiP_updated_set` folder.

Although this workflow provides an end-to-end process for dataset curation and initial benchmarking, we are unable to distribute the full dataset in this repository. Instead, we include demo CSV files and sample experimental structures organized to mirror the BioLiP directory structure. We also include a sample set of AF3 outputs only containing the top-ranked mmcif files.


## ⚙️ Environment Setup

Install [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) before proceeding.

This work mainly leverages [Biopython](https://biopython.org/) (an open-source bioinformatics toolkit) and [RDKit](https://www.rdkit.org/) (an open-source cheminformatics toolkit).

Run the following commands to create and configure the environment:
```bash
conda create -n af3-benchmarking python=3.10
conda activate af3-benchmarking
conda install -c conda-forge biopython rdkit
conda install pandas
```


## 📂 Directory Structure

```
out_of_sample/
├── preparation/
│   ├── 0_DatasetConstruction.ipynb
│   └── 1_create_json_inputs.py
├── 2_copy_biolip_files.py
├── 3_find_pocket_residues.py
├── 4_add_pocket.sh
├── 5_run_apoc.sh
├── 6_coord_transform.sh
├── 7_save_mol2_array.sh
├── 8_run_dockrmsd.sh
├── 9_save_metrics.py
├── finished_outputs/
└── BioLiP_updated_set/
    ├── ligand/
    └── receptor/
```


## 📝 Workflow and Script Descriptions

### 1️⃣ Preparation

#### `0_DatasetConstruction.ipynb`
Processes data downloaded from [RCSB](https://www.rcsb.org/) and cross-references with the BioLip2 database.  
Saves a CSV file containing information for each protein–ligand complex, including receptor sequence and ligand SMILES strings.

#### `1_create_json_inputs.py`
Generates AlphaFold3 JSON input files from a CSV. Outputs are saved in the `af_input` folder.

**Usage:**
```bash
python preparation/1_create_json_inputs.py -c preparation/af3_inputs_demo.csv -o af_input
```


### 2️⃣ Analysis

> Once you run AF3 and make sure all predictions ran successfully (or filter out failed ones), rename the directory containing these AF3 outputs to `finished_outputs`.

#### `2_copy_biolip_files.py`
Copies the corresponding experimental protein and ligand structures from BioLip into the `finished_outputs` directory.

- Update the hardcoded paths to your BioLip directories before running if needed.

**Usage:**
```bash
python 2_copy_biolip_files.py
```

#### `3_find_pocket_residues.py`
Identifies binding pocket residues within 5 Å of the experimental ligand and maps them to the AF3 predicted structure.

**Usage:**
```bash
python 3_find_pocket_residues.py
```

#### `4_add_pocket.sh`
Appends user-defined pocket selections to create a combined PDB file (used as APoc input).

**Usage:**
```bash
./4_add_pocket.sh
```

#### `5_run_apoc.sh`
Runs the **APoc** executable for all AF3 outputs and generates alignment reports.

- Update the `APOC_BIN` variable to the path to your APoc installation.

**Usage:**
```bash
./5_run_apoc.sh
```

#### `6_coord_transform.sh`
Applies transformation matrices from APoc to align predicted and experimental structures (pocket-aligned). References `align_structures.py`.

**Usage:**
```bash
./6_coord_transform.sh
```

#### `7_save_mol2_array.sh`
Extracts and saves ligand `.mol2` files for RMSD computation. References `extract_ligand.py`.

**Usage:**
```bash
./7_save_mol2_array.sh
```

#### `8_run_dockrmsd.sh`
Runs **DockRMSD** to calculate ligand RMSD between predicted and reference structures.

- Update the `DOCKRMSD_BIN` variable to the path to your DockRMSD installation.

**Usage:**
```bash
./8_run_dockrmsd.sh
```

#### `9_save_metrics.py`
Aggregates all computed metrics into a summary CSV file.

**Usage:**
```bash
python 9_save_metrics.py
```

### 📊 Output

The final results include:
- Pocket-aligned protein-ligand structures 
- A summary CSV of all performance metrics
