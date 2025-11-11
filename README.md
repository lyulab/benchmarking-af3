# benchmarking-af3
This repository contains the code used for our paper titled "Benchmarking AlphaFold3 for structure-based ligand discovery"

## DUDEz and LTS
The `DUDEz` and `LTS` folers contain example templates to run Alphafold3.

## Enrichment analysis
The `logAUC_calculation` folder contains pipeline to calculate enrichment metrics like logAUC.

## Out-of-sample Analysis

### Installation

The following software is needed before running our analyses:
- APoc (Alignment of Pockets): https://sites.gatech.edu/cssb/apoc
- DockRMSD (Ligand RMSD): https://aideepmed.com/DockRMSD
- Biopython (Python package): https://biopython.org

The BioLip2 database will also need to be downloaded: https://zhanggroup.org/BioLiP

The `out_of_sample` folder contains two subfolders: `1_preparation` and `2_analysis`. 

In the `1_preparation` folder, there's:
- `0_DatasetConstruction.ipynb`: this notebook can be used to process data downloaded from RCSB and cross-reference with the BioLip2 database. A CSV file can be created containing information about the protein-ligand complex with the necessary information for AlphaFold3, including but not limited to the receptor sequence and ligand SMILES.
- `1_create_json_inputs.py`: script that takes in a CSV and generates AF3 JSON inputs into a folder.

In the `2_analysis` folder, there's:
- `2_copy_biolip_files.py`: copies the necessary experimental protein and ligand structures from the BioLiP database into the corresponding AF3 output directory.
- `3_find_pocket_residues.py`: performs a binding pocket residue search 5Å around the experimental ligand, and then performs a mapping to the AF3 mmcif file to obtain the corresponding pocket residues.
- `4_add_pocket.sh`: appends the user-defined pocket selection into a new pdb containing both the full structure and the added pocket. This file serves as an input to perform APoc.
- `5_run_apoc.sh`: runs the APoc executable in every AF3 output directory and saves a txt file.
- `6_coord_transform.sh`: uses the transformation coordinates generated from APoc and aligns the experimental and predicted structures (pocket-aligned). References the align_structures.py file.
- `7_save_mol2_array.sh`: creates and saves ligand mol2's, to be used for ligand RMSD calculations. This script has been parallelized to perform multiple operations at once. References the extract_ligand.py file.
- `8_run_dockrmsd.sh`: runs the DockRMSD executable and saves a txt file containing the ligand RMSD between the AF3 and reference structure.
- `9_save_metrics.py`: saves all metrics into a CSV file.

## Perspective
The `perspective` folder contains example to run Boltz prediction and experimental data correlation calculation.
