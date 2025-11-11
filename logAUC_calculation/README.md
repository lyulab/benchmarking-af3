# Calculate logAUC

## Installation
Install two conda environments `py3d` and `dock37_py27` via yaml files `./envs/dock37_py27.yml` and `./envs/dock37_py27.yml`

## Run analysis code
0. (optional) Set up conda `source "$CONDA_SH"`
1. [under `py3d`] Run `split_running_sum.py` to separate the raw data file containing `metric` values for all targets into individual folders organized by target. The output follows the structure `scoring_metric_name/target_name/split.csv`. For example, running it on `DUDEZ/running_sum/iptm_complex_af3_running_sum.csv` from the data deposit will produce `scoring_metric_name/target_name/split.csv` files for each DUDEZ target.
2. Run `compute_auc_from_splits.sh` to calculate AUC and logAUC values for each DUDEZ target by setting root to the output of step 1.
For example, `bash 'compute_auc_from_splits.sh' --dry-run --root ./DUDEZ_PB_valid/docking_score` where `docking_score` contains subfolders `target_name1`, `target_name2`...
3. [under `py3d`] Run `collect_auc_summary.py` to generate a summary `.csv` file containing AUC, logAUC, number of ligands, and number of decoys for each target.
For example,
```bash
cd ./DUDEZ_PB_valid
python collect_auc_summary.py
```
