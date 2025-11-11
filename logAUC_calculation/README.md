# Calculate logAUC
1. Run `split_running_sum.py` to separate the raw data file containing `metric` values for all targets into individual folders organized by target. The output follows the structure `scoring_metric_name/target_name/split.csv`. For example, running it on `DUDEZ/running_sum/iptm_complex_af3_running_sum.csv` from the data deposit will produce `scoring_metric_name/target_name/split.csv` files for each DUDEZ target.
2. Run `compute_auc_from_splits.sh` to calculate AUC and logAUC values for each DUDEZ target.
3. Run `collect_auc_summary.py` to generate a summary `.csv` file containing AUC, logAUC, number of ligands, and number of decoys for each target.
