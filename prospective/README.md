# Prospective

## Run Boltz

You can run Boltz via command `boltz predict input_path` where input_path is `boltz_input_example.yaml`.

See further Boltz documentation [here](https://github.com/jwohlwend/boltz)

## Correlation
1. Create a conda environment using `pip install -r requirements.txt`
2. Load `sigma2_prospective/cpm_ic50.csv` which is stored in the separate data deposit and run the Spearman correlation between the experimental log IC50 value and predicted affinity value and affinity probability.
