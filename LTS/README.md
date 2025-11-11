# Large-Scale Experimental Testing Sets

## Run Alphafold3

You can run Alphafold3 via command 
```
python run_alphafold.py \
   --input_dir=/root/af_inputs \
   --model_dir=/root/models \
   --output_dir=/root/af_output
```
where input_path is folder that contains the json configurations like `./input_with_template_examples`.

Please note that you need to set up your docker environment. See further Alphafold3 documentation [here](https://github.com/google-deepmind/alphafold3) about how to set it up.