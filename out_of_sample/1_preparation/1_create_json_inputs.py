import pandas as pd
import json
import os
import argparse

def create_af3_json_inputs(csv_file: str, output_dir: str):
    """
    Reads a CSV file (with at least columns: 'pdb_id', 'receptor_sequence', 'ligand_id', 'ligand_smiles')
    and generates a JSON file for each row in the AlphaFold 3 input format.

    The resulting JSON structure for each entry is:
    {
        "name": <pdb_id>_<ligand_id>,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": <receptor_sequence>
                }
            },
            {
                "ligand": {
                    "id": "B",
                    "smiles": <ligand_smiles>
                }
            }
        ],
        "modelSeeds": [4,9,8,31,20],
        "dialect": "alphafold3",
        "version": 1
    }

    Each file is named <pdb_id>_<ligand_id>_input.json, saved under the specified output directory.
    """

    # Read CSV into a DataFrame
    df = pd.read_csv(csv_file)

    # Create output directory (if it doesn't exist)
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each row in the DataFrame
    for idx, row in df.iterrows():
        entry_id = row['pdb_id']
        protein_sequence = row['receptor_sequence']
        ligand_smiles = row['ligand_smiles']
        ligand_id = row['ligand_id'].lower()

        # Construct the JSON dictionary
        af3_input = {
            "name": entry_id+"_"+ligand_id,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": protein_sequence
                    }
                },
                {
                    "ligand": {
                        "id": "B",
                        "smiles": ligand_smiles
                    }
                }
            ],
            "modelSeeds": [4,9,8,31,20],
            "dialect": "alphafold3",
            "version": 1
        }

        # Define the output JSON filename
        output_file = os.path.join(output_dir, f"{entry_id}_{ligand_id}_input.json")

        # Write the JSON file
        with open(output_file, 'w') as f:
            json.dump(af3_input, f, indent=2)

        # print(f"Created: {output_file}")

if __name__ == "__main__":
    # Use argparse to handle command-line arguments
    parser = argparse.ArgumentParser(
        description="Generate AlphaFold 3 JSON input files from a CSV containing protein-ligand data."
    )
    parser.add_argument(
        "-c",
        "--csv",
        required=True,
        help="Path to the input CSV file (must contain columns: 'pdb_id', 'receptor_sequence', 'ligand_id', and 'ligand_smiles')."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="af_input",
        help="Directory where the JSON files will be saved (default: af_input)."
    )

    args = parser.parse_args()

    create_af3_json_inputs(args.csv, args.outdir)
