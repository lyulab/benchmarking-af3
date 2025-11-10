import argparse
import pandas as pd
import shutil
import glob
from pathlib import Path
import sys

def main():
    p = argparse.ArgumentParser(
        description="Copy BioLiP receptor & ligand PDBs into AF3 folders"
    )
    p.add_argument(
        "-c", "--csv", required=True,
        help="Path to CSV file (with columns: pdb_id, ligand_chain, ligand_id, receptor_chain, …)"
    )
    args = p.parse_args()

    # PARAMETERS
    receptor_dir = Path(
        "/lustre/fs6/lyu_lab/scratch/adavasam/"
        "BioLip2_DB/BioLiP_updated_set/receptor"
    )
    ligand_dir = Path(
        "/lustre/fs6/lyu_lab/scratch/adavasam/"
        "BioLip2_DB/BioLiP_updated_set/ligand"
    )
    af3_root = Path.cwd() / "finished_outputs"
    no_exp_dir = af3_root.parent / "no_experimental_structures"
    no_exp_dir.mkdir(parents=True, exist_ok=True)

    # Load the CSV
    try:
        df = pd.read_csv(args.csv, dtype=str)
    except Exception as e:
        sys.exit(f"Error reading CSV: {e}")

    for i, row in df.iterrows():
        pdb_id         = row["pdb_id"].lower()
        ligand_id      = row["ligand_id"]
        ligand_chain   = row["ligand_chain"]
        receptor_chain = row["receptor_chain"]

        subdir_name = f"{pdb_id}_{ligand_id.lower()}"
        subdir = af3_root / subdir_name
        if not subdir.is_dir():
            print(f"[SKIP] AF3 output missing: {subdir_name}")
            continue

        # find receptor
        receptor_file = receptor_dir / f"{pdb_id}{receptor_chain}.pdb"
        if not receptor_file.is_file():
            print(f"[MISSING] {receptor_file.name}")
            missing = True
        else:
            missing = False

        # find ligand (first match)
        ligand_pattern = ligand_dir / f"{pdb_id}_{ligand_id}_{ligand_chain}_*.pdb"
        matches = sorted(glob.glob(str(ligand_pattern)))
        if matches:
            ligand_file = Path(matches[0])
        else:
            print(f"[MISSING] ligand globs to {ligand_pattern.name}")
            missing = True

        # if anything is missing → move entire folder
        if missing:
            target = no_exp_dir / subdir_name
            print(f"[MOVE] {subdir_name} → no_experimental_structures/")
            shutil.move(str(subdir), str(target))
            continue

        # copy receptor & ligand into place
        dst_receptor = subdir / "ref_prot.pdb"
        dst_ligand   = subdir / "ref_lig.pdb"
        shutil.copy2(str(receptor_file), str(dst_receptor))
        shutil.copy2(str(ligand_file),   str(dst_ligand))
        # print(f"[COPY] {receptor_file.name} → {subdir_name}/ref_prot.pdb")
        # print(f"[COPY] {ligand_file.name}   → {subdir_name}/ref_lig.pdb")

if __name__ == "__main__":
    main()

