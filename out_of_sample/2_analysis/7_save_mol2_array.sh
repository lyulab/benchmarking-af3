#!/bin/bash
set -euo pipefail

BASE_DIR="${1:-finished_outputs}"
JOBS="${JOBS:-$(nproc)}"   # override with: JOBS=8 ./7_save_mol2_array.sh [BASE_DIR]

export BASE_DIR

find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -print0 |
xargs -0 -I{} -P "$JOBS" bash -c '
  set -euo pipefail
  SUBDIR="{}"
  [[ "${SUBDIR: -1}" == "/" ]] || SUBDIR="${SUBDIR}/"

  name=$(basename "$SUBDIR")
  resname="${name#*_}"
  # truncate resname to 3 chars if ligand is 5 chars long
  if [[ ${#resname} -eq 5 ]]; then
      resname="${resname:0:3}"
  fi
  resname="${resname^^}"

  python extract_ligand.py --input "${SUBDIR}${name}_model.cif" --resname LIG_B --output "${SUBDIR}af3_lig.mol2"

  if [[ -f "${SUBDIR}ref_complex_pocket_aligned.pdb" ]]; then
      python extract_ligand.py \
          --input  "${SUBDIR}ref_complex_pocket_aligned.pdb" \
          --resname "$resname" \
          --output "${SUBDIR}ref_lig_pocket_aligned.mol2"
  fi

  echo "[DONE] $name"
'

echo "Done"

