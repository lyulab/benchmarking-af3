#!/bin/bash

BASE_DIR="finished_outputs"
DOCKRMSD_BIN="/lustre/fs6/lyu_lab/scratch/kmenon01/software/bin/DockRMSD"

for SUBDIR in "$BASE_DIR"/*/; do
    [ -d "$SUBDIR" ] || continue

	# run pocket-aligned ligand rmsd only if that file exists
	if [[ -f "${SUBDIR}ref_lig_pocket_aligned.mol2" ]]; then
		$DOCKRMSD_BIN "${SUBDIR}af3_lig.mol2" "${SUBDIR}ref_lig_pocket_aligned.mol2" > "${SUBDIR}dockrmsd_pocket_output.txt"
	fi
done
echo "Done"
