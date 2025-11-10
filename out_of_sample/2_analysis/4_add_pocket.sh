#!/bin/bash

set -euo pipefail

BASE_DIR="finished_outputs"

for SUBDIR in "$BASE_DIR"/*/; do
	[ -d "$SUBDIR" ] || continue

	af3_model="${SUBDIR}af3_model.pdb"
	af3_pocket="${SUBDIR}af3_pocket.pdb"
	af3_out="${SUBDIR}af3_model_pocket_added.pdb"

	if [ ! -f $af3_pocket ]; then
		continue
	fi

	ref_complex="${SUBDIR}ref_complex.pdb"
	ref_pocket="${SUBDIR}ref_pocket.pdb"
	ref_out="${SUBDIR}ref_complex_pocket_added.pdb"

	{
		grep -v '^END' "$af3_model"
	
		echo 'TER'
		echo 'PKT          20                1000         pkt'

		cat "$af3_pocket"
	} > "$af3_out"

    {
        grep -v '^END' "$ref_complex"
    
        echo 'TER'
        echo 'PKT          20                1000         pkt'

        cat "$ref_pocket"
    } > "$ref_out"
done

echo "Done"
