# align_structures.py
#
# Usage:
#     python align_structures.py /path/to/finished_outputs/subdir
#
# Notes:
# - Reads APoc text output and extracts ONLY the "Pocket alignment" rotation/translation.
# - Writes ref_complex_pocket_aligned.pdb if a pocket transform is found.

import os, sys, re
import numpy as np
from Bio.PDB import PDBParser, PDBIO

HEADER_RE = re.compile(r'^\s*i\s+t\(i\)\s+u\(i,1\)\s+u\(i,2\)\s+u\(i,3\)')

def parse_pocket_transformation(apoc_path):
    """
    Look for the 'Pocket alignment' block and parse the 3x3 rotation and 3x1 translation.
    Returns (U_pocket, t_pocket) or (None, None) if no pocket transform is present.
    """
    in_pocket_section = False

    with open(apoc_path, 'r') as f:
        for line in f:
            # Enter the pocket section
            if 'Pocket alignment' in line:
                in_pocket_section = True
                continue

            # Once in pocket section, look for the rotation-matrix header and grab next 3 rows
            if in_pocket_section and HEADER_RE.match(line):
                rows = []
                for _ in range(3):
                    try:
                        rows.append(next(f).split())
                    except StopIteration:
                        return None, None  # incomplete table

                # rows[i] looks like: [i, t(i), u(i,1), u(i,2), u(i,3)]
                try:
                    t_pocket = np.array([float(r[1]) for r in rows], dtype=float)
                    U_pocket = np.array([[float(r[j]) for j in (2, 3, 4)] for r in rows], dtype=float)
                except (IndexError, ValueError):
                    return None, None

                return U_pocket, t_pocket

    # If we never found a pocket alignment block with a rotation table
    return None, None

def transform_with_biopython(ref_pdb, out_pdb, U, t):
    """
    Applies x' = U x + t to every atom in ref_pdb and writes out_pdb.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("ref_pdb", ref_pdb)
    for atom in struct.get_atoms():
        coord = atom.get_coord()
        atom.set_coord(U.dot(coord) + t)

    io = PDBIO()
    io.set_structure(struct)
    io.save(out_pdb)

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    subdir = sys.argv[1]
    apoc_txt = os.path.join(subdir, "apoc_output.txt")
    ref_pdb = os.path.join(subdir, "ref_complex.pdb")
    out_pdb_pocket = os.path.join(subdir, "ref_complex_pocket_aligned.pdb")

    if not os.path.isfile(apoc_txt):
        print(f"[WARN] APoc output not found: {apoc_txt}")
        sys.exit(0)

    if not os.path.isfile(ref_pdb):
        print(f"[WARN] Reference PDB not found: {ref_pdb}")
        sys.exit(0)

    U_pocket, t_pocket = parse_pocket_transformation(apoc_txt)

    if U_pocket is None or t_pocket is None:
        print("[INFO] No successful APoc pocket alignment found; nothing to write.")
        sys.exit(0)

    transform_with_biopython(ref_pdb, out_pdb_pocket, U_pocket, t_pocket)
    # print(f"[OK] Wrote pocket-aligned PDB: {out_pdb_pocket}")

if __name__ == "__main__":
    main()

