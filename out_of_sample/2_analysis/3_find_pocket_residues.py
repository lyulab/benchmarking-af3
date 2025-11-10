import os
import sys
from typing import List, Dict, Set, Tuple
import warnings

from Bio.PDB import (
    PDBParser,
    MMCIFParser,
    PDBIO,
    Select,
    NeighborSearch
)
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


# ------------------------------ parsing & IO helpers ------------------------------

def _parse_structure(path: str, struct_id: str = "s"):
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    structure = parser.get_structure(struct_id, path)
    model = structure[0]  # first model
    return structure, model


def _write_selected_residues(structure, residues: Set, out_path: str):
    """Write only residues in `residues` from `structure`."""
    io = PDBIO()
    io.set_structure(structure)

    class PocketSelect(Select):
        def accept_residue(self, residue):
            return residue in residues

    io.save(out_path, PocketSelect())


# ------------------------------ pocket on reference ------------------------------

def _get_all_protein_chains(model):
    """Get all protein chains with their residues."""
    chains = {}
    ppb = PPBuilder()
    
    for chain in model:
        peptides = ppb.build_peptides(chain)
        if peptides:
            # Combine all peptides in this chain
            all_residues = []
            for pep in peptides:
                all_residues.extend(list(pep))
            if all_residues:
                chains[chain.id] = all_residues
    
    return chains


def _protein_residues(model, chain_id=None) -> List:
    """Return amino-acid Residues, optionally for specific chain."""
    ppb = PPBuilder()
    
    if chain_id:
        # Get residues from specific chain
        if chain_id in model:
            peptides = ppb.build_peptides(model[chain_id])
        else:
            print(f"Warning: Chain {chain_id} not found", file=sys.stderr)
            return []
    else:
        # Get all peptides from all chains
        peptides = ppb.build_peptides(model)
    
    if not peptides:
        return []
    
    # If multiple chains, combine all or take the longest
    if len(peptides) > 1:
        print(f"Note: Found {len(peptides)} polypeptide segments", file=sys.stderr)
    
    # Combine all peptides (maintains order)
    all_residues = []
    for pep in peptides:
        all_residues.extend(list(pep))
    
    return all_residues


def _compute_reference_pocket(res_complex_path: str, ligand_resname: str, 
                             cutoff: float = 5.0, verbose: bool = False) -> Tuple[Set, object, object, Dict]:
    """
    Build protein-ligand 'ref_complex' from file, then gather all AA residues within `cutoff` Å
    of any ligand atom (by residue name).
    Returns: (ref_pocket_residues_set, structure, model, stats_dict)
    """
    structure, model = _parse_structure(res_complex_path, "ref_complex")
    
    # Collect statistics
    stats = {
        'ligand_atoms': 0,
        'pocket_residues': 0,
        'chains_involved': set()
    }

    # gather ligand atoms (by residue name match)
    ligand_atoms = [
        atom
        for res in model.get_residues()
        if res.get_resname() == ligand_resname
        for atom in res.get_atoms()
    ]
    stats['ligand_atoms'] = len(ligand_atoms)
    
    if not ligand_atoms:
        print(f"Error: Ligand '{ligand_resname}' not found in {res_complex_path}", file=sys.stderr)
        return set(), structure, model, stats

    # spatial neighbor search over all atoms
    ns = NeighborSearch(list(model.get_atoms()))
    pocket: Set = set()
    
    for latom in ligand_atoms:
        for neighbor in ns.search(latom.coord, cutoff):
            res = neighbor.get_parent()
            if res.get_id()[0] == " " and is_aa(res):
                pocket.add(res)
                # Track which chain this residue belongs to
                chain = res.get_parent()
                stats['chains_involved'].add(chain.id)
    
    stats['pocket_residues'] = len(pocket)
    
    if verbose:
        print(f"Reference pocket: {stats['pocket_residues']} residues from chain(s) {stats['chains_involved']}")
        print(f"  Found {stats['ligand_atoms']} ligand atoms")
    
    return pocket, structure, model, stats


# ------------------------------ sequence alignment mapping ------------------------------

def _sequence_and_residues(model, chain_id=None) -> Tuple[str, List]:
    """Return (protein sequence as string, list of Residue objects in that sequence order)."""
    residues = _protein_residues(model, chain_id)
    if not residues:
        return "", []
    
    THREE_TO_ONE = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Alternative names sometimes seen
        'ASX': 'B',  # ASP or ASN
        'GLX': 'Z',  # GLU or GLN
        'XLE': 'J',  # LEU or ILE
        'SEC': 'U',  # Selenocysteine
        'PYL': 'O',  # Pyrrolysine
    }

    seq_parts = []
    valid_residues = []
    
    for res in residues:
        resname = res.get_resname().strip().upper()
        
        if resname in THREE_TO_ONE:
            # Known amino acid
            aa = THREE_TO_ONE[resname]
            seq_parts.append(aa)
            valid_residues.append(res)
        else:
            # Unknown/non-standard amino acid - use 'X'
            seq_parts.append('X')
            valid_residues.append(res)
            print(f"Warning: Non-standard residue {resname} at {res.get_id()}", file=sys.stderr)
    
    seq = ''.join(seq_parts)
    return seq, valid_residues


def _build_ref_to_model_index_map(ref_seq: str, model_seq: str, verbose: bool = False) -> Dict[int, int]:
    """
    Map residue indices (0-based positions in ungapped ref_seq) to indices in model_seq,
    using a global alignment (PairwiseAligner + BLOSUM62).
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(ref_seq, model_seq)
    
    if not alignments:
        print("Error: No alignment found between sequences", file=sys.stderr)
        return {}
    
    alignment = alignments[0]  # best alignment
    
    if verbose:
        print(f"\nAlignment score: {alignment.score}")
        print(f"Reference length: {len(ref_seq)}")
        print(f"Model length: {len(model_seq)}")
        
        # Calculate identity
        aligned_ref = str(alignment).split('\n')[0]
        aligned_model = str(alignment).split('\n')[2]
        matches = sum(1 for a, b in zip(aligned_ref, aligned_model) 
                     if a == b and a != '-')
        identity = matches / max(len(ref_seq), len(model_seq)) * 100
        print(f"Sequence identity: {identity:.1f}%")
    
    coords = alignment.coordinates
    mapping: Dict[int, int] = {}

    for c in range(1, coords.shape[1]):
        ra0, ra1 = int(coords[0, c - 1]), int(coords[0, c])
        rb0, rb1 = int(coords[1, c - 1]), int(coords[1, c])
        
        # Only map positions that align to actual residues on both sides
        if (ra1 - ra0) > 0 and (rb1 - rb0) > 0:
            block_len = min(ra1 - ra0, rb1 - rb0)
            for i in range(block_len):
                mapping[ra0 + i] = rb0 + i
    
    return mapping


def _map_pocket_residues_to_model(
    ref_pocket: Set,
    ref_model,
    model_model,
    verbose: bool = False
) -> Tuple[Set, Dict]:
    """
    Using sequence alignment, map reference pocket residues to the equivalent residues in model.
    Returns: (model_pocket_set, mapping_stats_dict)
    """
    # Get chains involved in the pocket
    pocket_chains = {}
    for res in ref_pocket:
        chain_id = res.get_parent().id
        if chain_id not in pocket_chains:
            pocket_chains[chain_id] = []
        pocket_chains[chain_id].append(res)
    
    if verbose and len(pocket_chains) > 1:
        print(f"\nPocket spans {len(pocket_chains)} chain(s): {list(pocket_chains.keys())}")
    
    model_pocket: Set = set()
    stats = {
        'ref_pocket_size': len(ref_pocket),
        'mapped_residues': 0,
        'unmapped_residues': [],
        'chain_mapping': {}
    }
    
    # Try to map each chain separately
    for chain_id, chain_residues in pocket_chains.items():
        ref_seq, ref_res_seq = _sequence_and_residues(ref_model, chain_id)
        
        # Try to find corresponding chain in model
        # First try same chain ID, then try all chains
        model_seq, model_res_seq = _sequence_and_residues(model_model, chain_id)
        
        if not model_seq:
            # Try to find best matching chain
            if verbose:
                print(f"Chain {chain_id} not found in model, searching for best match...")
            
            best_chain = None
            best_score = -float('inf')
            
            for test_chain in model_model:
                test_seq, test_res = _sequence_and_residues(model_model, test_chain.id)
                if test_seq:
                    aligner = PairwiseAligner()
                    aligner.mode = "global"
                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    score = aligner.score(ref_seq, test_seq)
                    if score > best_score:
                        best_score = score
                        best_chain = test_chain.id
                        model_seq = test_seq
                        model_res_seq = test_res
            
            if best_chain and verbose:
                print(f"  Best match: chain {best_chain} (score: {best_score:.1f})")
        
        if not ref_seq or not model_seq:
            print(f"Warning: Could not build sequences for chain {chain_id}", file=sys.stderr)
            continue
        
        # Build residue to index map for this chain
        ref_idx_map = {res: i for i, res in enumerate(ref_res_seq)}
        
        # Build index mapping via alignment
        r2m = _build_ref_to_model_index_map(ref_seq, model_seq, verbose)
        
        # Map pocket residues
        chain_mapped = 0
        chain_unmapped = []
        
        for res in chain_residues:
            if res in ref_idx_map:
                ref_idx = ref_idx_map[res]
                if ref_idx in r2m:
                    midx = r2m[ref_idx]
                    if 0 <= midx < len(model_res_seq):
                        model_pocket.add(model_res_seq[midx])
                        chain_mapped += 1
                    else:
                        chain_unmapped.append(f"{res.get_resname()}{res.get_id()[1]}")
                else:
                    chain_unmapped.append(f"{res.get_resname()}{res.get_id()[1]}")
        
        stats['chain_mapping'][chain_id] = {
            'mapped': chain_mapped,
            'unmapped': len(chain_unmapped),
            'unmapped_list': chain_unmapped
        }
        stats['mapped_residues'] += chain_mapped
        stats['unmapped_residues'].extend(chain_unmapped)
    
    if verbose:
        print(f"\nMapping summary:")
        print(f"  Reference pocket: {stats['ref_pocket_size']} residues")
        print(f"  Successfully mapped: {stats['mapped_residues']} residues")
        print(f"  Could not map: {len(stats['unmapped_residues'])} residues")
        if stats['unmapped_residues']:
            print(f"  Unmapped: {', '.join(stats['unmapped_residues'][:10])}")
            if len(stats['unmapped_residues']) > 10:
                print(f"    ... and {len(stats['unmapped_residues']) - 10} more")
    
    return model_pocket, stats


# ------------------------------ main folder processor ------------------------------

def process_finished_outputs(base_dir="finished_outputs", cutoff=5.0, verbose=False):
    """Process all subdirectories with optional verbose output for debugging."""
    
    summary = []
    
    for entry in os.scandir(base_dir):
        if not entry.is_dir():
            continue
        subdir = entry.path
        name = entry.name
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing: {name}")
            print('='*60)

        # parse pdb and ligand IDs
        try:
            pdbid, ligandid_raw = name.split('_', 1)
        except ValueError:
            print(f"Skipping '{name}': not in <pdbid>_<ligandid> format", file=sys.stderr)
            continue

        # ligand residue name normalization
        ligandid = ligandid_raw[:3] if len(ligandid_raw) == 5 else ligandid_raw
        true_ligand_resname = ligandid.upper()

        # ---------- (A) predicted AF3 model ----------
        model_fname = f"{pdbid}_{ligandid_raw}_model.cif"
        model_path = os.path.join(subdir, model_fname)
        af3_structure = af3_model = None
        
        if os.path.isfile(model_path):
            af3_structure, af3_model = _parse_structure(model_path, "af3_struct")
            out_pdb_copy = os.path.join(subdir, "af3_model.pdb")
            io = PDBIO()
            io.set_structure(af3_structure)
            io.save(out_pdb_copy)
        else:
            print(f"Warning: '{model_fname}' not found in {subdir}", file=sys.stderr)

        # ---------- (B) reference complex ----------
        ref_prot = os.path.join(subdir, "ref_prot.pdb")
        ref_lig = os.path.join(subdir, "ref_lig.pdb")

        if os.path.isfile(ref_prot) and os.path.isfile(ref_lig):
            combined = os.path.join(subdir, "ref_complex.pdb")
            with open(combined, 'w') as w:
                for fname in (ref_prot, ref_lig):
                    with open(fname) as r:
                        w.write(r.read())

            # 1) define pocket on the reference complex
            ref_pocket, ref_structure, ref_model, ref_stats = _compute_reference_pocket(
                combined, true_ligand_resname, cutoff=cutoff, verbose=verbose
            )
            
            if ref_pocket:
                # write reference pocket
                out_ref = os.path.join(subdir, "ref_pocket.pdb")
                _write_selected_residues(ref_structure, ref_pocket, out_ref)

                # 2) map that same pocket onto the AF3 model
                if af3_structure and af3_model:
                    af3_pocket, mapping_stats = _map_pocket_residues_to_model(
                        ref_pocket, ref_model, af3_model, verbose=verbose
                    )
                    
                    if af3_pocket:
                        out_af3 = os.path.join(subdir, "af3_pocket.pdb")
                        _write_selected_residues(af3_structure, af3_pocket, out_af3)
                        
                        # Save summary
                        summary.append({
                            'name': name,
                            'ref_pocket_size': ref_stats['pocket_residues'],
                            'af3_pocket_size': len(af3_pocket),
                            'difference': ref_stats['pocket_residues'] - len(af3_pocket),
                            'unmapped': len(mapping_stats['unmapped_residues'])
                        })
                        
                        # Write mapping report
                        report_path = os.path.join(subdir, "pocket_mapping_report.txt")
                        with open(report_path, 'w') as f:
                            f.write(f"Pocket Mapping Report for {name}\n")
                            f.write("="*50 + "\n\n")
                            f.write(f"Reference pocket: {ref_stats['pocket_residues']} residues\n")
                            f.write(f"AF3 pocket: {len(af3_pocket)} residues\n")
                            f.write(f"Successfully mapped: {mapping_stats['mapped_residues']} residues\n")
                            f.write(f"Could not map: {len(mapping_stats['unmapped_residues'])} residues\n\n")
                            
                            if mapping_stats['unmapped_residues']:
                                f.write("Unmapped residues:\n")
                                for res in mapping_stats['unmapped_residues']:
                                    f.write(f"  - {res}\n")
                    else:
                        print(f"Warning: could not map pocket to AF3 model in {subdir}", file=sys.stderr)
            else:
                print(f"Warning: empty pocket or ligand '{true_ligand_resname}' not found in {subdir}", file=sys.stderr)
        else:
            print(f"Warning: missing ref_prot.pdb or ref_lig.pdb in {subdir}", file=sys.stderr)
    
    # Print summary
    if summary and verbose:
        print(f"\n{'='*60}")
        print("SUMMARY")
        print('='*60)
        for item in summary:
            if item['difference'] != 0:
                print(f"{item['name']}: Ref={item['ref_pocket_size']}, "
                      f"AF3={item['af3_pocket_size']}, "
                      f"Diff={item['difference']}, "
                      f"Unmapped={item['unmapped']}")
    
    return summary


if __name__ == "__main__":
    # Add command line argument parsing
    import argparse
    
    parser = argparse.ArgumentParser(description='Map protein pockets between structures')
    parser.add_argument('--dir', default='finished_outputs', help='Base directory')
    parser.add_argument('--cutoff', type=float, default=5.0, help='Distance cutoff in Angstroms')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    results = process_finished_outputs(args.dir, args.cutoff, args.verbose)
