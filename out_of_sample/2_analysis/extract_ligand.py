import argparse
from pymol import cmd, finish_launching

def main():
	parser = argparse.ArgumentParser(
		description = 'Extract a ligand from a PDB complex and save as mol2 without common metals and halogens'
	)
	parser.add_argument(
		'-i', '--input', required=True,
		help='Input PDB filename (e.g. my_complex.pdb)'
	)
	parser.add_argument(
		'-r', '--resname', required=True,
		help='Ligand residue name (e.g. LIG)'
	)
	parser.add_argument(
		'-o', '--output', required=True,
		help='Output MOL2 filename (e.g. lig.mol2)'
	)
	args = parser.parse_args()

	# launch PyMOL quietly, no GUI
	finish_launching(['pymol', '-qc'])

	# load PDB
	cmd.load(args.input, 'complex')

	# extract ligand
	cmd.extract('lig', f'resn {args.resname}')

	# remove common metals and halogens
	metals = ['ZN', 'MG', 'CA', 'FE', 'MN', 'CO', 'NI', 'CU', 'NA', 'K', 'PD', 'CD', 'HG', 'I', 'MO', 'RU', 'AG', 'CO', 'PT', 'AU', 'CL', 'F', 'BR']
	metal_sel = ' or '.join(f'elem {m}' for m in metals)
	cmd.remove(f'lig and ({metal_sel})')

	# save cleaned ligand as mol2
	cmd.save(args.output, 'lig', format='mol2')

	# exit
	cmd.quit()

if __name__ == '__main__':
	main()
