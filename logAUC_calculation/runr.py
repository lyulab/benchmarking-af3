import csv
import sys
import subprocess

if len(sys.argv) < 2:
    print("Usage: python convert_csv_to_txt.py <input_csv_file>")
    sys.exit(1)
csv_file = sys.argv[1]
output_file = 'extract_all.txt'

with open(csv_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.reader(infile)
    # Skip header
    next(reader)

    for row in reader:
        compound_id = row[1]
        docking_score = row[2]

        # Create output row. Column 3: compound_id; Column 22: docking_score
        output_row = ['0', '0', compound_id] + ['0'] * 18 + [docking_score]
        outfile.write('\t'.join(output_row) + '\n')

print(f"Conversion complete, output saved to {output_file}")

print("Generating extract_all.sort.uniq.txt...")

subprocess.run(
    #"sort -k3,3 -k22,22nr extract_all.txt | awk '!seen[$3]++' | sort -k22,22nr > extract_all.sort.uniq.txt",
    "sort -k22,22nr extract_all.txt > extract_all.sort.uniq.txt",
    shell=True,
    check=True
)

print("File extract_all.sort.uniq.txt generated successfully.")
