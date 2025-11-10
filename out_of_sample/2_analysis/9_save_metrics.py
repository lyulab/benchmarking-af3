import os
import re
import pandas as pd

BASE_DIR = "finished_outputs"

def parse_apoc_pocket_metrics(txt_path):
    """Return pocket RMSD, Seq identity, and PS-Score as a tuple of strings; empty strings if missing."""
    if not os.path.exists(txt_path):
        return "", "", ""
    
    with open(txt_path, encoding="utf-8", errors="ignore") as fh:
        text = fh.read()

    # Find the "Pocket alignment" section (robust to lots of > and spaces)
    idx = text.find("Pocket alignment")
    if idx == -1:
        return "", "", ""

    pocket_sec = text[idx:]
    
    # Parse RMSD
    rmsd_match = re.search(r"RMSD\s*=\s*([0-9]*\.?[0-9]+)", pocket_sec)
    rmsd = rmsd_match.group(1) if rmsd_match else ""
    
    # Parse Seq identity
    seq_id_match = re.search(r"Seq identity\s*=\s*([0-9]*\.?[0-9]+)", pocket_sec)
    seq_identity = seq_id_match.group(1) if seq_id_match else ""
    
    # Parse PS-Score
    ps_score_match = re.search(r"PS-score\s*=\s*([0-9]*\.?[0-9]+)", pocket_sec)
    ps_score = ps_score_match.group(1) if ps_score_match else ""
    
    return rmsd, seq_identity, ps_score

def parse_dockrmsd(txt_path):
    """Return RMSD score from DockRMSD output; blank if missing."""
    if not os.path.exists(txt_path):
        return ""
    with open(txt_path, encoding="utf-8", errors="ignore") as fh:
        text = fh.read()
    m = re.search(r"Calculated Docking RMSD:\s*([\d\.]+)", text)
    return m.group(1) if m else ""

def saveall():
    # Updated columns to include new APoc metrics
    column_names = ['complex_name', 'apoc_pocket', 'apoc_seq_identity', 'apoc_ps_score', 'dockrmsd_pocket']
    data_to_append = []

    for sub in sorted(os.listdir(BASE_DIR)):
        metrics_path = os.path.join(BASE_DIR, sub, 'metrics.dat')
        if not os.path.isfile(metrics_path):
            continue
        with open(metrics_path, 'r', encoding="utf-8", errors="ignore") as f:
            line = f.read().strip().split(',')
            # Ensure we have exactly the expected number of fields
            if len(line) != len(column_names):
                continue
            metrics_dict = {column_names[i]: line[i] for i in range(len(column_names))}
            data_to_append.append(metrics_dict)

    all_results_df = pd.DataFrame(data_to_append, columns=column_names)
    all_results_df.to_csv('all_metrics_new.csv', index=False)

def main():
    for sub in sorted(os.listdir(BASE_DIR)):
        subdir = os.path.join(BASE_DIR, sub)
        if not os.path.isdir(subdir):
            continue

        # APOC pocket metrics (RMSD, Seq identity, PS-Score)
        apoc_path = os.path.join(subdir, "apoc_output.txt")
        apoc_rmsd, apoc_seq_identity, apoc_ps_score = parse_apoc_pocket_metrics(apoc_path)

        # DockRMSD pocket
        dockrmsd_pocket_path = os.path.join(subdir, "dockrmsd_pocket_output.txt")
        dockrmsd_pocket = parse_dockrmsd(dockrmsd_pocket_path) if os.path.exists(dockrmsd_pocket_path) else ""

        # Write metrics including new APoc fields
        metrics_path = os.path.join(subdir, "metrics.dat")
        with open(metrics_path, "w", encoding="utf-8") as f:
            line = ",".join([sub, apoc_rmsd, apoc_seq_identity, apoc_ps_score, dockrmsd_pocket])
            f.write(line)

    # Combine metrics from all systems
    saveall()

if __name__ == "__main__":
    main()
