#!/usr/bin/env python3
import re, csv, sys
from pathlib import Path
import pandas as pd

"""
Outputs:
  - auc_summary.csv (long form)
  - pivot_auc.csv, pivot_log_auc.csv, pivot_n_ligands.csv, pivot_n_decoys.csv
    (rows=recp_name, columns=metric)
"""

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/mnt/data")

NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
ROC_HEADER_PAT = re.compile(
    rf'\bAUC\b[^0-9\-+]*({NUM}).*?\blogAUC\b[^0-9\-+]*({NUM})',
    re.IGNORECASE | re.DOTALL
)
AUC_FALLBACK_PATS = [
    re.compile(rf'\bAUC\b[^0-9\-+]*({NUM})', re.IGNORECASE),
    re.compile(rf'\bROC\b[^0-9\-+]*({NUM})', re.IGNORECASE),
]
LOGAUC_FALLBACK_PATS = [
    re.compile(rf'\blogAUC\b[^0-9\-+]*({NUM})', re.IGNORECASE),
    re.compile(rf'\blog\s*AUC\b[^0-9\-+]*({NUM})', re.IGNORECASE),
]
FALLBACK_FILES = ("enrich.out", "enrich.log", "enrich.txt", "roc.txt", "plots.out", "summary.txt")

def parse_roc_own_header(p: Path):
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            text = ''.join([f.readline() for _ in range(5)])
        m = ROC_HEADER_PAT.search(text)
        if m:
            return float(m.group(1)), float(m.group(2))
    except Exception:
        pass
    return None, None

def parse_fallback(text: str):
    auc = None; logauc = None
    for pat in AUC_FALLBACK_PATS:
        m = pat.search(text)
        if m:
            try:
                auc = float(m.group(1)); break
            except ValueError:
                pass
    for pat in LOGAUC_FALLBACK_PATS:
        m = pat.search(text)
        if m:
            try:
                logauc = float(m.group(1)); break
            except ValueError:
                pass
    return auc, logauc

def count_lines(p: Path) -> int:
    try:
        with p.open("r", encoding="utf-8") as f:
            return sum(1 for _ in f)
    except Exception:
        return 0

def collect_rows(root: Path):
    rows = []
    for split_csv in root.rglob("split.csv"):
        recp_dir = split_csv.parent
        metric_dir = recp_dir.parent
        metric = metric_dir.name
        recp_name = recp_dir.name

        ligs = recp_dir / "ligands.name"
        decs = recp_dir / "decoys.name"
        n_ligs = count_lines(ligs) if ligs.exists() else 0
        n_decs = count_lines(decs) if decs.exists() else 0

        auc = None; logauc = None; source_file = ""

        roc_own = recp_dir / "roc_own.txt"
        if roc_own.exists():
            a, l = parse_roc_own_header(roc_own)
            if a is not None or l is not None:
                auc, logauc, source_file = a, l, "roc_own.txt"

        if source_file == "":
            for fname in FALLBACK_FILES:
                fpath = recp_dir / fname
                if not fpath.exists(): continue
                try:
                    text = fpath.read_text(encoding="utf-8", errors="ignore")
                except Exception:
                    continue
                a, l = parse_fallback(text)
                if a is not None or l is not None:
                    auc = a if a is not None else auc
                    logauc = l if l is not None else logauc
                    source_file = fname
                    break

        rows.append({
            "metric": metric,
            "recp_name": recp_name,
            "auc": auc if auc is not None else "",
            "log_auc": logauc if logauc is not None else "",
            "n_ligands": n_ligs,
            "n_decoys": n_decs,
            "source_file": source_file,
            "recp_path": str(recp_dir),
        })
    return rows

def write_summary_csv(root: Path, rows):
    out_csv = root / "auc_summary.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["metric","recp_name","auc","log_auc","n_ligands","n_decoys","source_file","recp_path"])
        w.writeheader()
        w.writerows(rows)
    print(f"[INFO] Wrote {out_csv} with {len(rows)} rows.")
    return out_csv

def make_pivots(root: Path, summary_csv: Path):
    df = pd.read_csv(summary_csv, dtype={"metric":str, "recp_name":str}, keep_default_na=False)

    for col in ["auc", "log_auc"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["n_ligands", "n_decoys"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    agg = (df
           .groupby(["metric","recp_name"], as_index=False)
           .agg({
               "auc":"mean",
               "log_auc":"mean",
               "n_ligands":"max",
               "n_decoys":"max"
           }))

    pivots = {
        "pivot_auc.csv":       agg.pivot(index="recp_name", columns="metric", values="auc"),
        "pivot_log_auc.csv":   agg.pivot(index="recp_name", columns="metric", values="log_auc"),
        "pivot_n_ligands.csv": agg.pivot(index="recp_name", columns="metric", values="n_ligands"),
        "pivot_n_decoys.csv":  agg.pivot(index="recp_name", columns="metric", values="n_decoys"),
    }

    for name, pvt in pivots.items():
        pvt = pvt.sort_index(axis=0).sort_index(axis=1)
        outp = root / name
        pvt.to_csv(outp, float_format="%.6g")
        print(f"[INFO] Wrote {outp} (rows=recp_name, cols=metric)")

def main():
    rows = collect_rows(ROOT)
    summary_csv = write_summary_csv(ROOT, rows)
    make_pivots(ROOT, summary_csv)

if __name__ == "__main__":
    main()
