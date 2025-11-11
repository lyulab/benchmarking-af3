#!/usr/bin/env python3
"""
Split *_running_sum.csv files by metric and recp_name.

For each input file named <metric>_running_sum.csv, this script creates:
  <out_dir>/<metric>/<recp_name>/split.csv
  <out_dir>/<metric>/<recp_name>/ligands.name      (compound_id where is_active == 1)
  <out_dir>/<metric>/<recp_name>/decoys.name       (compound_id where is_active == 0)

Designed to be safe for large files via chunked processing.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd

REQUIRED_COLUMNS = {"recp_name", "compound_id", "is_active"}

def parse_args():
    p = argparse.ArgumentParser(description="Split *_running_sum.csv files by metric and recp_name.")
    p.add_argument(
        "-i", "--input-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory to search for *_running_sum.csv files (default: current dir)"
    )
    p.add_argument(
        "-o", "--out-dir",
        type=Path,
        default=Path.cwd(),
        help="Base output directory (default: current dir)"
    )
    p.add_argument(
        "-c", "--chunksize",
        type=int,
        default=200_000,
        help="Rows per chunk when reading CSV (default: 200000)"
    )
    p.add_argument(
        "--glob",
        type=str,
        default="*_running_sum.csv",
        help="Glob pattern to find input files (default: '*_running_sum.csv')"
    )
    p.add_argument(
        "--sep",
        type=str,
        default=",",
        help="CSV separator (default: ',')"
    )
    p.add_argument(
        "--encoding",
        type=str,
        default="utf-8",
        help="CSV encoding (default: utf-8)"
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="List actions without writing files"
    )
    return p.parse_args()


def metric_from_filename(path: Path) -> str:
    name = path.name
    # Expect "<metric>_running_sum.csv"
    if name.endswith("_running_sum.csv"):
        return name[:-len("_running_sum.csv")]
    # Fallback: drop suffix and return stem
    return path.stem


def ensure_required_columns(df: pd.DataFrame, src: Path):
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{src} is missing required columns: {sorted(missing)}")


def append_dataframe(out_csv: Path, df: pd.DataFrame, header_written_cache: set, dry: bool):
    """
    Append DataFrame to CSV, writing header only once per file path.
    """
    if dry:
        return
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    write_header = out_csv not in header_written_cache and not out_csv.exists()
    df.to_csv(out_csv, mode="a", index=False, header=write_header)
    header_written_cache.add(out_csv)


def append_lines(out_path: Path, lines, dry: bool):
    if dry:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("a", encoding="utf-8") as f:
        for ln in lines:
            f.write(str(ln) + "\n")


def process_file(csv_path: Path, out_base: Path, chunksize: int, sep: str, encoding: str, dry_run: bool) -> None:
    metric = metric_from_filename(csv_path)
    metric_root = out_base / metric

    print(f"[INFO] Processing: {csv_path}  -> metric='{metric}'")

    header_written_cache = set()  # track which split.csv have header written

    chunk_iter = pd.read_csv(
        csv_path,
        sep=sep,
        encoding=encoding,
        chunksize=chunksize,
        low_memory=False
    )

    total_rows = 0
    for chunk_idx, chunk in enumerate(chunk_iter, start=1):
        total_rows += len(chunk)

        ensure_required_columns(chunk, csv_path)

        # Drop obviously empty recp_name entries to avoid creating "nan" folders
        chunk = chunk[chunk["recp_name"].notna()].copy()

        # Iterate groups by recp_name inside this chunk
        for recp_name, sub in chunk.groupby("recp_name", dropna=True, sort=False):
            recp_dir = metric_root / str(recp_name)
            split_csv = recp_dir / "split.csv"
            # Append this group's rows to split.csv (header once globally per file)
            append_dataframe(split_csv, sub, header_written_cache, dry_run)

            # Ligands and decoys
            lig_ids = sub.loc[sub["is_active"] == 1, "compound_id"].dropna().drop_duplicates().tolist()
            dec_ids = sub.loc[sub["is_active"] == 0, "compound_id"].dropna().drop_duplicates().tolist()

            if lig_ids:
                append_lines(recp_dir / "ligands.name", lig_ids, dry_run)
            if dec_ids:
                append_lines(recp_dir / "decoys.name", dec_ids, dry_run)

        print(f"  - chunk {chunk_idx}: {len(chunk)} rows -> groups: {chunk['recp_name'].nunique()}")

    print(f"[DONE] {csv_path.name}: {total_rows} rows processed.")


def main():
    args = parse_args()

    input_dir: Path = args.input_dir
    out_dir: Path = args.out_dir
    pattern: str = args.glob
    chunksize: int = args.chunksize

    if not input_dir.exists():
        print(f"[ERROR] Input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(2)

    files = sorted(input_dir.glob(pattern))
    if not files:
        print(f"[WARN] No files matched pattern '{pattern}' under {input_dir}")
        sys.exit(0)

    print(f"[INFO] Found {len(files)} files to process under {input_dir} with pattern '{pattern}'")
    print(f"[INFO] Output base: {out_dir}")
    if args.dry_run:
        print("[INFO] DRY RUN (no files will be written)")

    for f in files:
        try:
            process_file(f, out_dir, chunksize, args.sep, args.encoding, args.dry_run)
        except Exception as e:
            print(f"[ERROR] Failed on {f}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
