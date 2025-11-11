#!/usr/bin/env bash
set -euo pipefail

JOBS=1
DRYRUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="${2:-$PWD}"; shift 2 ;;
    --jobs|-j) JOBS="${2:-1}"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

CONVERT_PY="convert.py"
RUNR_SH="runr.py"

for f in "$CONVERT_PY" "$RUNR_SH"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing required file: $f" >&2; exit 2; }
done

is_highest_better_metric() {
  case "$1" in
    iptm_complex_af3|ptm_ligand_only_af3) return 0 ;; 
    *) return 1 ;;
  esac
}

process_one_dir() {
  local recp_dir="$1"
  local metric_dir="$(dirname "$recp_dir")"
  local metric="$(basename "$metric_dir")"
  local split_csv="$recp_dir/split.csv"
  local ligs="$recp_dir/ligands.name"
  local decs="$recp_dir/decoys.name"

  [[ -f "$split_csv" ]] || { echo "[SKIP] no split.csv: $recp_dir"; return 0; }
  [[ -s "$ligs" && -s "$decs" ]] || { echo "[SKIP] missing ligands/decoys: $recp_dir"; return 0; }

  echo "[INFO] >>> $metric / $(basename "$recp_dir")"

  if [[ "$DRYRUN" -eq 1 ]]; then
    if is_highest_better_metric "$metric"; then
      echo "  Would rank descending: (cd $recp_dir && conda activate dock37_py27 && conda activate py3d && python $RUNR_SH split.csv && conda deactivate)"
    else
      echo "  Would rank ascending:  (cd $recp_dir && conda activate dock37_py27 && conda activate py3d && python $CONVERT_PY split.csv && conda deactivate)"
    fi
    echo "  Would AUC/logAUC: (cd $recp_dir && python ./logAUC_calculation/dockbase_files -i . -l ./ligands.name -d ./decoys.name)"
    echo "                     (cd $recp_dir && python ./logAUC_calculation/plots.py  -i . -l ./ligands.name -d ./decoys.name)"
    return 0
  fi

  if is_highest_better_metric "$metric"; then
    ( cd "$recp_dir" && conda activate dock37_py27 && conda activate py3d && python "$RUNR_SH" "split.csv" && conda deactivate )
  else
    ( cd "$recp_dir" && conda activate dock37_py27 && conda activate py3d && python "$CONVERT_PY" "split.csv" && conda deactivate )
  fi

  ( cd "$recp_dir" && conda activate dock37_py27 && python ./logAUC_calculation/dockbase_files -i . -l ./ligands.name -d ./decoys.name )
  ( cd "$recp_dir" && conda activate dock37_py27 && python ./logAUC_calculation/plots.py  -i . -l ./ligands.name -d ./decoys.name )
}

export CONVERT_PY RUNR_SH DRYRUN
export -f process_one_dir
export -f is_highest_better_metric

mapfile -t TARGET_DIRS < <(find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "split.csv" -printf '%h\n' | sort -u)
echo "[INFO] Found ${#TARGET_DIRS[@]} target directories under $ROOT"

if [[ "$DRYRUN" -eq 1 ]]; then
  for d in "${TARGET_DIRS[@]}"; do process_one_dir "$d"; done
  exit 0
fi

if [[ "$JOBS" -le 1 ]]; then
  for d in "${TARGET_DIRS[@]}"; do process_one_dir "$d"; done
else
  printf "%s\0" "${TARGET_DIRS[@]}" | xargs -0 -n1 -P "$JOBS" -I{} bash -c 'process_one_dir "$@"' _ {}
fi

echo "[DONE] All done."
