#!/usr/bin/env bash
# pull.sh — rsync results from cluster to local (dry run unless -y)
set -euo pipefail
source "$(dirname "$0")/config.env"

DRY="-n"; [[ "${1:-}" == "-y" ]] && { DRY=""; shift; }
projects=("$@")
[[ ${#projects[@]} -eq 0 ]] && projects=($(ls "$LOCAL_ROOT/projects"))

for p in "${projects[@]}"; do
  mkdir -p "$LOCAL_ROOT/projects/$p/results"
  echo ">>> $p"
  rsync -avh $DRY --max-size=200m \
    --exclude-from="$(dirname "$0")/excludes.txt" \
    "$CLUSTER_HOST:$REMOTE_ROOT/projects/$p/results/" \
    "$LOCAL_ROOT/projects/$p/results/"
done

[[ -n "$DRY" ]] && echo ">> Dry run only — re-run with -y to apply"
