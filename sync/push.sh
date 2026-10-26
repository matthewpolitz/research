#!/usr/bin/env bash
# push.sh — FALLBACK ONLY: rsync code/slurm to cluster
# Normal workflow: git push hpc main
# Use this only if you prefer not to use git.
set -euo pipefail
source "$(dirname "$0")/config.env"

DRY="-n"; [[ "${1:-}" == "-y" ]] && { DRY=""; shift; }
projects=("$@")
[[ ${#projects[@]} -eq 0 ]] && projects=($(ls "$LOCAL_ROOT/projects"))

for p in "${projects[@]}"; do
  for d in code slurm; do
    if [[ -n "$DRY" ]]; then
      rsync -avh -n --exclude-from="$(dirname "$0")/excludes.txt" \
        "$LOCAL_ROOT/projects/$p/$d/" \
        "$CLUSTER_HOST:$REMOTE_ROOT/projects/$p/$d/"
    else
      ssh "$CLUSTER_HOST" "mkdir -p '$REMOTE_ROOT/projects/$p/$d'"
      rsync -avh --exclude-from="$(dirname "$0")/excludes.txt" \
        "$LOCAL_ROOT/projects/$p/$d/" \
        "$CLUSTER_HOST:$REMOTE_ROOT/projects/$p/$d/"
    fi
  done
done

# sync shared templates/envs
if [[ -n "$DRY" ]]; then
  rsync -avh -n --exclude-from="$(dirname "$0")/excludes.txt" \
    "$LOCAL_ROOT/_shared/" "$CLUSTER_HOST:$REMOTE_ROOT/_shared/"
else
  ssh "$CLUSTER_HOST" "mkdir -p '$REMOTE_ROOT/_shared'"
  rsync -avh --exclude-from="$(dirname "$0")/excludes.txt" \
    "$LOCAL_ROOT/_shared/" "$CLUSTER_HOST:$REMOTE_ROOT/_shared/"
fi

[[ -n "$DRY" ]] && echo ">> Dry run only — re-run with -y to apply"
