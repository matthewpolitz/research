#!/bin/bash
#SBATCH --job-name=mm_02_peaks_and_normalise
#SBATCH --account=vaquerizas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Submit from projects/mouse_multiome on the cluster (logs/ must exist): sbatch slurm/02_peaks_and_normalise.sbatch [sample]
set -euo pipefail
mkdir -p logs
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mouse_multiome
export MM_COMMIT="$(cat ../../.deployed_commit 2>/dev/null || echo unknown)"
echo "Started $(date) on $(hostname); deployed commit: $MM_COMMIT"
Rscript code/02_peaks_and_normalise.R "$@"
