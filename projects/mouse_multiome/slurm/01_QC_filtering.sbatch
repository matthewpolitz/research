#!/bin/bash
#SBATCH --job-name=mm_01_QC
#SBATCH --account=vaquerizas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# One array task per sample (order = config/samples_QC.csv). Submit from projects/mouse_multiome: sbatch slurm/01_QC_filtering.sbatch
set -euo pipefail
mkdir -p logs
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate mouse_multiome
export MM_COMMIT="$(cat ../../.deployed_commit 2>/dev/null || echo unknown)"
SAMPLE=$(tail -n +2 config/samples_QC.csv | cut -d, -f1 | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
echo "Started $(date) on $(hostname); sample $SAMPLE; deployed commit: $MM_COMMIT"
Rscript code/run_01_QC_filtering.R "$SAMPLE"
