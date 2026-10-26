#!/bin/bash
#SBATCH --job-name=CHANGEME
#SBATCH --account=vaquerizas
#SBATCH --partition=gpu            # up to 7 days; 4 nodes x 4 NVIDIA A5000, 56 CPUs and ~500 GB per node
#SBATCH --gres=gpu:A5000:1
#SBATCH --cpus-per-task=12         # fair share of a 56-CPU, 4-GPU node
#SBATCH --mem=96G                  # fair share of ~500 GB per node
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# QOS is set by the partition; no --qos line needed. Avoid --exclusive unless a whole node is genuinely needed.
# Submit from the project root on the cluster (projects/<project>), with logs/ already created.
set -euo pipefail
mkdir -p logs
echo "Started $(date) on $(hostname)"
nvidia-smi -L
echo "deployed commit: $(cat ../../.deployed_commit 2>/dev/null || echo unknown)"
# module load / conda activate here
# For multi-GPU torchrun, derive the process count from what the job was given:
#   NGPUS="${SLURM_GPUS_ON_NODE:-$(nvidia-smi -L | wc -l)}"
#   torchrun --nproc_per_node=${NGPUS} your_script.py
# your command here
