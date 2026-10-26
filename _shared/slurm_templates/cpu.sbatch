#!/bin/bash
#SBATCH --job-name=CHANGEME
#SBATCH --account=vaquerizas
#SBATCH --partition=cpu            # cpu: up to 3 days, 256 GB per node. Use 'hmem' (up to 7 days, ~4 TB) for very large memory.
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# QOS is set by the partition; no --qos line needed.
# Submit from the project root on the cluster (projects/<project>), with logs/ already created.
set -euo pipefail
mkdir -p logs
echo "Started $(date) on $(hostname)"
echo "deployed commit: $(cat ../../.deployed_commit 2>/dev/null || echo unknown)"
# module load / conda activate here
# your command here
