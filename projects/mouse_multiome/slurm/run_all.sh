#!/bin/bash
# Chain steps 01-04 with Slurm dependencies. Run from projects/mouse_multiome on the cluster.
# Before: check inputs, put the GENCODE vM38 GTF in data/reference/, set annotation columns in config.
set -euo pipefail
mkdir -p logs
j1=$(sbatch --parsable slurm/01_QC_filtering.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 slurm/02_peaks_and_normalise.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 slurm/03_integrate_cluster.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j3 slurm/04_annotate.sbatch)
echo "Submitted: $j1 $j2 $j3 $j4"
