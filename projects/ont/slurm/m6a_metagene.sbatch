#!/bin/bash
#SBATCH --job-name=m6a_metagene
#SBATCH --account=vaquerizas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# QOS is set by the partition; no --qos line needed.
# Submit from the project root on the cluster (projects/ont) with logs/ already created:
#   mkdir -p logs && sbatch slurm/m6a_metagene.sbatch
set -euo pipefail

# ---- inputs -----------------------------------------------------------------
BW=/home/mp2924/mnt/storage/vaquerizas/irina/spermatogenesis/ONT_direct_RNA_seq_mouse_testis_bulk/analysis/IGV/M7T1_discover/M7T1.mods.m6A.bw
GTF=/home/mp2924/mnt/storage/vaquerizas/Matt/mouse_multiome/ref_genome/gencode.vM38.annotation.gtf.gz
CONDA_ENV=deeptools

# ---- metagene geometry (nt; must be multiples of BIN) ----------------------
BIN=10
FLANK=500        # upstream of TSS and downstream of TES
UTR_LEN=200      # 5'UTR and 3'UTR each scaled to this length
CDS_LEN=600      # CDS (gene body) scaled to this length

OUT=results/m6a_metagene
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p logs "$OUT"
echo "Started $(date) on $(hostname)"
echo "deployed commit: $(cat ../../.deployed_commit 2>/dev/null || echo unknown)"

# ---- environment ------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then module load conda 2>/dev/null || module load miniconda 2>/dev/null || true; fi
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
# ignore packages in ~/.local (a numpy 2.x there shadows the env's numpy and breaks deepTools: np.NAN)
export PYTHONNOUSERSITE=1
python -c "import numpy; print('numpy', numpy.__version__, numpy.__file__)"
computeMatrix --version

# ---- 1. feature BEDs (5'UTR / CDS / 3'UTR, spliced BED12, one canonical transcript per gene)
python code/make_feature_beds.py --gtf "$GTF" --outdir "$OUT/beds"

# ---- 2. sanity check: chromosome naming agrees between bigWig and GTF-derived BEDs
python - "$BW" "$OUT/beds/cds.bed" <<'PY'
import sys, pyBigWig
bw = pyBigWig.open(sys.argv[1]); chroms = set(bw.chroms())
bed = {l.split("\t", 1)[0] for l in open(sys.argv[2])}
common = chroms & bed
print(f"bigWig chroms: {len(chroms)}, BED chroms: {len(bed)}, shared: {len(common)}")
if not common:
    sys.exit("ERROR: no shared chromosome names (chr1 vs 1?)")
PY

# ---- 3. computeMatrix per feature (--metagene joins exons of each BED12 and skips introns)
COMMON="-S $BW --binSize $BIN -p $THREADS --metagene"
computeMatrix scale-regions $COMMON -R "$OUT/beds/utr5.bed" -b $FLANK -a 0 --regionBodyLength $UTR_LEN \
    -o "$OUT/utr5.mat.gz"
computeMatrix scale-regions $COMMON -R "$OUT/beds/cds.bed"  -b 0 -a 0 --regionBodyLength $CDS_LEN \
    -o "$OUT/cds.mat.gz"
computeMatrix scale-regions $COMMON -R "$OUT/beds/utr3.bed" -b 0 -a $FLANK --regionBodyLength $UTR_LEN \
    -o "$OUT/utr3.mat.gz"

# ---- 4. stitched, annotated plot
python code/plot_metagene.py --utr5 "$OUT/utr5.mat.gz" --cds "$OUT/cds.mat.gz" --utr3 "$OUT/utr3.mat.gz" \
    --out-prefix "$OUT/M7T1_m6A_metagene" --title "M7T1 m6A (mouse testis, ONT dRNA-seq)"

echo "Finished $(date)"
