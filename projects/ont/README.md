# ont

Oxford Nanopore direct RNA-seq analysis of bulk mouse testis (m6A from modification-called bigWigs).

## Layout

- `code/`    — analysis scripts (tracked by git)
- `slurm/`   — Slurm job scripts (tracked by git)
- `results/` — outputs pulled from the cluster with `../../sync/pull.sh` (not tracked)
- `docs/`    — Markdown notes (tracked)
- `data/`    — data on the cluster only (not tracked)

## Analyses

### 1. m6A metagene (`slurm/m6a_metagene.sbatch`)

Mean m6A bigWig signal across transcripts, labelled upstream | TSS | 5'UTR | start codon | CDS (gene body) | stop codon | 3'UTR | TES | downstream.

- Input bigWig: `M7T1.mods.m6A.bw` (path set at the top of the sbatch file)
- Annotation: GENCODE vM38 GTF; one Ensembl_canonical protein-coding transcript per gene, complete CDS, 5'UTR and 3'UTR >= 20 nt
- `code/make_feature_beds.py` builds spliced, strand-aware BED12s per feature; deepTools `computeMatrix scale-regions --metagene` (env `deeptools`) is run per feature; `code/plot_metagene.py` stitches them (500 nt flanks, 10 nt bins)
- Outputs in `results/m6a_metagene/`: `M7T1_m6A_metagene.{png,pdf}`, `M7T1_m6A_metagene_profile.tsv`, matrices, BEDs, `beds/transcripts.tsv`

Run on the cluster:

```bash
ssh hpc
cd /home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/research/projects/ont
mkdir -p logs
sbatch slurm/m6a_metagene.sbatch
```

Then locally: `../../sync/pull.sh` (dry run first, then `-y`).

### Caveats

- Untested on the cluster; the `deeptools` conda activation lines may need adjusting.
- The bigWig is read from Irina's folder via the shared storage mount (read access needed).
- If the bigWig is unstranded, signal from overlapping antisense genes is included.
- Check the bigWig's units (modified fraction vs counts) when labelling the y-axis; edit `--ylabel` in `plot_metagene.py` call if needed.
