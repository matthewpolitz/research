# te_spermatogenesis

Transposable elements in human spermatogenesis.

## Layout

- `code/`    — analysis scripts (tracked by git)
- `slurm/`   — Slurm job scripts (tracked by git)
- `results/` — output files pulled from cluster with `../../../sync/pull.sh` (not tracked)
- `docs/`    — Markdown notes (tracked); PDFs and exports stay local
- `data/`    — raw/processed data on the cluster only (not tracked)

## Submitting jobs

Always submit from the project root on the cluster. Make sure `logs/` exists there first:

```bash
ssh hpc
cd /home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/research/projects/te_spermatogenesis
mkdir -p logs
sbatch slurm/cpu.sbatch
```
