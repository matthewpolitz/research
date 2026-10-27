# Setup notes (written 2026-10-26)

Owner: Matthew Politz (mp2924), PhD student, Vaquerizas lab. Git identity: Matthew Politz <mp2924@ic.ac.uk>.

## What was created

- Local tree: `research/` with `_admin/`, `_shared/{genomes,envs,slurm_templates}`, `projects/te_spermatogenesis/{code,slurm,results,docs,data}`, `sync/`.
- `.gitignore`, `README.md` (rules), project README, sync scripts (`pull.sh`, `push.sh` fallback, `excludes.txt`, `config.env`).
- Slurm templates: `_shared/slurm_templates/cpu.sbatch`, `gpu.sbatch`, `interactive.md`.
- Existing folders in `TE_spermatogenesis_project/` (MSCI, alignment, ...) are untouched and sit outside `research/`.

## Values chosen

| Item | Value |
|---|---|
| SSH alias | `hpc` (jex.lms.mrc.ac.uk, login01), key `~/.ssh/id_ed25519_hpc` |
| Local root | `~/Documents/TE_spermatogenesis_project/research` |
| Remote root | `/home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/research` |
| Bare repo | `/home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/git/research.git` (remote `hpc`, branch `main`) |
| Shared genomes | `/opt/resources/genomes` (symlinked at `_shared/genomes` on the cluster) |
| Slurm account | `vaquerizas` |
| Partitions | nice (default, 6 h), cpu (3 d), hmem (7 d), gpu (7 d, A5000), int (16 h interactive), ctrl (30 d, light jobs) |

## Verified by pasted output

SSH login, cluster research root and genomes symlink, bare repo and hook, first push (deployed d993e62), Slurm account/partition info, CPU smoke test (job 1331438), GPU smoke test (one NVIDIA RTX A5000 listed, deployed commit recorded in the log; GPU nodes were temporarily unavailable earlier, which delayed the job).

## Still to do yourself

1. Commit and push the Slurm templates (`git add . && git commit && git push`).
2. Dry-run `sync/pull.sh` (no `-y`) once a project has results on the cluster; then run with `-y`.
3. Smoke tests done for CPU and GPU via `sbatch --wrap`; the full `cpu.sbatch`/`gpu.sbatch` files themselves have not yet been submitted, so run one short real job with each before a long run.
4. Before using a fasta from `/opt/resources/genomes`, check a `.fai` exists; if not, build it once by hand before multi-worker jobs.
5. Check chromosome naming (chr1 vs 1) agrees between fasta, GTF and any cooler/bigwig inputs.
6. Peer-review material (manuscripts, grants under review) must not go into any AI tool unless the journal/funder allows it.

## Where to ask for help

Ask lab members or the LMS HPC/IT team; keep this file updated when paths change.
