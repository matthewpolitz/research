# research/

This folder is version-controlled with git and deployed to the cluster via `git push`.

## What lives here

| Folder | Tracked? | Notes |
|--------|----------|-------|
| `projects/*/code/` | ✅ git | All analysis scripts |
| `projects/*/slurm/` | ✅ git | Job submission scripts |
| `projects/*/docs/*.md` | ✅ git | Notes in Markdown |
| `_shared/` | ✅ git | Templates, env specs |
| `sync/` | ✅ git | Sync scripts |
| `_admin/` | local only | Drafts, personal notes |
| `projects/*/results/` | ❌ never | Comes back via rsync (pull.sh) |
| `projects/*/data/` | ❌ never | Stays on cluster |
| `_shared/genomes/` | ❌ never | Symlink to shared storage |

## Rules

1. **Edit code locally, deploy with `git push`.** Do not edit tracked files directly on the cluster — the next push overwrites them.
2. **Results flow one way: cluster → local** via `sync/pull.sh`.
3. **Data, checkpoints, logs stay on the cluster.** Never commit large or binary files.
4. **Peer-review material** (manuscripts or grants under review) must not go into any AI tool unless the journal/funder explicitly allows it.

## Key paths

| | Local | Cluster |
|---|---|---|
| Research root | `~/Documents/TE_spermatogenesis_project/research/` | `/home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/research/` |
| Cluster alias | — | `hpc` |
