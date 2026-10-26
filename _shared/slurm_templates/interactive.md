# Interactive sessions (account: vaquerizas)

Partitions on jex (QOS is applied by the partition):

| Partition | Max time | Notes |
|---|---|---|
| nice (default) | 6 h | short jobs; max 64 GB per node; CPU or GPU nodes |
| cpu | 3 days | 7 CPU nodes, 64 CPUs, up to 256 GB per node |
| hmem | 7 days | 1 node, 64 CPUs, ~4 TB RAM |
| gpu | 7 days | 4 nodes x 4 A5000 GPUs |
| int | 16 h | interactive/debugging (qos_int), nodes can be shared |
| ctrl | 30 days | very light long-running controller jobs only (qos_ctrl, max ~8 GB per node) |

## Debugging sessions

CPU:
    srun --account=vaquerizas --partition=int --cpus-per-task=4 --mem=16G --time=04:00:00 --pty bash

GPU:
    srun --account=vaquerizas --partition=int --gres=gpu:A5000:1 --cpus-per-task=8 --mem=48G --time=04:00:00 --pty bash

## Reminders

- Submit batch jobs from the project root on the cluster: cd .../research/projects/<project>
- logs/ must exist there before you submit.
- Run a short smoke test (one epoch or a few minutes on a small input) before any long run, and check the first output/checkpoint is actually written.
