# mouse_multiome

10x Genomics Multiome (snRNA + snATAC) of adult mouse testis, two samples (Mouse1, Mouse2). cellranger-arc already run.
Pipeline is R: Seurat v5 + Signac (+ Harmony, scDblFinder, MACS2 via Signac::CallPeaks).

## Inputs (cluster only)

- cellranger-arc: `/home/mp2924/mnt/storage/vaquerizas/Matt/mouse_multiome/cellranger-arc/{Mouse1,Mouse2}`
- Markers: `/home/mp2924/mnt/storage/vaquerizas/Matt/mouse_multiome/pre_processing/Ernst_mouse_cluster_marker_genes.xlsx`
  (Ernst et al., "Staged developmental mapping and X chromosome transcriptional dynamics during mouse spermatogenesis")

## Steps

| Step | Script | What it does | Output |
|---|---|---|---|
| 00 | `00_check_inputs.R`, `00_inspect_markers.R` | verify cellranger-arc files; show marker xlsx structure | console |
| 01 | `01_QC_filtering.Rmd` (run via `run_01_QC_filtering.R`) | per-sample QC from your QC Rmd: sequential TSS / nucleosome signal / FRiP / scDblFinder (RNA+ATAC) / MAD counts+features / %mt filters, then SoupX | `data/qc/<s>/<s>_qc_filtered.rds`, html report, tables, `data/qc/summary/` |
| 02 | `02_peaks_and_normalise.R` | MACS2 peaks per sample, merged peak set, requantify, SCT/log-norm, TF-IDF + SVD | `data/<s>.norm.rds` |
| 03 | `03_integrate_cluster.R` | merge, Harmony (RNA PCA + ATAC LSI), WNN, Leiden clustering at several resolutions | `data/mouse_multiome.integrated.rds` |
| 04 | `04_annotate.R` | module-score clusters against Ernst et al. markers + canonical marker check | `data/mouse_multiome.annotated.rds` |
| 05 | `05_summary_report.R` | composition tables/plots | `results/05_summary/` |

QC thresholds, doublet rule and SoupX settings are the `params:` block at the top of `code/01_QC_filtering.Rmd`
(copied from `Documents/mouse_multiome/scripts/1_QC_filtering_claude.Rmd`; only default paths changed).
Everything else is in `config/config.yml`. Samples/paths for step 01: `config/samples_QC.csv`.
Genome: GRCm39 / GENCODE vM38 (as in the QC Rmd).

## Running (cluster)

```bash
ssh hpc
cd /home/mp2924/mnt/storage/vaquerizas/Matt/TE_spermatogenesis_project/research/projects/mouse_multiome
mkdir -p logs
conda env create -f environment.yml          # once
conda activate mouse_multiome
Rscript code/00_check_inputs.R               # confirm paths
Rscript code/00_inspect_markers.R            # then set annotation$cluster_col / gene_col in config
# QC step needs the GENCODE vM38 GTF (and optionally the cached Signac annotation rds) in data/reference/
mkdir -p data/reference   # copy/symlink gencode.vM38.annotation.gtf.gz (+ gencode.vM38_signac_annotation.rds) here
bash slurm/run_all.sh                        # steps 01-04 chained (01 = one array task per sample)
```

Later steps can be re-run alone after editing the config (changing QC params means re-running 01 onwards).
After step 04, review `results/04_annotation/`, optionally write `docs/manual_annotation.csv` (`cluster,cell_type`) and re-run step 04.
Pull results with `sync/pull.sh -y mouse_multiome`.

## Status / caveats

Written but NOT yet run on the cluster. Untested assumptions: cellranger-arc `outs/` layout, marker xlsx column names, the cellranger-arc reference matching GRCm39 / GENCODE vM38 (your local QC data used `.../outs` directly under each run; here I assume `cellranger-arc/<sample>/outs`).
