# Analysis notes: mouse_multiome

## Design choices
- Per-sample QC (Rmd from Documents/mouse_multiome) then merge; peaks re-called per sample with MACS2 and merged so both samples share one feature set.
- Harmony on RNA PCA and ATAC LSI separately, then WNN. Two samples only, so batch effect is expected to be small; compare `integration_before_after.png`.
- Germ cell QC: late spermatids have low RNA/transcription and condensed chromatin (low TSS enrichment, altered nucleosome signal). Check whether filters remove whole spermatid clusters (compare `results/02_qc` flags vs UMAP).
- Ambient RNA is typically substantial in testis nuclei; consider SoupX/decontX if somatic markers appear in germ clusters.

## TODO
- [ ] confirm sample sheet paths (config/samples_QC.csv) and that data/reference has the GTF
- [ ] set annotation columns after `00_inspect_markers.R`
- [ ] smoke test on one sample
