# Data for `scripts/make_ext_data_fig1.R`

Inputs consumed by `scripts/make_ext_data_fig1.R`, organized for reproducibility.

## Contents

### Top-level reference files
- `clusters.csv` — cell-type ➜ colour + observed/predicted bigwig paths for Gviz tracks. Bigwig paths are overridden at runtime via the `BIGWIG_DIR` environment variable (see script header).
- `DMS_300bp_tiles.bed` — genomic coordinates of the five oCRE max-activity tiles (Gata4, Epas1, Lama1, Sparc, Bend5), used to highlight regions in Gviz tracks.
- `overlap.peaks.resolved.500.filtered_cmb_peaks.bed` — non-overlapping 500 bp consensus peak set across the four mEB cell types, used to sample background regions when computing bigwig track ymax values.
- `mm10_biomart_GeneExport.txt` — Ensembl/BioMart gene export, used to filter refGene transcripts to protein-coding loci for the Gviz `GeneRegionTrack`.
- `endo_modisco_motifs.html` — TF-MoDISco motif discovery report (endoderm model), parsed for the per-motif `num_seqlets` count displayed in `ext_data_fig4.pdf`.
- `CIS-BP_Mus_musculus.meme` — CIS-BP v2.00 *Mus musculus* motif MEME file, used only to map MoDISco `match*` IDs to human-readable TF names.

### `consensus_peak_preds/`
Per-cell observed vs predicted log-count tables for the consensus peak set. Prepared from the large chromBPNet prediction HDF5 files (`predictions/logcounts` column only) by `scripts/utils/preprocess_ext_data_fig1_h5.R`.

Files: `{endo,pluri,meso,ecto}_cons_peaks_obs_pred.txt.gz`

Columns: `chr, start, end, logcts, pred_logcts, V6, V7, V8, qval, summit`

### `da_peak_preds/`
Same format as above but restricted to per-cell differentially accessible (DA) peak sets.

Files: `{endo,pluri,meso,ecto}_da_peaks_obs_pred.txt.gz`

### `gkm_svm/`
Files for the chromBPNet vs gkm-SVM comparison on a length-matched endoderm holdout (positive = consensus peaks, negative = flanking regions).

- `gkm_length_{pos,neg}_pred.txt` — gkm-SVM prediction per held-out peak.
- `gkm_length_{pos,neg}_test_acc.txt` — observed accessibility (bedtools `bigWigAverageOverBed` output) per peak.
- `gkm_length_pos_test.DA_overlap.bed` — peaks overlapping DA regions.
- `chrombpnet_length_{pos,neg}_test_acc.txt` — observed accessibility for chromBPNet-matched peaks.
- `chrombpnet_length_pos_test.DA_overlap.bed` — DA overlap for chromBPNet set.
- `endo_{pos,neg}_chrombpnet_logcounts.txt.gz` — chromBPNet `predictions/logcounts` extracted from the per-peak HDF5 files.

### `insilico_tf_profiles/`
Per-cell TF injection profiles from the in-silico marginal footprinting experiment (TSV, gzipped). One file per injection cell type (`endo`, `pluri`, `meso`, `ecto`); each contains predicted average accessibility around an injected motif across all four cell-type models, with and without the motif.

Files: `insilico_inj_{cell}TF_profiles_across_cell_models_random_seq.tsv.gz`

## Regenerating the preprocessed files

The per-peak HDF5 files used as source are ~1.5 GB each (four cells × two peak sets) and are therefore not suitable for the repo. To rebuild the `consensus_peak_preds/`, `da_peak_preds/`, and `gkm_svm/endo_{pos,neg}_chrombpnet_logcounts.txt.gz` tables from those HDF5 files, run:

```bash
Rscript scripts/utils/preprocess_ext_data_fig1_h5.R
```

This reads only the `predictions/logcounts` dataset from each HDF5 file and writes the slim TSVs used by the figure script.

## Missing (external) bigwig tracks

The Gviz browser-track panel (`ext_data_fig3_pII.pdf`) requires per-cell bigwigs not included in the repo:

```
endo.bw, meso.bw, ecto.bw, pluri.bw                                   (observed)
endo_chrombpnet_nobias.bw, meso_chrombpnet_nobias.bw, …               (predicted)
```

Set `BIGWIG_DIR` to a directory containing these files before running the figure script:

```bash
BIGWIG_DIR=/path/to/bigwigs Rscript scripts/make_ext_data_fig1.R
```

Without `BIGWIG_DIR` (or with missing files), the figure script skips that panel automatically with a warning.
