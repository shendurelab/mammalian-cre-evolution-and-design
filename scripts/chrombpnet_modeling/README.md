# ChromBPNet Modeling

Scripts for training cell-type-specific chromBPNet models from scATAC-seq data. These models predict base-resolution chromatin accessibility profiles from DNA sequence and are used by the oligo design scripts (`scripts/oligo_design/`) for in silico CRE optimization.

## Overview

ChromBPNet (chromatin BPNet) is a deep learning model that learns the relationship between DNA sequence and chromatin accessibility (ATAC-seq signal). The pipeline trains one model per cell type from mouse embryoid body (mEB) scATAC-seq data.

**Cell types modeled**: meso (mesoderm), endo (endoderm), pluri (pluripotent), ecto (ectoderm)

## Pipeline

Training proceeds in three stages:

```
scATAC BAMs (per cell type)
        |
        v
  Snakefile: BAM merging + peak calling + background set generation
        |
        v
  train_gpu.sh: Train Tn5 bias model (per cell type)
        |
        v
  train_factorized_gpu.sh: Train factorized chromBPNet model (per cell type)
        |
        v
  Snakefile: Predictions, contribution scores, TF footprints, motif discovery
```

## Scripts

### `Snakefile`

Snakemake workflow for preprocessing and downstream analysis. Run via `sm_script.sh`.

**Preprocessing rules:**
1. `generate_coarse_replicate_bams` - Subsets scATAC BAMs by cell barcode per replicate
2. `generate_merge_coarse_bams` - Merges replicate BAMs into one per cell type
3. `generate_mouse_blacklist` - Extends mm10 blacklist regions, generates chromosome splits for train/test/validation
4. `generate_coarse_peak_calling` - MACS2 peak calling (shift -100, extsize 200, p < 0.01)
5. `generate_overlaps_peaks` - Merges peaks across cell types into non-overlapping consensus set (500bp)
6. `generate_clean_peaks` - Removes blacklist and tested MPRA peaks, generates negative (background) regions

**Downstream rules:**
7. `generate_chrombpnet_predict_all` - Genome-wide predictions at peaks + background
8. `generate_chrombpnet_predictions` - Predictions at held-out gene body windows
9. `generate_chrombpnet_TSSpredictions` - Predictions at TSS regions
10. `score_heldout` - Score held-out predictions with bigWigAverageOverBed
11. `generate_chrombpnet_contrib_score` - DeepSHAP contribution scores at tested MPRA regions
12. `generate_denovo_motif_discovery` - TF-MoDISco de novo motif discovery
13. `generate_tf_footprinting` - Marginal TF footprinting with known motif PWMs

**Dependencies:**
- `chrombpnet` (conda environment)
- `samtools`, `bedtools`, `macs2`
- `subset-bam` (10x Genomics)
- Kent utilities (`bigWigAverageOverBed`, `bedsort`)
- `resolve.py` (peak resolution)
- `realign_footprints.py` (footprint visualization)

**Input data:**
- scATAC BAMs: `EB_dev_project/data/scATAC_bams/10x_bams/EBD21_crispri_{replicate}.bam`
- Cell barcodes: `cellBC_to_celltype_table_scATAC_mEB_20230201.txt`
- ArchR bigwigs: `scATAC_bigwigs/{cluster}-TileSize-50-normMethod-ReadsInTSS-ArchR.bw`
- Motif PWMs: `EBmotif_to_pwm.TF.tsv`
- mm10 genome: `/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa`
- mm10 blacklist: `mm10-blacklist.v2.bed`

---

### `train_gpu.sh`

SGE GPU job script that trains the **Tn5 bias model** for each cell type.

```bash
#$ -l mfree=100G,h_rt=168:00:00
#$ -l gpgpu=1,cuda=1
```

**What it does:**
For each cell type (meso, endo, pluri, ecto):
1. Runs `chrombpnet bias pipeline` with:
   - Input BAM (merged, sorted)
   - ATAC-seq assay type
   - Cleaned peaks (blacklist-removed)
   - Background (negative) regions
   - Fold 0 train/test split
   - Bias fraction = 0.5
2. Outputs bias model to `bias_models/{cell}/models/{cell}_bias.h5`

**Must be run before** `train_factorized_gpu.sh`.

---

### `train_factorized_gpu.sh`

SGE GPU job script that trains the **factorized chromBPNet model** for each cell type.

```bash
#$ -l mfree=100G,h_rt=168:00:00
#$ -l gpgpu=1,cuda=1
```

**What it does:**
For each cell type (meso, endo, pluri, ecto):
1. Runs `chrombpnet pipeline` with:
   - Same inputs as bias training
   - Pre-trained bias model from `train_gpu.sh`
2. Outputs:
   - `chrombpnet_models/{cell}/models/chrombpnet.h5` (full model with bias)
   - `chrombpnet_models/{cell}/models/chrombpnet_nobias.h5` (sequence-only model, used by oligo design scripts)

**Must be run after** `train_gpu.sh`.

---

### `sm_script.sh`

Snakemake submission wrapper for SGE cluster:
```bash
snakemake --cluster "qsub -l mfree={params.memory}G -l h_rt={params.run_time} ..." -j 500
```

---

### `resolve.py`

Resolves overlapping narrowPeak entries into a non-overlapping set. Takes peaks sorted by significance and greedily selects non-overlapping 500bp windows centered on summits. From the Kundaje lab.

---

### `realign_footprints.py`

Visualizes motif footprints from HDF5 files. Generates per-motif PNG plots of the predicted Tn5 insertion profile around the motif center.

## Training order

```bash
# 1. Preprocessing (run from EB_chrombpnet/ directory)
./sm_script.sh   # or: snakemake -np to dry-run first

# 2. Train bias models (GPU required)
qsub train_gpu.sh

# 3. Train factorized models (GPU required, after bias training completes)
qsub train_factorized_gpu.sh

# 4. Downstream analysis (predictions, footprints, motifs)
# Uncomment desired rules in Snakefile 'rule all' input, then:
./sm_script.sh
```

## Model architecture

- **Input**: 2114bp DNA sequence (one-hot encoded)
- **Output**: 1000bp base-resolution ATAC-seq profile (multinomial) + total count (scalar)
- **Architecture**: Factorized = sequence model + bias model. The `chrombpnet_nobias.h5` model captures sequence-intrinsic accessibility, removing Tn5 insertion bias.

## Chromosome splits (fold_0)

- **Training**: all chromosomes except test and validation
- **Test**: chr1, chr4, chr6
- **Validation**: chr10, chr17
