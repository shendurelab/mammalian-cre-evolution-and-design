# ChromBPNet Modeling

Train cell-type-specific [chromBPNet](https://github.com/kundajelab/chrombpnet) models on mouse embryoid-body scATAC-seq. Output models predict base-resolution accessibility from sequence and are consumed by `scripts/oligo_design/`.

**Cell types**: `meso`, `endo`, `pluri`, `ecto`. Input: 2114 bp sequence → 1000 bp profile + total count. Fold 0 splits: test = chr1,4,6; val = chr10,17; train = rest.

## Run order

```bash
# 1) Preprocessing: merge BAMs, call peaks, make background sets
./sm_script.sh                    # wraps snakemake, submits SGE

# 2) Bias models (GPU required)
qsub train_gpu.sh

# 3) Factorized nobias models (GPU, after step 2)
qsub train_factorized_gpu.sh

# 4) Downstream (predictions, contribs, MoDISco, footprints):
#    uncomment the relevant rule in Snakefile's `rule all`, then ./sm_script.sh
```

## Scripts

| File | What |
|---|---|
| `Snakefile` | Snakemake workflow: BAM subsetting by cell barcode, merging, MACS2 peak calling, consensus-peak resolution, blacklist cleaning, background sampling, and downstream (genome-wide predictions, TSS predictions, DeepSHAP contrib scores, MoDISco motifs, marginal TF footprinting). |
| `sm_script.sh` | Snakemake → SGE qsub wrapper. |
| `train_gpu.sh` | SGE GPU job: `chrombpnet bias pipeline` for each cell type → `bias_models/{cell}/models/{cell}_bias.h5`. |
| `train_factorized_gpu.sh` | SGE GPU job: `chrombpnet pipeline` (needs bias model) → `chrombpnet_models/{cell}/models/chrombpnet_nobias.h5` (sequence-only model used by oligo design). |
| `resolve.py` | Peak overlap resolution (Kundaje lab): sorts by significance, greedy non-overlapping 500 bp windows at summits. |
| `realign_footprints.py` | Per-motif footprint PNGs from the `chrombpnet footprints` HDF5. |

## Dependencies

- `chrombpnet` (conda env), `samtools`, `bedtools`, `macs2`, `subset-bam` (10x), Kent utilities (`bigWigAverageOverBed`, `bedsort`), snakemake.
- Inputs expected at training time: per-replicate scATAC BAMs, a cell-barcode→cluster table, ArchR cluster bigwigs, `EBmotif_to_pwm.TF.tsv`, mm10 FASTA + chrom.sizes, mm10 blacklist v2.
