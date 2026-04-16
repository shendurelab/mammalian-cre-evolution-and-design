# Compact CRE optimization

Iterative 1 bp-deletion sequence compaction of mouse CRE tiles using
chromBPNet marginal footprinting. For each CRE tile, `compact_CREs.py`
greedily removes the single base whose deletion **maximizes** predicted
endoderm accessibility, then repeats on the shortened sequence until
nothing is left. The trajectory reveals which bases can be dropped with
minimal activity loss and is used to design compact synthetic CREs.

This is the same algorithm described in the paper; the script here is
self-contained: small inputs ship under `data/`, the CRE tile FASTA is
reused from the neighboring `pCREs_v2/data/` directory, and the two large
external inputs (the mm10 reference genome and the trained chromBPNet
model) are supplied via environment variables.

## Contents

| File | Purpose |
| --- | --- |
| `compact_CREs.py` | Main driver. Loads the chromBPNet endoderm nobias model, runs greedy 1 bp-deletion compaction on each CRE tile, writes per-step marginal footprints to a TSV. |
| `data/fold_0.json` | Chrom split used during chromBPNet training (test = chr1, chr4, chr6). The script samples 1000 **test-chromosome** negative regions as the background for marginal footprinting, matching the training split. |
| `data/endo_negatives_10ksubset.bed` | 10 000-region subset of endoderm negative (non-peak) genomic regions. 1000 are sampled (seed = 42) to build the background sequence pool. |

## Reused inputs (no duplication)

The CRE tile FASTA is **not** copied here — it lives in the sibling pCRE
library directory and is loaded by relative path:

```
scripts/oligo_design/pCREs_v2/data/DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa
```

If that file is moved, update `CRE_FASTA` at the top of `compact_CREs.py`.

## Inputs you must supply

Two files are too large to ship and must be provided via environment
variables (or by editing the defaults at the top of the script):

| Env var | What it points to | Typical size |
| --- | --- | --- |
| `MM10_FASTA` | UCSC mm10 reference genome FASTA (indexed; `.fai` must sit next to it) | ~3 GB |
| `CHROMBPNET_NOBIAS_MODEL` | Trained chromBPNet endoderm **nobias** model (`chrombpnet_nobias.h5`) | ~25 MB |

The `chrombpnet_nobias.h5` model is the sequence-only factorized model
produced by the training pipeline in `scripts/chrombpnet_modeling/` (see
that directory's README for how to train it from scATAC data).

Output location is controlled by `$OUTPUT_TSV` (default
`compact_marginal_fp.txt` in the current working directory).

## Requirements

- Python ≥ 3.8 with a working `chrombpnet` environment. The script pulls
  helper modules directly from the installed package
  (`chrombpnet.training.utils.losses`, `...data_utils.get_seq`,
  `...one_hot`).
- TensorFlow / Keras (via `chrombpnet`'s conda env).
- `pandas`, `numpy`, `biopython`, `pyfaidx`, `pyBigWig`, `deepdish`,
  `matplotlib`.
- A GPU is strongly recommended — marginal footprinting on 1000
  background sequences × hundreds of positions per CRE is compute-heavy.
  `CUDA_VISIBLE_DEVICES=0` is hardcoded; set it to whichever GPU you have.

## Usage

```bash
export MM10_FASTA=/path/to/mm10.fa
export CHROMBPNET_NOBIAS_MODEL=/path/to/endo/chrombpnet_nobias.h5
# optional:
# export OUTPUT_TSV=my_compact_run.txt

python compact_CREs.py
```

Output is a tab-separated table with one row per (CRE, step) including:
- `id` — `ref` for the wild-type, or `del_<pos>` for a deletion at
  1-indexed position `pos`
- `length` — sequence length after the deletion(s)
- `num_del` — cumulative number of deletions (0 for the reference)
- `seq` — current sequence
- `motif_footprint`, `norm_footprint` — predicted marginal footprint and
  its value normalized to the empty-motif control
- `CRE` — CRE tile ID

## How it works

1. **Background model**: loads the chromBPNet endoderm nobias model.
   Samples 1000 negative regions from the test-fold chromosomes
   (`fold_0.json`), one-hot encodes a 2114 bp window around each, and
   computes an "empty motif" (control) marginal footprint. All
   subsequent footprints are normalized to this control to remove
   baseline accessibility drift.
2. **Wild-type pass**: for each CRE tile, predicts the wild-type
   marginal footprint by inserting the CRE sequence at the center of the
   1000 background windows and averaging the model output.
3. **Greedy 1 bp compaction** (inner loop, one step = one base removed):
   - Enumerate all *N* single-base deletions of the current sequence.
   - Score each by marginal footprint (in parallel via an 8-thread
     executor — `concurrent.futures.ThreadPoolExecutor`).
   - Keep the deletion with the **highest** normalized footprint as the
     new current sequence.
   - Record that step.
   - Repeat until the sequence is fully consumed (N steps for an N bp
     CRE).
4. **Output**: concatenate per-CRE trajectories and write TSV.

## Reproducibility notes

- The 1000 background-region subsample is seeded (`random_state=42`), so
  the control footprint is deterministic given the same input BED.
- The deletion search is deterministic — at each step the deletion with
  the strictly highest `norm_footprint` is chosen (ties are broken by
  pandas' stable sort on insertion order, i.e. the earliest position).
- chromBPNet model inference is itself deterministic on a given GPU
  with fixed TensorFlow versions, so repeated runs produce the same
  trajectory.

## Downstream

The resulting compaction trajectories were filtered to pick sequence
lengths that retain ≥ some fraction of predicted wild-type activity,
and those compacted sequences were included in the Twist oligo order
alongside the pCRE v2 library. See
`scripts/oligo_design/README.md` for the full oligo-design context.
