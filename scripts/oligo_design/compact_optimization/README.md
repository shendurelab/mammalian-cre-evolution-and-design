# Compact CRE optimization

`compact_CREs.py` — iterative greedy 1 bp-deletion compaction of mouse CRE tiles using chromBPNet marginal footprinting. At each step it evaluates every single-base deletion and keeps the one with the highest predicted endoderm accessibility; repeats until the sequence is empty. The per-step trajectory reveals which bases can be dropped with minimal activity loss and is used to design compact synthetic CREs.

## Inputs

Small inputs ship here:

- `data/fold_0.json` — chromBPNet train/test split (test = chr1, chr4, chr6). 1000 test-chromosome negatives are sampled for the background pool.
- `data/endo_negatives_10ksubset.bed` — 10 k endoderm non-peak regions; 1000 are sampled (seed 42).

The CRE tile FASTA is reused from the neighboring derivatization library (no duplication): `../derivatization_CRE/data/DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa`.

Two external inputs come from env vars (or edit the defaults at the top of the script):

| Env var | What |
|---|---|
| `MM10_FASTA` | Indexed mm10 FASTA (`.fai` next to it). ~3 GB. |
| `CHROMBPNET_NOBIAS_MODEL` | Endoderm nobias model (`chrombpnet_nobias.h5`, ~25 MB). Defaults to the repo copy under `data/chrombpnet_models/endo/`. |

## Usage

```bash
export MM10_FASTA=/path/to/mm10.fa
# CHROMBPNET_NOBIAS_MODEL defaults to the repo-local model; override if needed.
python compact_CREs.py                        # → compact_marginal_fp.txt
```

Override the output path with `$OUTPUT_TSV`. GPU strongly recommended; `CUDA_VISIBLE_DEVICES=0` is hardcoded.

## Output columns

`id` (`ref` or `del_<pos>`), `length`, `num_del`, `seq`, `motif_footprint`, `norm_footprint` (normalized to empty-motif control), `CRE`.

Deterministic given fixed inputs: subsample seed 42, strictly highest-`norm_footprint` deletion at each step (ties → earliest position by pandas stable sort), and deterministic chromBPNet inference.

## Requirements

Python ≥ 3.8 with `chrombpnet` installed (pulls `chrombpnet.training.utils.{losses,data_utils.get_seq,one_hot}`), TensorFlow/Keras, `pandas`, `numpy`, `biopython`, `pyfaidx`, `pyBigWig`, `deepdish`, `matplotlib`.
