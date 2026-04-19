# Oligo Design Scripts

Scripts that use a trained chromBPNet model to design synthetic CRE oligo libraries for MPRA testing. The trained endoderm model ships in the repo at `data/chrombpnet_models/endo/chrombpnet_nobias.h5` (see `scripts/chrombpnet_modeling/` for how it was trained).

## Dependencies

```
tensorflow chrombpnet pyfaidx pyBigWig deepdish biopython
pandas numpy matplotlib python-Levenshtein sequence_align
```

The mm10 FASTA is not shipped — set `$MM10_FASTA` before running any script.

`evolvability_optimization/evolvability.py` optionally imports `predict_chrombpnet` from the companion repo <https://github.com/zeitlingerlab/Brennan_Zelda_2023/tree/master/analysis/scripts/py>. Clone it and either set `$BRENNAN_ZELDA_2023_PY` or edit the fallback in `evolvability.py`. If unavailable, `predict_chrombpnet` falls back to `None` and callers that only need `load_model_wrapper` / `get_seq` still work.

## Scripts

### `evolvability_optimization/`

In silico evolutionary analysis using chromBPNet marginal footprinting.

| Script | What |
|---|---|
| `dms_marginal_footprinting.py` | Full DMS of mouse CRE tiles (all 3N single-base mutations) → log2 FC vs WT. |
| `pairwise_ism_snp_rescue_stepwise.py` | Greedy stepwise SNP rescue along the ancestral mouse lineage for 3 CREs — picks the highest-impact SNP at each step. |
| `pairwise_ism_snp_rescue_random_walk.py` | Random-walk null for the stepwise analysis — 50 random SNP orderings per CRE pair; 17 pairs across 5 CREs. |

Helpers: `evolvability.py` (model loading, ISM, generation selection), `iterative_evo_mutations.py` (mutation tracking across indels), `pairwise_extract_func.py` (NW alignment, mutation extraction, phyloP mapping).

### `compact_optimization/`

`compact_CREs.py` — iterative greedy 1 bp-deletion compaction. At each step evaluates every single-base deletion and keeps the one with the highest predicted footprint; repeats until the sequence is empty. Produces per-step activity trajectories used to pick minimal-length CRE designs. See `compact_optimization/README.md`.

### `derivatization_CRE/`

`shuffleCREFasta_chromTiles.py` — generate the derivatized CRE library testing TFBS sufficiency / arrangement / context. Per CRE tile, emits dinucleotide shuffles (whole and TFBS-only), CRE thripsis (5/10/20-break random re-assembly), TFBS reconstitutions into ~300 backgrounds at original spacing, random-position depositions, high-replicate depositions into 6 curated backgrounds, and ±1/±2 bp TFBS-flank variants. Output: `tiles_300bp_subtile_shuffled.tsv`. See `derivatization_CRE/README.md`.
