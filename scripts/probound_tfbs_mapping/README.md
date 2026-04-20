# ProBound TFBS mapping

Regenerates the full (unslimmed) ProBound TFBS prediction tables used in the paper. The slim versions shipped in `data/*_slim.txt.gz` are what the figure scripts actually read; this pipeline is for reproducing the complete, per-k-mer ProBound affinity table from scratch.

## Files

| File | Purpose |
|---|---|
| `ProBound_TFBS_mapping.py` | Per-k-mer TFBS scoring. For each FASTA record × each TF model, enumerates all k-mers at the TF's model size, scores them (forward + reverse) via ProBound, normalizes to each TF's `max_aff`, picks the stronger orientation, and writes a gzipped TSV. |
| `run_probound.sh` | Thin driver; passes the shipped TF model table to the Python script. |

## Requirements

- **Java 8+** and a built ProBound jar. Clone [RubeLab/ProBound](https://github.com/RubeLab/ProBound) and build the `ProBound-jar-with-dependencies.jar` (Maven). Expose its path as `$PROBOUND_JAR`.
- Python 3 with `biopython`, `pandas`, `numpy`, `seaborn`.
- ParEndo marker TF models — reuses the repo-local `data/ParEndo_markerTF_cuttoffs.txt` (17 TFs; superset of the 12 in the published runs).

## Usage

```bash
export PROBOUND_JAR=/path/to/ProBound-jar-with-dependencies.jar

# oCRE v2 sequences -> oCRE TFBS predictions (used by Fig 1/2 loaders)
bash scripts/probound_tfbs_mapping/run_probound.sh \
     /path/to/all_oCRE_cactusv2_reoriented.fa \
     data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz 8

# CRE-optimization trajectory sequences -> Fig 5 loader
bash scripts/probound_tfbs_mapping/run_probound.sh \
     /path/to/twist_in_silico_evolution_model_driven_order.fa \
     data/CRE_optimization_ProBound_ParEndo_TFBS_predictions.txt.gz 8
```

After regenerating either full table, rerun the corresponding slim preprocessor if you want to refresh the committed slim file:

```bash
Rscript scripts/utils/preprocess_oCRE_TFBS.R              # oCRE
Rscript scripts/utils/preprocess_CRE_optimization_TFBS.R  # CRE optimization
```

## Output schema

Gzipped TSV with columns: `start_pos, end_pos, kmer_seq, affinity_FOR, affinity_REV, TF_name, substring_l, CRE, norm_affinity_FOR, norm_affinity_REV, TFBS_orientation, norm_affinity`.

`norm_affinity` is the max of forward/reverse normalized affinities and is the single value used by the figure code.

## Input FASTAs

These are not shipped in this directory (several MB each; use the original lab paths or regenerate from the cactus alignment + oCRE coordinates):

- `all_oCRE_cactusv2_reoriented.fa` — reoriented oCRE sequences across 241 mammalian species
- `twist_CRE_activity_trajectory_order.fa`, `twist_in_silico_evolution_{model_driven,phyloP_walk,random_walk}_order.fa` — Twist oligo order sequences for the CRE-optimization trajectories

The sequence headers used become the `CRE` column in the output.
