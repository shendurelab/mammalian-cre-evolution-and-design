# Oligo Design Scripts

Scripts for designing synthetic CRE (cis-regulatory element) oligo libraries using chromBPNet model predictions. These scripts use a trained chromBPNet model to predict chromatin accessibility (marginal footprint) and iteratively design or perturb CRE sequences for MPRA testing.

All scripts require a trained chromBPNet model (see `scripts/chrombpnet_modeling/`).

## Dependencies

### Python packages

```
tensorflow
chrombpnet
pyfaidx
pyBigWig
deepdish
biopython
pandas
numpy
matplotlib
python-Levenshtein
sequence_align
```

### External helper module

`evolvability_optimization/evolvability.py` optionally imports `predict_chrombpnet` from the companion Brennan_Zelda_2023 repo:
<https://github.com/zeitlingerlab/Brennan_Zelda_2023/tree/master/analysis/scripts/py>

Clone that repo and either set

```bash
export BRENNAN_ZELDA_2023_PY=/path/to/Brennan_Zelda_2023/analysis/scripts/py
```

or edit the fallback path inside `evolvability.py`. If the module cannot be imported, `predict_chrombpnet` falls back to `None` — callers that only need `load_model_wrapper` / `get_seq` still work.

### Data files

All scripts depend on:
- **chromBPNet model**: `data/chrombpnet_models/endo/chrombpnet_nobias.h5` (included in this repo)
- **mm10 genome**: `/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa`
- **Fold splits**: `EB_chrombpnet/splits/fold_0.json`
- **Background regions**: `cactus_CREs/bg_beds/endo_negatives_10ksubset.bed` (negative genomic regions for normalizing marginal footprints)

## Scripts

### evolvability_optimization/

Evaluates the effect of evolutionary mutations on CRE activity using chromBPNet marginal footprinting.

#### `dms_marginal_footprinting.py`

**Purpose**: In silico deep mutational scanning (DMS) of mouse CRE tiles. For each 300bp max-activity CRE tile, generates all possible single-nucleotide mutations and predicts their effect on chromatin accessibility using chromBPNet marginal footprinting.

**How it works**:
1. Loads the chromBPNet endoderm model
2. Computes a control (empty motif) marginal footprint on 500 background genomic sequences
3. For each CRE tile from the input FASTA:
   - Predicts wild-type marginal footprint (normalized to control)
   - Generates all 3N single-base mutations (N = tile length)
   - Predicts marginal footprint for each mutant
   - Computes log2 fold-change vs wild-type
4. Outputs a table of all mutations with predicted effects

**Input**: `DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa`

**Output**: `dms_max_300bptile_selected_CREs_marginal_footprint_v3.txt`

**Local dependencies**: `evolvability.py`

---

#### `pairwise_ism_snp_rescue_stepwise.py`

**Purpose**: Stepwise SNP rescue analysis along the ancestral lineage. For each pair of ancestral CRE sequences along the mouse lineage, identifies all pairwise SNPs/indels and applies them one at a time in order of predicted impact (greedy stepwise approach, highest-predicted-footprint-first).

**How it works**:
1. For each consecutive ancestral pair along the mouse lineage (e.g., fullTreeAnc239 -> fullTreeAnc55 -> ... -> Mus_musculus):
   - Aligns the reference and target sequences
   - Extracts all SNP/indel events
   - Predicts marginal footprint for each individual SNP rescue
   - Applies the highest-impact SNP first, re-aligns, and repeats
   - Tracks cumulative edit distance and predicted activity at each step
2. This reveals the order of mutations that most efficiently restore CRE activity

**CREs analyzed**: Gata4_chr14_5729, Epas1_chr17_10063, Sparc_chr11_7211

**Output**: `ancestral_lineage_mouse_snps_marginal_fp_stepwise.txt`

**Local dependencies**: `evolvability.py`, `iterative_evo_mutations.py`, `pairwise_extract_func.py`

---

#### `pairwise_ism_snp_rescue_random_walk.py`

**Purpose**: Random walk SNP rescue analysis between pairs of CRE orthologs. Same as stepwise, but instead of picking the best SNP at each step, randomly selects one mutation per step (repeated 50 times per pair). This provides a null distribution for comparison with the stepwise approach.

**How it works**:
1. For each CRE pair (e.g., Mus_musculus vs Homo_sapiens, or Mus_musculus vs fullTreeAnc239):
   - Aligns reference and target sequences
   - For 50 random iterations:
     - Randomly selects one SNP/indel, applies it
     - Re-aligns, randomly selects next, repeats until sequences converge
   - Adds phyloP conservation scores when reference is Mus_musculus
2. Predicts marginal footprint at each step

**CRE pairs**: 17 pairs across 5 CREs (Gata4, Epas1, Sparc, Bend5, Lama1) between mouse, rat, human, and ancestral sequences

**Output**: `ancestral_lineage_random_walk_snps_marginal_fp_v2.txt`

**Local dependencies**: `evolvability.py`, `iterative_evo_mutations.py`, `pairwise_extract_func.py`

---

#### Helper modules

- **`evolvability.py`** - Core module for chromBPNet model loading, in silico mutagenesis, and evolutionary optimization functions (maximize/minimize/neutral next generation selection)
- **`iterative_evo_mutations.py`** - Handles iterative mutation tracking and coordinate updating after insertions/deletions
- **`pairwise_extract_func.py`** - Sequence alignment (Needleman-Wunsch), mutation extraction from alignments, phyloP score mapping, and DMS score mapping

---

### compact_optimization/

#### `compact_CREs.py`

**Purpose**: Iterative sequence compaction by greedy 1bp deletion. For each CRE tile, systematically evaluates all possible single-base deletions and removes the base whose deletion maximizes predicted activity. Repeats until the sequence is fully consumed.

**How it works**:
1. Loads the chromBPNet endoderm model
2. For each CRE tile:
   - Predicts reference (full-length) marginal footprint
   - Generates all possible 1bp deletions of the current sequence
   - Predicts marginal footprint for each deletion variant
   - Selects the deletion that maximizes predicted activity
   - Uses the shortened sequence as input for the next round
   - Repeats for all positions (N rounds for N-bp sequence)
3. Tracks the predicted activity trajectory as the sequence shrinks

**Use case**: Identifies which bases can be removed with minimal activity loss, enabling design of compact synthetic CREs.

**Output**: `compact_marginal_fp.txt`

**Local dependencies**: `evolvability.py` (via import of shared functions in the script)

---

### pCREs_v2/

#### `shuffleCREFasta_chromTiles.py`

**Purpose**: Generate permuted CRE (pCRE) oligo libraries by systematically rearranging CRE sequence features. Creates a large combinatorial library of synthetic sequences that test the sufficiency and arrangement of TFBS motifs.

**How it works**:
For each CRE tile, generates multiple classes of permuted sequences:
1. **Dinucleotide shuffles** (200x): Full-sequence shuffle preserving dinucleotide frequency
2. **TFBS shuffles** (200x): Shuffle only within TFBS regions
3. **CRE thripsis** (200x each at 5, 10, 20 break points): Cut sequence at random points, shuffle and reverse-complement segments while preserving TFBS
4. **Reconstitutions** (per background, 300+ backgrounds): Implant TFBS motifs at their original spacing into diverse background sequences
5. **Random depositions** (per background, 2 seeds): Implant TFBS at random non-overlapping positions in backgrounds
6. **Fixed-background depositions** (6 selected backgrounds x 100 seeds): High-replicate random depositions in permissive/intermediate/recalcitrant backgrounds
7. **Flanking reconstitutions** (+1bp and +2bp flanking around TFBS): Test the role of immediate TFBS context

**Input data**:
- CRE tile FASTA: `DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa`
- TFBS predictions from ProBound (Gata6, Sox17, Foxa2, Klf4 in ParEndo): `TFBS_ProBound_*_ParEndo_maxTileDMS_20240409.txt`
- Background sequences: `endo_bg_2k.fa` and previous pCRE v1 backgrounds

**Output**: `JB_tiles_300bp_subtile_shuffled_v4.tsv` (tab-delimited: header, sequence, tile_id)

**Local dependencies**: `SequencePremutationTools.py`

**`SequencePremutationTools.py`** - Class implementing dinucleotide-preserving shuffling (Altschul & Erickson algorithm via Eulerian path), interval-aware TFBS extraction, TFBS implantation into backgrounds, and random segment shuffling.
