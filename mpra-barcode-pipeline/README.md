# MPRA Barcode Pipeline

End-to-end pipeline for MPRA (Massively Parallel Reporter Assay) barcode dictionary construction and quantification. Takes designed CRE sequences and raw sequencing FASTQs as input, and outputs per-CRE MPRA activity measurements.

Configured for the **PYS2 oCREv2 / CRE trajectory MPRA** experiment.

## Overview

```
Step 1:  Build Bowtie2 index from CRE FASTA(s)
Step 2:  Align subassembly reads -> build BC-CRE dictionary
Step 3:  Filter dictionary by read count
Step 4:  Quantify mBC + UMI from MPRA expression FASTQs
Step 5:  Merge dictionaries + join with quantification
Step 5b: Add duplicate CRE sequences to count table
Step 6:  QC, controls, replicate correlation, oCRE activity analysis
```

```
                  +----------------+
                  |  CRE FASTA(s)  |
                  +-------+--------+
                          |
               +----------v-----------+
               |  01_build_index.sh   |   (supports SGE array for multiple refs)
               +----------+-----------+
                           |
    Dictionary FASTQs      |
    (R1,R2=BC,R3)          |
         |    +------------v--------------+
         +--->|  02_build_dictionary.sh   |
              +------------+--------------+
                           |
              +------------v---------------+
              |  03_filter_dictionary.R    |
              +------------+---------------+
                           |
                      BC Dictionary
                           |
    MPRA FASTQs            |
    (R1=BC,R2=RC_BC,       |
     R3=UMI)               |
         |    +------------v------------------+
         +--->|  04_quantify_barcodes.sh      |   (SGE array: 7 libraries)
              +------------+------------------+
                           |
              +------------v---------------+
              |   05_merge_and_count.R     |   (merges 6 BC dictionaries)
              +------------+---------------+
                           |
              +------------v-------------------+
              |  05b_add_dups_in_oCREs.py/.ipynb|  (handles duplicate CREs)
              +------------+-------------------+
                           |
              +------------v---------------+
              |   06_analyze_mpra.R        |   (QC + oCRE activity + phylo)
              +----------------------------+
```

## PYS2 configuration

The pipeline is configured for the PYS2 oCREv2 MPRA experiment:

- **7 MPRA libraries**: oCREv2_DNA_1-3, oCREv2_RNA_1-3, oCREv2_plasmid
- **7 Bowtie2 indexes**: all_oCRE_cactus, compact, CRE_trajectory, evo_model, evo_phyloP, oCREs_Gata4_Human, pCREs_v2
- **6 BC dictionaries**: oCRE, CRE_traj, full_CRE, Engreitz_control, neg_control (minP), pos_control (EEF1A1p)
- **Sample naming**: `{lib}_{material}_{rep}_S{n}` (e.g., `oCREv2_DNA_1_S4`)

## Quick start

### 1. Configure

All paths are set in `config.sh`. Edit if data locations change.

### 2. Run

```bash
# Step 1: Build bowtie2 indexes (SGE array for all 7 references)
qsub sge_01_build_indexes.sh
# Or single index:
bash 01_build_index.sh

# Step 2: Build BC-CRE dictionary from subassembly sequencing
bash 02_build_dictionary.sh

# Step 2b (compact library only): length-stratified BC-CRE dictionary
# (Iterates lengths 40..300; parallelize via sge_02b_build_dictionary_compact.sh)
bash 02b_build_dictionary_compact.sh

# Step 3: Filter dictionary (min 10 reads)
Rscript 03_filter_dictionary.R \
    output/dict_pileups/sample1_full_pileup_k1_*.txt.gz \
    10 \
    output/PYS2_BC_dictionary.txt

# Step 4: Quantify barcodes (SGE array for 7 libraries)
qsub sge_04_quantify_barcodes.sh
# Or sequential:
bash 04_quantify_barcodes.sh

# Step 5: Merge 6 dictionaries and join with BC counts
Rscript 05_merge_and_count.R \
    output/bc_quant \
    output/PYS2_oCREv2_CRE_traj_MPRA

# Step 5b: Add duplicate CRE sequences
python 05b_add_dups_in_oCREs.py \
    output/PYS2_oCREv2_CRE_traj_MPRA_count_table_wo_zeroes.txt.gz \
    /net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/all_oCRE_cactusv2_reoriented_20240417_marked_dups.txt \
    output/PYS2_oCREv2
# Or interactively via the Jupyter notebook:
#   jupyter notebook 05b_add_dups_in_oCREs.ipynb

# Step 6: First-pass analysis (QC, controls, replicate correlation, oCRE activity)
Rscript 06_analyze_mpra.R \
    output/PYS2_oCREv2_count_table_wo_zeroes_w_dups.txt.gz \
    plots/
```

## Input file formats

### CRE FASTA (`CRE_FASTA`)
Standard FASTA with one entry per designed element:
```
>Bend5_chr4_8201__Mus_musculus_20240416.fa__chr4:8201-8500
ACGTACGT...
```

### Subassembly FASTQs (dictionary building)
Three-read sequencing of the CRE-BC library:
- **R1**: CRE read 1 (forward)
- **R2**: Barcode read (16bp default)
- **R3**: CRE read 2 (reverse)

Files named: `{sample}_R1_001.fastq.gz`, `{sample}_R2_001.fastq.gz`, `{sample}_R3_001.fastq.gz`

### MPRA expression FASTQs (quantification)
Three-read sequencing of MPRA RNA/DNA samples:
- **R1**: mBC forward (15bp)
- **R2**: mBC reverse (15bp)
- **R3**: UMI (10bp)

Files named: `{sample}_R1_001.fastq.gz`, `{sample}_R2_001.fastq.gz`, `{sample}_R3_001.fastq.gz`

PYS2 sample naming convention: `{lib}_{material}_{rep}_S{n}` (e.g., `oCREv2_DNA_1_S4`).
The pipeline parses `material` and `rep` from positions 2, 3 of the underscore-split name.

## Installation

### Option A: Conda (recommended)
```bash
conda env create -f environment.yml
conda activate mpra-pipeline
```

### Option B: Manual
```bash
# CLI tools
conda install -c bioconda bowtie2 samtools pear

# Python
pip install -r requirements.txt   # pysam, biopython

# R packages
Rscript install_r_packages.R
```

### Option C: Module load (SGE cluster)
```bash
module load samtools/1.19 bedtools/2.31.1 pear/0.9.11 bowtie2/2.5.3
```

## Dependencies

### CLI tools
| Tool | Version | Purpose |
|------|---------|---------|
| bowtie2 | >=2.5 | Align CRE reads to build BC dictionary |
| samtools | >=1.14 | BAM sorting, indexing, filtering |
| PEAR | >=0.9.11 | Error-correct overlapping BC reads |

### Python (>=3.6)
| Package | Purpose |
|---------|---------|
| pysam | Parse BAM alignments (dictionary step) |
| biopython | Parse FASTQ files (Bio.SeqIO) |
| pandas | Duplicate CRE handling (Step 5b) |

### R (>=4.3)
| Package | Used in | Purpose |
|---------|---------|---------|
| dplyr | Steps 3, 5, scripts | Data wrangling |
| tidyverse | Step 6, scripts | Data wrangling + plotting |
| ggplot2 | Steps 3, 6 | Plotting |
| cowplot | Step 6 | Plot themes |
| stringr | Step 5, scripts | String manipulation |
| reshape2 | Step 6 | Wide/long table reshaping |
| ggrastr | Step 6, mpra_tablemaker.R | Rasterized scatter plots |
| DescTools | mpra_tablemaker.R | Winsorization of MPRA activity |
| castor | Step 6 | Phylogenetic analysis |
| ape | Step 6 | Tree reading and distance calculation |
| phytools | Step 6 | Phylogenetic utilities |

## Outputs

| Step | Output | Description |
|------|--------|-------------|
| 1 | `bt_idx/*` | Bowtie2 index files (7 sets) |
| 2 | `output/dict_pileups/*_full_pileup_k1_*.txt.gz` | Raw BC-CRE read pileups |
| 3 | `output/*_BC_dictionary.txt` | Filtered BC-CRE mapping |
| 4 | `output/bc_quant/*_BC_quant_*.txt.gz` | Per-sample mBC quantification (UMI + reads) |
| 5 | `output/*_count_table.txt.gz` | BC x sample count matrix |
| 5 | `output/*_recovery_per_lib.txt.gz` | BC recovery rates per library |
| 5b | `output/*_count_table_wo_zeroes_w_dups.txt.gz` | Count table with duplicate CREs added |
| 5b | `output/*_duplicate_id_collapsed.txt` | CRE_id to quantification_id mapping |
| 6 | `plots/lib_saturation_calculations.txt` | Sequencing saturation per library |
| 6 | `plots/dna_barcode_recovery.pdf` | DNA BC recovery histogram |
| 6 | `plots/all_class_expression_hex_v2.pdf` | DNA vs RNA hex plots by class |
| 6 | `plots/control_expression.pdf` | Control QC scatter plots |
| 6 | `plots/cactus_biol_rep_corr_mpra_act.pdf` | Biological replicate correlation |
| 6 | `plots/oCRE_cre_hist.pdf` | oCRE activity by CRE |
| 6 | `plots/oCRE_cre_dist.pdf` | Activity vs phylogenetic distance |
| 6 | `plots/mouse_expression_bxplt.pdf` | Mouse element activity boxplot |
| 6 | `plots/PYS2_oCREv2_cactus_tiles_MPRA_mean_act.txt` | Mean oCRE activity table |

## SGE cluster usage

Pre-built SGE wrapper scripts for the PYS2 experiment:

```bash
# Build all 7 Bowtie2 indexes in parallel (array job)
qsub sge_01_build_indexes.sh

# Quantify all 7 MPRA libraries in parallel (array job, 108GB/task)
qsub sge_04_quantify_barcodes.sh
```

Steps 2 and 4 also support manual array jobs:

```bash
# In your SGE script with #$ -t 1-N:
bash 02_build_dictionary.sh config.sh $((SGE_TASK_ID - 1))
bash 04_quantify_barcodes.sh config.sh $((SGE_TASK_ID - 1))
```

## Directory structure

```
mpra-barcode-pipeline/
├── config.sh                        # Configuration (PYS2 paths & params)
├── README.md
├── environment.yml                  # Conda environment
├── requirements.txt                 # Python dependencies
├── install_r_packages.R             # R package installer
│
├── 01_build_index.sh                # Step 1: Build Bowtie2 index(es)
├── 02_build_dictionary.sh           # Step 2: BC-CRE dictionary (PE CREs)
├── 02b_build_dictionary_compact.sh  # Step 2b: BC-CRE dictionary (compact,
│                                    #          length-stratified, SE+PEAR)
├── 03_filter_dictionary.R           # Step 3: Filter dictionary
├── 04_quantify_barcodes.sh          # Step 4: MPRA BC quantification
├── 05_merge_and_count.R             # Step 5: Merge dictionaries + counts
├── 05b_add_dups_in_oCREs.py         # Step 5b: Duplicate CRE handling (script)
├── 05b_add_dups_in_oCREs.ipynb      # Step 5b: Duplicate CRE handling (notebook)
├── 06_analyze_mpra.R                # Step 6: QC + oCRE activity analysis
│
├── sge_01_build_indexes.sh          # SGE wrapper: array job for Step 1
├── sge_02b_build_dictionary_compact.sh  # SGE wrapper: array job for Step 2b
├── sge_04_quantify_barcodes.sh      # SGE wrapper: array job for Step 4
│
└── scripts/
    ├── subassembly_1BC_2readsCRE.py       # BC extraction from PE alignments
    ├── subassembly_1BC_SE_CRE.py          # BC extraction from SE alignments (compact)
    ├── pileup_BC_pair_CRE.R               # BC-CRE pileup (PE, keeps tagmentation)
    ├── pileup_BC_pair_CRE_SE.R            # BC-CRE pileup (SE, used by compact)
    ├── reformat_pear_outputs.py           # PEAR output reformatting
    ├── parse_fastqs_w_reformatted_pear.py # UMI linking
    ├── condense_MPRA_BC_file.R            # BC+UMI condensing
    ├── pileup_MPRA_mBC_noUMI_correction.R # mBC pileup
    └── mpra_tablemaker.R                  # MPRA activity helper functions
```

## Compact library (Step 2b)

The **compact CRE library** is produced by greedy 1 bp-deletion compaction
of mouse CRE tiles (see `scripts/oligo_design/compact_optimization/`). Each
compacted CRE ends up at a different final length, typically 40–300 bp.
Dictionary association for this library uses `02b_build_dictionary_compact.sh`
instead of `02_build_dictionary.sh` because:

- A single pooled Bowtie2 index would let a short compact CRE spuriously
  align to a prefix of a longer CRE. Solution: build one index **per exact
  CRE length** and align reads against only length-matched references.
- Upstream, CRE reads are **PEAR-assembled** (R1+R3 overlapping paired-end
  reads merged into a single assembled read) and adapter-trimmed, so
  Bowtie2 runs in single-end `--end-to-end` mode (not PE `--local`).
- Reads are subset to exactly one length at a time using `seqtk comp` so
  the PEAR-assembled CRE read FASTQ and the raw BC (R2) FASTQ stay aligned
  by read ID.

For each length *L* ∈ 40..300 the script:

1. `seqtk subseq` on `COMPACT_REF_FASTA` → CREs of length *L*; build
   `bowtie2-build` index on just those.
2. `seqtk subseq` on `COMPACT_PEAR_FASTQ` + `COMPACT_BC_FASTQ` → reads of
   length *L*.
3. `bowtie2 --end-to-end` → per-length BAM.
4. `scripts/subassembly_1BC_SE_CRE.py` (single-end variant, drops the
   `is_proper_pair` check that the PE version enforces) → extract 16 bp BC
   from R2 and join to the alignment row by read ID.
5. `scripts/pileup_BC_pair_CRE_SE.R` (configurable `min_mapq`, simpler
   groupby — no tagmentation position info because PEAR-assembled reads
   are full CRE length) → per-length `BC1 × CRE_id` count table.
6. Intermediate FASTQs, BAMs, and indexes are deleted; only the gzipped
   per-length pileup is kept under `output/compact/pileups/`.

The 261 per-length pileups are then concatenated and filtered in Step 3
exactly like the PE pipeline.

**Inputs (set in `config.sh`)**

| Variable | What |
|---|---|
| `COMPACT_REF_FASTA` | FASTA of all compact CRE sequences (varying lengths). |
| `COMPACT_PEAR_FASTQ` | PEAR-assembled + adapter-trimmed CRE read FASTQ. |
| `COMPACT_BC_FASTQ` | Raw R2 (barcode) FASTQ. Read IDs must match `COMPACT_PEAR_FASTQ`. |
| `COMPACT_MIN_MAPQ` | Min Bowtie2 MAPQ to keep in the pileup (default 2). |
| `COMPACT_LENGTHS` | Optional: override default 40..300 range. |

**Prerequisite: PEAR + trim_galore**

If you only have raw PE CRE reads, run PEAR to merge R1+R3 and trim_galore
to remove adapters before invoking this step:

```bash
pear -f compact_S3_R1_001.fastq.gz -r compact_S3_R3_001.fastq.gz \
     -o compact_S3_pear -j 8
trim_galore compact_S3_pear.assembled.fastq \
     -o trim_galore_fastqs/ --length 40
```

**Usage**

```bash
# Sequential (261 lengths, slow; use the SGE wrapper instead on a cluster)
bash 02b_build_dictionary_compact.sh

# Single length (debugging or SGE array)
bash 02b_build_dictionary_compact.sh config.sh 150

# SGE array — one task per length
qsub sge_02b_build_dictionary_compact.sh
```

**Output**

```
output/compact/pileups/compact_<L>bp_pileup_e2e_<date>.txt.gz
```

Each per-length pileup has columns: `BC1`, `CRE_id`, `n_count`, `n_prct`.
Concatenate and feed to `03_filter_dictionary.R`.

