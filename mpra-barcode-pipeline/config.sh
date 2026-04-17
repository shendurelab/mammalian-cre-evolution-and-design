#!/usr/bin/env bash
# ============================================================================
# MPRA Barcode Pipeline — Configuration (PYS2 oCREv2 MPRA)
# ============================================================================
# Customized for the PYS2 oCREv2 / CRE trajectory MPRA experiment.
# All pipeline scripts source this file.
#
# Raw subassembly FASTQs (Step 2 input) are available from GEO accession
# GSE328309: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328309
# Download them into the DICT_FASTQ_DIR location below.
# MPRA expression FASTQs (Step 4) are not part of this GEO submission.
# ============================================================================

# --- Genome / references ---
GENOME_FA="/net/gs/vol1/home/tli824/bin/references/human/UCSC/hg38/hg38.fa"

# --- CRE library FASTAs (one entry per designed element) ---
# Multiple reference sets used for index building (Step 1).
# The primary index for oCRE dictionary alignment:
CRE_FASTA="/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/all_oCRE_cactusv2_reoriented_20240417_dedupped.fa"

# All reference FASTAs and their index names (used by 01_build_index.sh array mode)
CRE_FASTA_ARRAY=(
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/all_oCRE_cactusv2_reoriented_20240417_dedupped.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/compact.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/CRE_activity_trajectory.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/evo_model.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/evo_phyloP.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/evo_random.fa"
    "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/pCREs_v2.fa"
)
BT_INDEX_NAMES=(
    "all_oCRE_cactus"
    "compact"
    "CRE_trajectory"
    "evo_model"
    "evo_phyloP"
    "evo_random"
    "pCREs_v2"
)

# --- Bowtie2 index prefix (built in Step 1) ---
BT_INDEX="bt_idx/all_oCRE_cactus"

# --- Barcode dictionary FASTQs (subassembly sequencing) ---
# Paired-end: R1 = CRE read 1, R2 = barcode read, R3 = CRE read 2
DICT_FASTQ_DIR="/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615"
# Space-separated list of sample prefixes (files expected: {prefix}_R1_001.fastq.gz, etc.)
DICT_SAMPLES=("sample1")

# Barcode extraction parameters
BC_START=0
BC_END=15       # 0-indexed; extracts 16bp barcode (positions 0..15)

# Minimum read count to keep a BC-CRE pair in the dictionary
DICT_MIN_READS=10

# --- MPRA quantification FASTQs ---
# R1 = mBC forward (15bp), R2 = mBC reverse (15bp), R3 = UMI (10bp)
QUANT_FASTQ_DIR="/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/fastqs"
# PYS2 sample names: {lib}_{material}_{rep}_S{n}
QUANT_SAMPLES=("oCREv2_DNA_1_S4" "oCREv2_DNA_2_S5" "oCREv2_DNA_3_S6" \
               "oCREv2_RNA_1_S1" "oCREv2_RNA_2_S2" "oCREv2_RNA_3_S3" \
               "oCREv2_plasmid_S7")

# PEAR barcode length (used for -v, -m, -n, -t flags)
BC_LEN=15

# --- BC dictionary files (for Step 5: merge) ---
# Pre-built BC dictionaries with columns: BC (or BC1), CRE_id[, class]
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DICT_FILES=(
    "${PIPELINE_DIR}/data/dictionaries/CRE_traj_final_BC_list.txt.gz"
    "${PIPELINE_DIR}/data/dictionaries/oCRE_v2_final_BC_list.txt.gz"
    "${PIPELINE_DIR}/data/dictionaries/Engreitz_controls_final_BC_list.txt"
    "${PIPELINE_DIR}/data/dictionaries/minP_final_BC_list.txt.gz"
    "${PIPELINE_DIR}/data/dictionaries/EEF1aP_final_BC_list.txt.gz"
    "${PIPELINE_DIR}/data/dictionaries/fullCRE_BC_asso_10percent_complexity_lib_20231126.txt.gz"
)
# Class labels corresponding to each dictionary file above
DICT_CLASSES=("CRE_traj" "oCRE" "Engreitz_control" "neg_control" "pos_control" "full_CRE")

# --- Duplicate CRE reference (for Step 5b: add_dups_in_oCREs) ---
DUP_REF="${PIPELINE_DIR}/data/dictionaries/all_oCRE_cactusv2_reoriented_20240417_marked_dups.txt"

# --- Control dictionaries (optional; leave empty arrays if none) ---
CONTROL_DICT_FILES=()

# --- Compact library (02b_build_dictionary_compact.sh) -------------------
# The compact CRE library has length-variable entries (40-300 bp) from the
# greedy 1bp-deletion compaction pipeline. Dictionary building is
# length-stratified: one bowtie2 index + alignment per exact CRE length.
#
# Reference FASTA of all compact CRE sequences (varying lengths):
COMPACT_REF_FASTA="/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_fastqs_20240615/ref_fasta/compact.fa"
# PEAR-assembled + adapter-trimmed CRE read FASTQ (upstream PEAR + trim_galore):
COMPACT_PEAR_FASTQ="/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_libs_20240708/data/trim_galore_fastqs/compact_S3_pear_20240708.assembled_trimmed.fq.gz"
# Matching raw R2 (barcode) FASTQ (read IDs must match COMPACT_PEAR_FASTQ):
COMPACT_BC_FASTQ="/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/data/twist_v3_libs_20240708/data/fastqs/compact_S3_R2_001.fastq.gz"
# Lengths to process (space-separated). Defaults to 40..300 inside the script.
# COMPACT_LENGTHS="40 41 42 ... 300"
# Minimum bowtie2 MAPQ to keep in BC-CRE pileup:
COMPACT_MIN_MAPQ=2

# --- Output directories (auto-created) ---
OUT_DIR="output"

# --- Threads ---
THREADS=8

# --- Path to scripts/ directory (auto-detected) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
