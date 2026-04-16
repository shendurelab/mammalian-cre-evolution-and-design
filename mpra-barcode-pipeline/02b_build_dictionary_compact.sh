#!/usr/bin/env bash
# Step 2b: Build BC -> compact-CRE dictionary (length-stratified)
#
# The compact CRE library is produced by greedy 1bp-deletion compaction of
# mouse CRE tiles (see scripts/oligo_design/compact_optimization/). Because
# different CREs are compacted to different final lengths (typically 40-300
# bp), a short compact CRE can spuriously map to a prefix of a longer one if
# everything is indexed together. To avoid this we do a length-stratified
# association:
#
#   for each length L in {40..300}:
#       - subset the compact CRE FASTA to entries of length exactly L
#       - build a bowtie2 index on just those
#       - subset the PEAR-assembled read FASTQ + the matching BC FASTQ to
#         reads of length exactly L
#       - bowtie2 --end-to-end alignment (reads are PEAR-assembled to full
#         CRE length, so end-to-end is the correct setting)
#       - extract barcodes from the BC read and join with the alignment
#         (subassembly_1BC_SE_CRE.py)
#       - pile up BC-CRE pairs (pileup_BC_pair_CRE_SE.R)
#
# Inputs must be prepared upstream:
#   - A PEAR-assembled + adapter-trimmed FASTQ of the CRE reads
#     (e.g. compact_S3_pear_*.assembled_trimmed.fq.gz from PEAR + trim_galore)
#   - The matching raw R2 (barcode) FASTQ
# The two FASTQs must have the same read IDs so that per-length subsets
# remain aligned between CRE and BC.
#
# Usage:
#   bash 02b_build_dictionary_compact.sh [config.sh] [length]
#     - config.sh: pipeline config (default: ./config.sh)
#     - length:    single length to process (for SGE array / debugging);
#                  if omitted, iterates over COMPACT_LENGTHS
#
# SGE array example:
#   #$ -t 40-300
#   bash 02b_build_dictionary_compact.sh config.sh "$SGE_TASK_ID"

set -euo pipefail
source "${1:-config.sh}"

# ---- Required config variables (define in config.sh) -----------------------
: "${COMPACT_REF_FASTA:?COMPACT_REF_FASTA must be set in config.sh}"
: "${COMPACT_PEAR_FASTQ:?COMPACT_PEAR_FASTQ must be set in config.sh}"
: "${COMPACT_BC_FASTQ:?COMPACT_BC_FASTQ must be set in config.sh}"
: "${BC_START:?}"
: "${BC_END:?}"
: "${OUT_DIR:?}"
: "${SCRIPT_DIR:?}"
: "${THREADS:=8}"

COMPACT_LENGTHS_DEFAULT=$(seq 40 300)
COMPACT_LENGTHS="${COMPACT_LENGTHS:-$COMPACT_LENGTHS_DEFAULT}"
COMPACT_MIN_MAPQ="${COMPACT_MIN_MAPQ:-2}"

COMPACT_DIR="${OUT_DIR}/compact"
IDX_DIR="${COMPACT_DIR}/bt_idx"
READS_DIR="${COMPACT_DIR}/subset_reads"
ALIGN_DIR="${COMPACT_DIR}/alignments"
MERGE_DIR="${COMPACT_DIR}/merged"
PILEUP_DIR="${COMPACT_DIR}/pileups"
mkdir -p "$IDX_DIR" "$READS_DIR" "$ALIGN_DIR" "$MERGE_DIR" "$PILEUP_DIR"

date_str=$(date "+%Y%m%d")

process_length() {
    local L="$1"
    echo "=== [compact] length = ${L} bp ==="

    local name_list="${READS_DIR}/${L}_name.list"
    local subset_ref="${IDX_DIR}/${L}/ref_${L}.fa"
    local subset_fq="${READS_DIR}/${L}_reads.fq.gz"
    local subset_bc_fq="${READS_DIR}/${L}_bc_reads.fq.gz"
    local idx_prefix="${IDX_DIR}/${L}/${L}"
    local bam="${ALIGN_DIR}/compact_${L}bp.bam"
    local merged="${MERGE_DIR}/compact_${L}bp_merged_e2e_${date_str}.txt.gz"
    local pileup="${PILEUP_DIR}/compact_${L}bp_pileup_e2e_${date_str}.txt"

    mkdir -p "${IDX_DIR}/${L}"

    # 1) Subset reference FASTA to entries of length L, build per-length index
    echo "  Building length-${L} bowtie2 index..."
    seqtk comp "$COMPACT_REF_FASTA" \
        | awk -v L="$L" '$2 == L {print $1}' \
        > "$name_list"

    if [[ ! -s "$name_list" ]]; then
        echo "  No CREs of length ${L} in reference; skipping."
        rm -f "$name_list"
        return 0
    fi

    seqtk subseq "$COMPACT_REF_FASTA" "$name_list" > "$subset_ref"
    bowtie2-build --threads "$THREADS" --seed 42 \
        "$subset_ref" "$idx_prefix" > /dev/null

    # 2) Subset reads (CRE + matching BC) to those of length L
    echo "  Subsetting reads of length ${L}..."
    seqtk comp "$COMPACT_PEAR_FASTQ" \
        | awk -v L="$L" '$2 == L {print $1}' \
        > "$name_list"

    if [[ ! -s "$name_list" ]]; then
        echo "  No reads of length ${L}; skipping."
        rm -rf "${IDX_DIR}/${L}" "$name_list"
        return 0
    fi

    seqtk subseq "$COMPACT_PEAR_FASTQ" "$name_list" | gzip > "$subset_fq"
    seqtk subseq "$COMPACT_BC_FASTQ"   "$name_list" | gzip > "$subset_bc_fq"

    # 3) End-to-end bowtie2 alignment
    echo "  Aligning (end-to-end)..."
    bowtie2 --threads "$THREADS" --end-to-end \
        -x "$idx_prefix" \
        -U "$subset_fq" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"

    # 4) Extract BCs and join with alignments
    echo "  Extracting barcodes..."
    python "${SCRIPT_DIR}/subassembly_1BC_SE_CRE.py" \
        --input_bam "$bam" \
        --input_BC1_fastq "$subset_bc_fq" \
        --output_file "$merged" \
        --barcode1_start "$BC_START" \
        --barcode1_end "$BC_END" \
        --search_seq1 no_seq

    # 5) Pileup BC-CRE pairs
    echo "  Piling up BC-CRE pairs (min_mapq=${COMPACT_MIN_MAPQ})..."
    Rscript --vanilla "${SCRIPT_DIR}/pileup_BC_pair_CRE_SE.R" \
        "$merged" "$pileup" "$COMPACT_MIN_MAPQ"
    gzip -f "$pileup"

    # 6) Cleanup intermediates (keep pileup)
    rm -rf "${IDX_DIR}/${L}"
    rm -f "$name_list" "$subset_fq" "$subset_bc_fq" "$bam" "${bam}.bai"

    echo "  Done: ${pileup}.gz"
}

if [[ -n "${2:-}" ]]; then
    process_length "$2"
else
    for L in $COMPACT_LENGTHS; do
        process_length "$L"
    done
fi
