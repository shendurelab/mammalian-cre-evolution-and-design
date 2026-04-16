#!/usr/bin/env bash
# Step 2: Build BC-CRE dictionary from subassembly sequencing
#
# For each sample: align paired CRE reads -> extract barcodes -> pileup BC-CRE pairs
#
# Usage: bash 02_build_dictionary.sh [config.sh] [sample_index]
#   sample_index: 0-based index into DICT_SAMPLES array (for SGE array jobs)
#                 If omitted, processes all samples sequentially.
set -euo pipefail
source "${1:-config.sh}"

ALIGN_DIR="${OUT_DIR}/dict_alignments"
MERGE_DIR="${OUT_DIR}/dict_merged"
PILEUP_DIR="${OUT_DIR}/dict_pileups"
mkdir -p "$ALIGN_DIR" "$MERGE_DIR" "$PILEUP_DIR"

date_str=$(date "+%Y%m%d")

process_sample() {
    local sample="$1"
    echo "=== Processing dictionary sample: $sample ==="

    local r1="${DICT_FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    local r2_bc="${DICT_FASTQ_DIR}/${sample}_R2_001.fastq.gz"
    local r3="${DICT_FASTQ_DIR}/${sample}_R3_001.fastq.gz"
    local bam="${ALIGN_DIR}/${sample}_paired_k1_${date_str}.bam"
    local unaligned="${ALIGN_DIR}/${sample}_unaligned_k1_${date_str}.fastq.gz"

    # 1) Align paired CRE reads (R1 + R3) to CRE index
    echo "  Aligning..."
    bowtie2 --threads "$THREADS" -k 1 \
        -x "$BT_INDEX" --local \
        -1 "$r1" -2 "$r3" \
        --un-gz "$unaligned" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"

    # 2) Extract barcodes from R2 and link to aligned CRE
    local merged="${MERGE_DIR}/${sample}_merged_k1_${date_str}.txt.gz"
    echo "  Extracting barcodes..."
    python "${SCRIPT_DIR}/subassembly_1BC_2readsCRE.py" \
        --input_bam "$bam" \
        --input_BC1_fastq "$r2_bc" \
        --output_file "$merged" \
        --barcode1_start "$BC_START" \
        --barcode1_end "$BC_END" \
        --search_seq1 no_seq

    # 3) Pileup BC-CRE associations
    local pileup="${PILEUP_DIR}/${sample}_full_pileup_k1_${date_str}.txt"
    echo "  Piling up BC-CRE pairs..."
    Rscript --vanilla "${SCRIPT_DIR}/pileup_BC_pair_CRE.R" \
        "$merged" "$pileup" 2
    gzip "$pileup"

    echo "  Done: ${pileup}.gz"
}

if [[ -n "${2:-}" ]]; then
    process_sample "${DICT_SAMPLES[$2]}"
else
    for sample in "${DICT_SAMPLES[@]}"; do
        process_sample "$sample"
    done
fi
