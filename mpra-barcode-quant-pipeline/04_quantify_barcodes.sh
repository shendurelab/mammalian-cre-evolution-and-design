#!/usr/bin/env bash
# Step 4: Quantify mBC+UMI from MPRA expression FASTQs
#
# For each sample: PEAR error-correct BC -> parse + link UMI -> condense -> pileup per mBC
#
# Usage: bash 04_quantify_barcodes.sh [config.sh] [sample_index]
#   sample_index: 0-based index into QUANT_SAMPLES (for SGE array jobs)
#                 If omitted, processes all samples sequentially.
#
# For SGE array usage (PYS2: 7 libraries):
#   #$ -t 1-7
#   #$ -l mfree=108G,h_rt=24:00:00:00
#   bash 04_quantify_barcodes.sh config.sh $((SGE_TASK_ID - 1))
set -euo pipefail
source "${1:-config.sh}"

PEAR_DIR="${OUT_DIR}/pear_outs"
MERGE_DIR="${OUT_DIR}/quant_merged"
QUANT_DIR="${OUT_DIR}/bc_quant"
mkdir -p "$PEAR_DIR" "$MERGE_DIR" "$QUANT_DIR"

date_str=$(date "+%Y%m%d")

process_sample() {
    local sample="$1"
    echo "=== Quantifying sample: $sample ==="

    local r1="${QUANT_FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    local r2="${QUANT_FASTQ_DIR}/${sample}_R2_001.fastq.gz"
    local r3="${QUANT_FASTQ_DIR}/${sample}_R3_001.fastq.gz"

    # 1) PEAR: error-correct overlapping BC reads (R1 fwd + R2 rev)
    local pear_out="${PEAR_DIR}/${sample}_pear_BC_${date_str}"
    echo "  PEAR merging..."
    pear -j "$THREADS" \
         -v "$BC_LEN" -m "$BC_LEN" -n "$BC_LEN" -t "$BC_LEN" \
         -f "$r1" -r "$r2" \
         -o "$pear_out"
    gzip "${pear_out}.assembled.fastq"

    # 2) Reformat PEAR output
    local pear_seqs="${PEAR_DIR}/${sample}_pear_BC_seqs_${date_str}.txt.gz"
    echo "  Reformatting PEAR output..."
    python "${SCRIPT_DIR}/reformat_pear_outputs.py" \
        -i "${pear_out}.assembled.fastq.gz" \
        -o "$pear_seqs"

    # 3) Parse FASTQs: link error-corrected BC with UMI
    local valid_out="${MERGE_DIR}/${sample}_BC_valid_pear_simple_${date_str}.txt.gz"
    local nopear_out="${MERGE_DIR}/${sample}_BC_no_pear_simple_${date_str}.txt.gz"
    echo "  Parsing and linking UMI..."
    python "${SCRIPT_DIR}/parse_fastqs_w_reformatted_pear.py" \
        --pear_file "$pear_seqs" \
        --out_valid "$valid_out" \
        --out_no_pear_BC "$nopear_out" \
        --in_R1 "$r1" --in_R2 "$r2" --in_UMI "$r3" \
        --R1_name BC --R2_name RC_BC --UMI_name UMI \
        --simple_out_bool 1

    # 4) Condense: count reads per BC+UMI combination
    local condensed="${QUANT_DIR}/${sample}_BC_pear_UMI_condensed_${date_str}.txt"
    echo "  Condensing BC+UMI..."
    Rscript --vanilla "${SCRIPT_DIR}/condense_MPRA_BC_file.R" \
        "$valid_out" "$condensed" "BC_pear" "UMI"
    gzip "$condensed"

    # 5) Pileup: UMI + read counts per mBC
    local bc_quant="${QUANT_DIR}/${sample}_BC_quant_${date_str}.txt"
    echo "  Final pileup per mBC..."
    Rscript --vanilla "${SCRIPT_DIR}/pileup_MPRA_mBC_noUMI_correction.R" \
        "${condensed}.gz" "$bc_quant"
    gzip "$bc_quant"

    echo "  Done: ${bc_quant}.gz"
}

if [[ -n "${2:-}" ]]; then
    process_sample "${QUANT_SAMPLES[$2]}"
else
    for sample in "${QUANT_SAMPLES[@]}"; do
        process_sample "$sample"
    done
fi
