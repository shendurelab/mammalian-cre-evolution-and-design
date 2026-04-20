#!/usr/bin/env bash
# Step 1: Build Bowtie2 index from CRE FASTA
#
# Supports two modes:
#   Single index:  bash 01_build_index.sh [config.sh]
#   Array (SGE):   bash 01_build_index.sh [config.sh] <task_index>
#     where task_index is 0-based into CRE_FASTA_ARRAY / BT_INDEX_NAMES
#
# For SGE array usage, submit with e.g.:
#   #$ -t 1-7
#   bash 01_build_index.sh config.sh $((SGE_TASK_ID - 1))
set -euo pipefail
source "${1:-config.sh}"

if [[ -n "${2:-}" ]]; then
    # Array mode: build a specific index
    idx="$2"
    fasta="${CRE_FASTA_ARRAY[$idx]}"
    index_name="bt_idx/${BT_INDEX_NAMES[$idx]}"
    mkdir -p "$(dirname "$index_name")"
    echo "Building Bowtie2 index [$idx]: $fasta -> $index_name"
    bowtie2-build --threads "$THREADS" --seed 42 "$fasta" "$index_name"
    echo "Done. Index files at: ${index_name}.*"
else
    # Single mode: build from CRE_FASTA
    mkdir -p "$(dirname "$BT_INDEX")"
    echo "Building Bowtie2 index from: $CRE_FASTA"
    echo "Index prefix: $BT_INDEX"
    bowtie2-build --threads "$THREADS" --seed 42 "$CRE_FASTA" "$BT_INDEX"
    echo "Done. Index files at: ${BT_INDEX}.*"
fi
