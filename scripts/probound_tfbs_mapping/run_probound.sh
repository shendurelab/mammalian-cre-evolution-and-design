#!/usr/bin/env bash
# Run ProBound TFBS mapping on a CRE FASTA.
#
# Inputs (positional):
#   $1: input FASTA of CRE sequences
#   $2: output TSV.gz path
#   $3: (optional) number of threads (default 8)
#
# Required env:
#   PROBOUND_JAR  path to ProBound-jar-with-dependencies.jar
#                 (build from https://github.com/RubeLab/ProBound)
#
# Uses data/ParEndo_markerTF_cuttoffs.txt as the TF model table (17 TFs).
# Downstream scripts filter to the 7 TFs actually plotted.
#
# Example:
#   PROBOUND_JAR=/path/to/ProBound-jar-with-dependencies.jar \
#     bash scripts/probound_tfbs_mapping/run_probound.sh \
#          my_seqs.fa my_TFBS_predictions.txt.gz 8

set -euo pipefail

FASTA="${1:?missing FASTA argument}"
OUT="${2:?missing output path argument}"
THREADS="${3:-8}"

: "${PROBOUND_JAR:?set PROBOUND_JAR to the ProBound jar path}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TF_MODEL="${REPO_ROOT}/data/ParEndo_markerTF_cuttoffs.txt"
SCRIPT="${REPO_ROOT}/scripts/probound_tfbs_mapping/ProBound_TFBS_mapping.py"

python "$SCRIPT" -f "$FASTA" -m "$TF_MODEL" -o "$OUT" -t "$THREADS"
