#!/usr/bin/env bash
#
# generate_lift_regions.sh
#
# Lift a mouse (Mus_musculus) BED of regions to one or more target species in
# the Zoonomia 241-mammalian Cactus alignment using halLiftover, then stitch
# the fragmented lifted intervals back together with stitchHalFrags_v2.R.
#
# This is a standalone bash version of the `generate_lift_regions` Snakemake
# rule from the original analysis pipeline, provided so that the per-species
# Cactus orthologs used in the paper can be regenerated independently of
# Snakemake.
#
# Requirements:
#   - halLiftover (from the HAL tools / Cactus conda env)
#   - Rscript with GenomicFeatures, GenomicRanges, rtracklayer, stringr, dplyr
#
# The 241-mammalian Cactus HAL file is NOT shipped in this repo (it is ~1TB);
# download it from:
#   https://cglgenomics.ucsc.edu/data/cactus/  (241-mammalian-2020v2.hal)
#
# This script ships with:
#   - stitchHalFrags_v2.R                       (post-processing helper)
#   - active_CRE_for_tiling.bed        (default source mm10 regions
#                                                used in the paper)
#
# Usage:
#   ./generate_lift_regions.sh \
#       -H /path/to/241-mammalian-2020v2.hal \
#       -s Homo_sapiens \
#       -o lifted_hits/
#
#   # Multiple species at once (space-separated list, or a file):
#   ./generate_lift_regions.sh -H ... -o ... -s "Homo_sapiens Rattus_norvegicus"
#   ./generate_lift_regions.sh -H ... -o ... -S species_list.txt
#
#   # To lift a different source BED, override with -b:
#   ./generate_lift_regions.sh -H ... -o ... -s Homo_sapiens -b my_regions.bed
#
# Options:
#   -H  Path to the Cactus .hal file                                  (required)
#   -o  Output directory for lifted BED files                        (required)
#   -s  Target species name, or quoted space-separated list of names
#   -S  File with one target species name per line
#   -b  Source BED file in the source genome
#         (default: active_CRE_for_tiling.bed next to this script)
#   -q  Source (query) species name in the HAL         (default: Mus_musculus)
#   -m  minMatch: min fraction of original interval covered    (default: 0.5)
#   -f  maxFrac:  max fraction of original interval covered    (default: 1.5)
#   -l  Path to halLiftover binary       (default: halLiftover on $PATH)
#   -r  Path to stitchHalFrags_v2.R      (default: next to this script)
#   -h  Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# defaults
SRC_SPECIES="Mus_musculus"
MIN_MATCH="0.5"
MAX_FRAC="1.5"
HAL_LIFTOVER_BIN="halLiftover"
STITCH_R="${SCRIPT_DIR}/stitchHalFrags_v2.R"
SRC_BED="${SCRIPT_DIR}/active_CRE_for_tiling.bed"
HAL_FILE=""
OUT_DIR=""
SPECIES_ARG=""
SPECIES_FILE=""

usage() { sed -n '2,55p' "$0"; exit 1; }

while getopts "H:b:o:s:S:q:m:f:l:r:h" opt; do
    case "${opt}" in
        H) HAL_FILE="${OPTARG}" ;;
        b) SRC_BED="${OPTARG}" ;;
        o) OUT_DIR="${OPTARG}" ;;
        s) SPECIES_ARG="${OPTARG}" ;;
        S) SPECIES_FILE="${OPTARG}" ;;
        q) SRC_SPECIES="${OPTARG}" ;;
        m) MIN_MATCH="${OPTARG}" ;;
        f) MAX_FRAC="${OPTARG}" ;;
        l) HAL_LIFTOVER_BIN="${OPTARG}" ;;
        r) STITCH_R="${OPTARG}" ;;
        h|*) usage ;;
    esac
done

# validate required args
if [[ -z "${HAL_FILE}" ]]; then
    echo "ERROR: -H (HAL file) is required" >&2
    usage
fi
if [[ -z "${OUT_DIR}" ]]; then
    echo "ERROR: -o (output directory) is required" >&2
    usage
fi
if [[ -z "${SPECIES_ARG}" && -z "${SPECIES_FILE}" ]]; then
    echo "ERROR: supply target species via -s or -S" >&2
    usage
fi

for f in "${HAL_FILE}" "${SRC_BED}" "${STITCH_R}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done

if ! command -v "${HAL_LIFTOVER_BIN}" >/dev/null 2>&1; then
    echo "ERROR: halLiftover not found: '${HAL_LIFTOVER_BIN}'" >&2
    echo "  Install HAL tools or activate your cactus conda env, or pass -l." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# build species list
SPECIES_LIST=()
if [[ -n "${SPECIES_ARG}" ]]; then
    # shellcheck disable=SC2206
    SPECIES_LIST+=( ${SPECIES_ARG} )
fi
if [[ -n "${SPECIES_FILE}" ]]; then
    if [[ ! -f "${SPECIES_FILE}" ]]; then
        echo "ERROR: species file not found: ${SPECIES_FILE}" >&2
        exit 1
    fi
    while IFS= read -r line; do
        line="${line// /}"
        [[ -n "${line}" ]] && SPECIES_LIST+=( "${line}" )
    done < "${SPECIES_FILE}"
fi

echo "[$(date)] Lifting ${#SPECIES_LIST[@]} target species from ${SRC_SPECIES}"
echo "  HAL         : ${HAL_FILE}"
echo "  source BED  : ${SRC_BED}"
echo "  output dir  : ${OUT_DIR}"
echo "  minMatch    : ${MIN_MATCH}"
echo "  maxFrac     : ${MAX_FRAC}"

for sp in "${SPECIES_LIST[@]}"; do
    if [[ "${sp}" == "${SRC_SPECIES}" ]]; then
        echo "[$(date)] Skipping source species ${sp}"
        continue
    fi
    out_bed="${OUT_DIR}/${sp}_active_CRE_lift.bed"
    echo "[$(date)] >>> ${sp}"

    "${HAL_LIFTOVER_BIN}" --noDupes --keepExtra \
        "${HAL_FILE}" "${SRC_SPECIES}" "${SRC_BED}" \
        "${sp}" "${out_bed}"

    Rscript --vanilla "${STITCH_R}" \
        "${SRC_BED}" "${out_bed}" "${out_bed}" \
        "${MIN_MATCH}" "${MAX_FRAC}"

    echo "[$(date)] <<< ${sp} -> ${out_bed}"
done

echo "[$(date)] Done."
