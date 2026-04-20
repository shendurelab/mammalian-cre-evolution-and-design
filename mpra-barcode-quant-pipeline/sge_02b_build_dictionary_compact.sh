#!/bin/bash
#$ -N compact_dict
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -o logs/
#$ -l mfree=16G,h_rt=8:00:00
#$ -pe serial 8
#$ -t 40-300
# SGE array job: one task per CRE length (40..300 bp).
# Each task builds the per-length bowtie2 index, subsets reads, aligns,
# extracts BCs, and piles up BC-CRE pairs.

mkdir -p logs

module load seqtk/1.4
module load bowtie2/2.5.3
module load samtools/1.14

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PIPELINE_DIR"

bash 02b_build_dictionary_compact.sh config.sh "$SGE_TASK_ID"
