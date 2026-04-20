#$ -S /bin/bash
#$ -cwd
#$ -N build_bt_idx
#$ -j yes
#$ -l mfree=6G,h_rt=96:00:00
#$ -pe serial 8
#$ -t 1-7

# SGE array job: build Bowtie2 indexes for all 7 reference FASTA files.
# Wraps 01_build_index.sh in array mode.

module load modules modules-init modules-gs
module load python/3.7.7
module load pysam
module load tbb/2020_U2
module load bowtie2/2.5.3
module load samtools/1.14

set -e

bash 01_build_index.sh config.sh $((SGE_TASK_ID - 1))

qstat -j $JOB_ID
