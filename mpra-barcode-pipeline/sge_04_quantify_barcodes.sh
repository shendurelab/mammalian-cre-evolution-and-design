#$ -S /bin/bash
#$ -cwd
#$ -N pileup_PYS2_oCREv2_MPRA_mBC
#$ -o log
#$ -e log
#$ -j yes
#$ -l mfree=108G,h_rt=24:00:00:00
#$ -t 1-7

# SGE array job: quantify mBC+UMI for 7 PYS2 oCREv2 MPRA libraries.
# Libraries: oCREv2_DNA_1-3, oCREv2_RNA_1-3, oCREv2_plasmid
# Wraps 04_quantify_barcodes.sh in array mode.

module load modules modules-init modules-gs
module load samtools/1.19
module load bedtools/2.31.1
module load gcc/8.1.0
module load pear/0.9.11
module load seqtk/1.4

set -e

mkdir -p log

bash 04_quantify_barcodes.sh config.sh $((SGE_TASK_ID - 1))

qstat -j $JOB_ID
