#$ -S /bin/bash
#$ -cwd    ## use current working directory
#$ -j yes  ## merge stdout and stderr
#$ -l mfree=100G,h_rt=168:00:00
#$ -l gpgpu=1,cuda=1

conda activate chrombpnet

genome='/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa'
chromSizes='/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.chrom.sizes'


mkdir -p bias_models
mkdir -p chrombpnet_models

cell_array=( meso endo pluri ecto )
for cell in "${cell_array[@]}"
do
echo "Training ${cell} model"
start=`date +%s`
    ### train factorized model
    chrombpnet pipeline \
        -ibam ${cell}.sorted.bam \
        -d 'ATAC' \
        -g ${genome} \
        -c ${chromSizes} \
        -p input_peaks/${cell}_peaks_no_blacklist.bed \
        -n background_sets/${cell}_negatives.bed \
        -fl splits/fold_0.json \
        -b bias_models/${cell}/models/${cell}_bias.h5 \
        -o chrombpnet_models/${cell} ;
end=`date +%s`

runtime=$( echo "$end - $start" | bc -l )
echo "Model training took ${runtime} seconds"
done

## End-of-job summary
qstat -j $JOB_ID



