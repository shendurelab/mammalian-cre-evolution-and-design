#$ -S /bin/bash
#$ -cwd    ## use current working directory
#$ -N sm_EB_bamMaking
#$ -j yes  ## merge stdout and stderr
#$ -l mfree=8G,h_rt=168:00:00

conda activate snakemake

# can switch to -P sage if shendure going slow, or add -now no
snakemake --keep-going -j 500 --cluster "qsub -N {params.job_name} -o {params.error_out_file}.out -e {params.error_out_file}.error -l h_rt={params.run_time} -pe serial {params.cores} -l h_vmem={params.memory}G"

echo "Done"
