#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N humann-norm
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 20
conda activate humann-env
humann_renorm_table --input Joined_pathabundance.tsv --output Joined_pathabundance_relab.tsv --units relab
