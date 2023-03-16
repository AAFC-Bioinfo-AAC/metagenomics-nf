#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N humann-split
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 20
conda activate humann-env
mkdir Split_humann_pathwayabundance_tables
humann_split_stratified_table --input Joined_pathabundance_relab.tsv --output Split_humann_pathwayabundance_tables
