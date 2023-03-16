#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N humann-join
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 20
conda activate humann-env
humann_join_tables -i Pathway_abundance_files -o Joined_pathabundance.tsv --file_name pathabundance
