#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N human-dwnl
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mamba activate humann-env
humann_databases --download chocophlan full  /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/cleaned_all/clean_fastq/ --update-config yes
humann_databases --download uniref uniref90_diamond /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/cleaned_all/clean_fastq/
humann_databases --download utility_mapping full /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/cleaned_all/clean_fastq/ --update-config yes