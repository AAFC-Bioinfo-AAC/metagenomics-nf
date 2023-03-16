#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate Bowtie2
mkdir cow_2
bowtie2-build /isilon/projects/J-002460_LLQ/riccis/Cow_PhiX_genomes.fna /isilon/projects/J-002460_LLQ/riccis/cow_2/cow
