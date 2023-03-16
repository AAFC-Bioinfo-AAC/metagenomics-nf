#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate pigz
pigz *fastq -p 2
