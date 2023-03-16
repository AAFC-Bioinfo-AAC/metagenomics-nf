#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N megahit
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate megahit
 megahit -1 test_AB2020-05_R1.fastq.gz -2 test_AB2020-05_R2.fastq.gz -o Megahit_coassembly_2 --out-prefix Coassembly -t 2 --min-contig-len 1000
