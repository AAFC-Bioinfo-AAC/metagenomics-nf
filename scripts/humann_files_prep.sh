#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N concatenate
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
for i in *clean_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
cat $i $j >${prefix}_cat.fastq.gz
done
