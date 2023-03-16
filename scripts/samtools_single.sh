#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N samtools_s
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
conda activate samtools
for i in *.sam
do
prefix=$(basename $i .sam)
samtools sort $i > ${prefix}_sorted.bam
done
