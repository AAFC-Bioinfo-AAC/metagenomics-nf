#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N sort
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate samtools
mkdir Sorted_BAM
for i in *.bam
do
prefix=$(basename $i .bam)
samtools sort -@ 2 -n $i > Sorted_BAM/${prefix}.sorted.bam
done
