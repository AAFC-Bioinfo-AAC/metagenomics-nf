#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N samtools
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 4
mkdir Contigs_mapped
conda activate samtools
for i in *.sam
do
prefix=$(basename $i .sam)
samtools sort $i > Contigs_mapped/${prefix}_sorted.bam
done
