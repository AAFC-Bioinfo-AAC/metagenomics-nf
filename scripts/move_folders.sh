#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N move
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 4
for i in *_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
mv ${prefix}/${prefix}.contigs.fa Indiv_assemblies/
rm -r ${prefix}
done