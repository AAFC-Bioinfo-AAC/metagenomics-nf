#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N megahit_single
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mkdir Indiv_assemblies
conda activate megahit
for i in *_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
megahit -1 $i -2 $j -o ${prefix} --out-prefix ${prefix} -t 40 --min-contig-len 1000
mv ${prefix}/${prefix}.contigs.fa Indiv_assemblies/
rm -r ${prefix}
done
