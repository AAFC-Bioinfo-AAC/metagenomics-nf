#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie_single
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
conda activate Bowtie2
for i in *.contigs.fa
do
prefix=$(basename $i .contigs.fa)
mkdir ${prefix}
bowtie2-build $i ${prefix}/${prefix}
done
