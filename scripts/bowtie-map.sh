#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie-map
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 4
mkdir bowtie2_output_for_metabat
conda activate Bowtie2
for i in *_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
bowtie2 -x coassembly/coassembly -1 $i -2 $j -S bowtie2_output_for_metabat/${prefix}.sam -p 4
done
