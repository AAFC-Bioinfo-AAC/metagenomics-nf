#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bedtools
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate bedtools
mkdir /isilon/projects/J-002460_LLQ/riccis/final_paired
for i in *.bam
do
    prefix=$(basename $i _unmapped.sorted.bam)
bedtools bamtofastq -i $i -fq /isilon/projects/J-002460_LLQ/riccis/final_paired/${prefix}_R1.fastq -fq2 /isilon/projects/J-002460_LLQ/riccis/final_paired/${prefix}_R2.fastq
done
