#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie_2
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 30
# ENTER COMMANDS HERE:
conda activate Bowtie2
mkdir bowtie2_output
for i in *_trimmed_R1.fastq.gz
do
    prefix=$(basename $i _trimmed_R1.fastq.gz)
    j=${prefix}_trimmed_R2.fastq.gz
    bowtie2 -x /isilon/projects/J-002460_LLQ/riccis/cow/cow -1 $i -2 $j -S ${prefix}.sam
--thread=$NSLOTS
done
