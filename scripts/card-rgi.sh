
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N rgi
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mkdir Card-rgi_output
conda activate rgi
for i in *clean_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
rgi bwt --read_one $i --read_two $j --aligner bowtie2 --output_file Card-rgi_output/${prefix}_CARD_RGI --clean --threads 40 --local
done
