#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N kaiju-tax-reads
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mkdir Kaiju_output_reads
conda activate kaiju
for i in *clean_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
kaiju -t /isilon/common/reference/databases/kaiju_db/latest/nodes.dmp -f /isilon/common/reference/databases/kaiju_db/latest/kaiju_db_nr_euk.fmi -i $i -j $j -v -o Kaiju_output_reads/${prefix}.out -z 40
done
