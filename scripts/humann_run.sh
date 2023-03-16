#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N humann
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mamba activate humann-env
for i in *_cat.fastq.gz
do
prefix=$(basename $i _cat.fastq.gz)
humann -i $i -o ${prefix}_humann3_output --threads 40 --remove-temp-output --nucleotide-database /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/cleaned_all/clean_fastq/chocophlan --protein-database /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/cleaned_all/clean_fastq/uniref
done