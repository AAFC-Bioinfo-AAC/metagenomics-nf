#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N copy_paths
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 30
conda activate humann-env
mkdir Pathway_abundance_files
for i in *_cat.fastq.gz
do
prefix=$(basename $i _cat.fastq.gz)
cp ${prefix}_humann3_output/${prefix}_cat_pathabundance.tsv Pathway_abundance_files/
done
