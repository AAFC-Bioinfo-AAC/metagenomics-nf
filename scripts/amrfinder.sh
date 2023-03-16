#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N amrfinder
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 30
mkdir AMRfinder_output
conda activate AMRFinder
for i in *_clean.contigs.fa
do
prefix=$(basename $i _clean.contigs.fa)
amrfinder -n $i -o AMRfinder_output/${prefix}_AMRFinder_results.txt --threads 30 --plus
done
