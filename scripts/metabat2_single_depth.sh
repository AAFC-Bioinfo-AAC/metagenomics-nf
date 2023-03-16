#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N metabat_s_d
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
conda activate metabat2
for i in *_sorted.bam
do
prefix=$(basename $i _sorted.bam)
jgi_summarize_bam_contig_depths --outputDepth ${prefix}_depth.txt $i
done
