#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N checkm
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 30
conda activate checkm
for i in *.bam
do
prefix=$(basename $i _sorted.bam)
checkm lineage_wf --tab_table -t 30 -x fa Coassembled_bins/${prefix} CheckM_output_coassembled_bin
done