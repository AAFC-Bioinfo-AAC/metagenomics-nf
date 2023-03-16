#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N checkm_s
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mkdir CheckM_output
conda activate checkm
for i in Indiv_assemblies/*.contigs.fa
do
prefix=$(basename $i .contigs.fa)
checkm lineage_wf --tab_table -t 40 -x fa Indiv_assembled_bins/${prefix} CheckM_output > CheckM_individual_bin_summaries
done
