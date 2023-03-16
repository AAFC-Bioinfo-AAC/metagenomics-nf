#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N tax-table_reads
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 20
conda activate kaiju
for i in *clean.out
do
prefix=$(basename $i .out)
kaiju2table -t /isilon/common/reference/databases/kaiju_db/latest/nodes.dmp -n /isilon/common/reference/databases/kaiju_db/latest/names.dmp -r family -o ${prefix}.family.summary.tsv $i
done
