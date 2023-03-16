#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N metabat_s
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
mkdir Indiv_assembled_bins
conda activate metabat2
for i in Indiv_assemblies/*.contigs.fa
do
prefix=$(basename $i .contigs.fa)
metabat2 -t 40 -m 2000 -i Indiv_assemblies/${prefix}.contigs.fa -a Bowtie2_indiv_assembly_mapped/${prefix}_depth.txt -o Indiv_assembled_bins/${prefix}/${prefix}.individ.bin
done
