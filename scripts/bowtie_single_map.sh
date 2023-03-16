#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie_single_map
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 40
conda activate Bowtie2
mkdir Bowtie2_indiv_assembly_mapped
for i in *clean_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
j=${prefix}_R2.fastq.gz
bowtie2 -x Indiv_assemblies/${prefix}/${prefix} -1 $i -2 $j -S Bowtie2_indiv_assembly_mapped/${prefix}.sam -p 40
done
