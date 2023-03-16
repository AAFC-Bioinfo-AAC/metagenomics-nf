#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N metabat_bins
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 4
conda activate metabat2
mkdir Coassembled_bins
for i in *_sorted.bam
do
prefix=$(basename $i _sorted.bam)
metabat2 -i /isilon/projects/J-002460_LLQ/riccis/final_paired/Megahit_coassembly/Coassembly.contigs.fa -a /isilon/projects/J-002460_LLQ/riccis/final_paired/bowtie2_output_for_metabat/Coassembly_depth.txt -o Coassembled_bins/${prefix}/${prefix}.bin -t 4 -m 1500 -v
done
