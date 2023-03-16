#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N metabat_depth
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 4
conda activate metabat2
jgi_summarize_bam_contig_depths --outputDepth Coassembly_depth.txt Contigs_mapped/*.bam
