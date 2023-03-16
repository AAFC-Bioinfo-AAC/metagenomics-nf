#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N bowtie-assembly
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 2
conda activate Bowtie2
mkdir coassembly
bowtie2-build Megahit_coassembly/Coassembly.contigs.fa coassembly/coassembly
