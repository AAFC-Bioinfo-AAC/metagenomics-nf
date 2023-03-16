#$ -S /bin/bash
#$ -V
#$ -N fastp 
#$ -cwd
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -pe smp 20

#Create folders
mkdir paired_trimmed
mkdir unpaired_trimmed

# load conda env 
conda activate fastp

# test_AB2020-05_R1.fastq.gz
for i in /isilon/projects/J-002460_LLQ/data/*_R1.fastq.gz
do
   prefix=$(basename $i _R1.fastq.gz)
   j=${prefix}_R2.fastq.gz

  fastp -i $i -I /isilon/projects/J-002460_LLQ/data/$j -o paired_trimmed/${prefix}_trimmed_R1.fastq.gz -O paired_trimmed/${prefix}_trimmed_R2.fastq.gz --unpaired1 unpaired_trimmed/${prefix}_unpaired_R1.fastq.gz --unpaired2 unpaired_trimmed/${prefix}_unpaired_R2.fastq.gz --length_required 100 --cut_tail --cut_front --cut_mean_quality 15 --qualified_quality_phred 15 cut_window_size 4 --detect_adapter_for_pe --report_title="${prefix}" --html=${prefix}.html --thread=$NSLOTS
done
