process BOWTIE2 {

  label 'mem_medium'
  label 'cpus_large'


  publishDir "$params.results/decontamination/genome/bowtie2"
  input:
    path genome, stageAs: "genome" //to rename the folder containing the .bt2 files
    tuple \
      val(datasetID), \
      path(trimmed_R1), \
      path(trimmed_R2), \
      path(unpaired_R1), \
      path(unpaired_R2), \
      path(html)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_unmapped.sorted.bam"), \
      path("${datasetID}.metrics.txt")

  script:
  """
  bt2_basename=\$(ls -1 genome | head -n 1 | cut -f 1 -d ".")

  bowtie2 -p $task.cpus -x "genome/\${bt2_basename}" -1 ${trimmed_R1} -2 ${trimmed_R2} | \

  samtools view -u -f 12 -F 256 --threads $task.cpus | \

  samtools sort -n -m 4G --threads $task.cpus > ${datasetID}_unmapped.sorted.bam

  # calculate the number of reads
  count4=\$(zcat ${trimmed_R1} | wc -l | cut -f 1 -d ' ')
  countR1=\$((count4/4))
  printf "${datasetID}\t\$countR1" >> ${datasetID}.metrics.txt

  count4=\$(zcat ${trimmed_R2} | wc -l | cut -f 1 -d ' ')
  countR2=\$((count4/4))
  printf "\t\$countR2" >> ${datasetID}.metrics.txt

  # Add the number of unmmaped_reads in a text file
  unmapped=\$(samtools view -c ${datasetID}_unmapped.sorted.bam)
  printf "\t\$unmapped" >> ${datasetID}.metrics.txt

  # bash does not natively support floating point arithmetics
   
  ratio=\$(echo "print(1-(\$unmapped / ( \$countR1 + \$countR2)))" | python3)

  printf "\t\$ratio\n" >> ${datasetID}.metrics.txt
  """
}