process BOWTIE2 {

  label 'mem_medium'
  label 'cpus_large'
  
  publishDir "$params.results/decontamination/genome/bowtie2"
  input:
    path genome, stageAs: "genome" //to rename the folder containing the bt2 files
    val(genome_basename)
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
      path("${datasetID}_unmapped.sorted.bam")
  
  script:
  """
  bowtie2 -p $task.cpus -x "genome/${genome_basename}" -1 ${trimmed_R1} -2 ${trimmed_R2} | \
  
  samtools view -u -f 12 -F 256 --threads $task.cpus | \

  samtools sort -n -m 4G --threads $task.cpus | \

  # step3 (old SORTBAM process)
  samtools sort -n -m 4G --threads $task.cpus > ${datasetID}_unmapped.sorted.bam
  """
  
}