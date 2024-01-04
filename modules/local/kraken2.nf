process KRAKEN2 {

  label 'cpus_xlarge'
  label 'mem_xlarge'
  publishDir "$params.results/kraken2"

  input:
    val confidence_threshold
    path db
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)
 
  output:   
    tuple \
      val(datasetID), \
      path("Kraken2_${datasetID}.report.txt")

  script:
  """
  kraken2 --use-names \
  --threads $task.cpus \
  --db $db \
  --paired ${final_R1} ${final_R2} \
  --report Kraken2_${datasetID}.report.txt \
  --confidence ${confidence_threshold} \
  --report-zero-counts > /dev/null
  # It is important to redirect the large Kraken2 output to /dev/null
  # Otherwise, massive info is written in .command.log
  # and .command.out Nextflow files
  """
  
}