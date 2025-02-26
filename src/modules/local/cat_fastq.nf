process CAT_FASTQ {

  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_cat.fastq.gz")
    
  script:
  """
  cat ${final_R1} ${final_R2} > ${datasetID}_cat.fastq.gz
  """
  
}