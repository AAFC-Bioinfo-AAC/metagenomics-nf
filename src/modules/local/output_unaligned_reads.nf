process OUTPUT_UNALIGNED_READS {

  publishDir "$params.results/prepared_reads"

  input: 
    tuple \
      val(datasetID), \
      path(aln), \
      path(txt)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_R1.fastq.gz"), \
      path("${datasetID}_R2.fastq.gz")
  
  script:
  """
  bedtools bamtofastq -i ${aln} -fq ${datasetID}_R1.fastq -fq2 ${datasetID}_R2.fastq &&
  gzip ${datasetID}_R1.fastq &&
  gzip ${datasetID}_R2.fastq
  """
  
}