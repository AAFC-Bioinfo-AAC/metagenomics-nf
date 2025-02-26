process JGI_SUMMARIZE {
  
  publishDir "$params.results/coassemblies/jgi"

  input: 
    tuple \
      val(datasetID), \
      val(type), \
      path(aln)
 
  output:
    tuple \
      val(type), \
      val(datasetID), \
      path("${datasetID}_depth.txt")
  
  script:
  """
  jgi_summarize_bam_contig_depths \
    --outputDepth ${datasetID}_depth.txt \
    $aln
  """
  
}