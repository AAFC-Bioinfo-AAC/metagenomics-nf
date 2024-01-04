process JGI_SUMMARIZE_SINGLE {

  publishDir "$params.results/indiv_assemblies/jgi"

  input: 
    tuple \
      val(datasetID), \
      path(aln)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt")
  
  script:
  """
  jgi_summarize_bam_contig_depths \
    --outputDepth ${datasetID}_depth.txt \
    $aln
  """
  
}