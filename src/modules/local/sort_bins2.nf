process SORT_BINS2 {

  label 'bins'
  publishDir "$params.results/sorted_bins"
  
  input:
    tuple \
      val(datasetID), \
      path (checkm, stageAs: "checkm2_out/*")
      
  output:   
      path("${datasetID}_checkM2_i_hq.tsv")
 
  script:
  """
  awk '{if (\$2 > 90 && \$3 < 5) {print}}' checkm2_out/quality_report.tsv >  ${datasetID}_checkM2_i_hq.tsv
  """
  
}