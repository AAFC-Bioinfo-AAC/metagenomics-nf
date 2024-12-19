process KEEP_HQ_BINS_2 {

  label 'bins'
  publishDir "$params.results/sorted_bins"
  
  input:
    tuple \
      val(datasetID), \
      path (metabat2_individ_outfiles, stageAs: "indiv_assembled_bins/*"), \
      path (checkm, stageAs: "checkm2_out/*")
      
  output:
     tuple \
      val(datasetID), \
      path("hq_bins"), \
      path("checkM2_i_hq.tsv")
 
  script:
  """

  mkdir hq_bins
  
  awk '{if (\$2 > 90 && \$3 < 5) {print}}' checkm2_out/quality_report.tsv >  checkM2_i_hq.tsv

  
  for i in `cut -f 1 checkM2_i_hq.tsv` ; do l=\$(readlink \$i.fa); ln -s \$l ./hq_bins ; done

  """
  
}