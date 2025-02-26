process KEEP_HQ_BINS_2 {

  label 'bins'
  publishDir "$params.results/sorted_bins_from_individual_assemblies"
  
  input:
    tuple \
      val(datasetID), \
      path (metabat2_individ_outfiles, stageAs: "bins/*"), \
      path (checkm, stageAs: "checkm2_out/*")
      
  output:
     tuple \
      val(datasetID), \
      path("hq_bins/*"), \
      path("checkM2_i_hq.tsv")
 
  script:
  """
  mkdir hq_bins
  mkdir all_bins

  # Trick to avoid a missing output file error
  touch hq_bins/${datasetID}.void  

  awk '{if (\$2 > 90 && \$3 < 5) {print}}' checkm2_out/quality_report.tsv >  checkM2_i_hq.tsv

  cd bins
  
  # do the hq bins
  for i in `cut -f 1 ../checkM2_i_hq.tsv` ; do l=\$(readlink \$i.fa); ln -s \$l ../hq_bins ; done

  # do all bins
  for i in `ls *.fa`; do l=\$(readlink \$i); ln -s \$l ../all_bins; done
  """
  
}
