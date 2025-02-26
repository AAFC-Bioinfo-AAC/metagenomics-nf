process GET_BINS2 {

  label 'HQ_bins'
  publishDir "$params.results/bins"
  
  input:
      path (tsv_files, stageAs: "checkM2_hq/*")
      path (individ_assembled_bins, stageAs: "bins/*")
      
  output:
      path("High_quality_bins.txt")
      path("hq_bins/*")
      path("all_bins/*")
      
  script:
  """
  mkdir hq_bins
  mkdir all_bins
  
  cat checkM2_hq/*.tsv > High_quality_bins.txt
  
  cd bins
  
  for i in `cut -f 1 ../High_quality_bins.txt` ; do l=\$(readlink \$i.fa); ln -s \$l ../hq_bins ; done
  
  cd ..
  
  # a bit crazy recopy the links in another folder..
  cd bins
  for i in `ls *.fa`; do l=\$(readlink \$i); ln -s \$l ../all_bins; done
  """
  
}