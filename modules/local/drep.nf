process DREP {

  label 'cpus_xxlarge'
  publishDir "$params.results/drep"

  input:
    path hq_bins_file
    path (hq_bins, stageAs: "hq_bins/*")
    path (all_bins, stageAs: "all_bins/*")
    
  output:
    path("dRep_output/*")

  script:
  """
  # Tweak the High quality bins file to be used by drep
  echo "genome,completeness,contamination,strain_heterogeneity" > header
  awk {'print \$1".fa,"\$2","\$3","0'} $hq_bins_file > corpus
  cat header corpus > checkM_results.csv
  
  # Define MPLCONFIGDIR env variable to speed up the import
  # of Matplotlib and to better support multiprocessing.
  # Otherwise, Matplotlib may create a temporary cache directory at /tmp/matplotlib-xxx
  mkdir tmp
  export MPLCONFIGDIR=$PWD/tmp

  # JSB add -centW 0 option to resolve a bug
  # https://github.com/MrOlm/drep/issues/120

  dRep dereplicate -g hq_bins/*.fa \
    -comp 90 -con 5 --processors $task.cpus \
    -strW 1 -pa 0.90 -sa 0.99 -centW 0 \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    --greedy_secondary_clustering \
    --run_tertiary_clustering \
    --genomeInfo checkM_results.csv dRep_output

  """
  
}