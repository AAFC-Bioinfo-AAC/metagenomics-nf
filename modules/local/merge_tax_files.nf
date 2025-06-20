process MERGE_TAX_FILES {

  publishDir "$params.results/kaiju/kaiju_merged"

  input:
    path (tsv_files, stageAs: "Species/*")   // to put input files in a folder parsed by R script

  output:
    path("kaiju_merged_species.csv")
  
  script:
  """
  $baseDir/src/merge_tax_files.R
  """
  
}