process MERGE_FULL_TAX_FILES {

  publishDir "$params.results/kaiju/kaiju_merged"

  input:
    path (tsv_files, stageAs: "kaiju_full_tax_table/*")   // to put input files in a folder parsed by R script

  output:
    path("kaiju_merged_full_tax.csv")
  
  script:
  """
  $baseDir/src/merge_full_tax_files.R
  """
  
}