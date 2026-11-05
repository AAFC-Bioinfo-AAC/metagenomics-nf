process RENAME_SEQUENCES {

  publishDir "$params.results/renamed_input_data"

  input:
    path (seq, stageAs: "data/*")
    path map_file
    
  output:
    tuple \
      path("log.txt"), \
      path("renamed/*")
  
  script:
  """
  mkdir renamed
  $baseDir/src/rename_sequences.py ${map_file} "./data" "renamed/" "fastq.gz" > log.txt
  """
  
}
