process RENAME_SEQUENCES {

  publishDir "$baseDir/data"

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
  $baseDir/src/rename_sequences.py $map_file \$PWD/data renamed/ fastq.gz > log.txt
  """
  
}