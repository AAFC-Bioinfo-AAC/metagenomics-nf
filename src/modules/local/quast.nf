process QUAST {

  label 'cpus_xlarge'
  
  publishDir "$params.results/quast"
  
  input:
      path (dRep, stageAs: "dRep_output/*")

  output:
    path ("QUAST/*")
  
  script:
  """
  quast.py dRep_output/*.fa \
    --threads $task.cpus \
    -o QUAST
  """
  
}