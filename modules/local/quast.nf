process QUAST {

  label 'cpus_xlarge'
  
  publishDir "$params.results/quast"
  
  input:
    path(dereplicated_genomes, stageAs: "dRep_output/*")

  output:
    path "QUAST_replicated_MAGs/*"
  
  script:
  """
  quast.py dRep_output/dereplicated_genomes/*.fa \
    --threads $task.cpus \
    -o QUAST_replicated_MAGs
  """
  
}