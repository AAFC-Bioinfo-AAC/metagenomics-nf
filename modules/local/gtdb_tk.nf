process GTDB_TK {

  label 'cpus_xxlarge'
  label 'mem_xxlarge'

  publishDir "$params.results/GTDB"

  input:
    path(db)
    path(dereplicated_genomes, stageAs: "dRep_output/*")
 
  output:
    path("GTDBtk_output/*")
  
  
  script:
  """
  export GTDBTK_DATA_PATH=$db
  
  gtdbtk classify_wf \
         --genome_dir \
         dRep_output/dereplicated_genomes \
         -x fa \
         --out_dir GTDBtk_output \
         --cpus $task.cpus
  """
  
}