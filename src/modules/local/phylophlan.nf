process PHYLOPHLAN {

 label 'cpus_xxlarge'
 
 publishDir "$params.results/phylophlan"

  input:
      path db
      path(dereplicated_genomes, stageAs: "dRep_output/*")
      
 
  output:   
      path("Phylophlan_output/*")
  
  script:
    """
    mkdir -p phylophlan_databases &&
    chmod a+rwx phylophlan_databases &&

    # for generating the four default configuration files
    phylophlan_write_default_configs.sh &&
  
    phylophlan -d phylophlan \
                 -i dRep_output \
                 -o Phylophlan_output \
                 --db_type a \
                -f supermatrix_aa.cfg \
                --nproc $task.cpus \
                --diversity low \
                --fast \
                --verbose \
                --genome_extension fa \
                --databases_folder $db
    """
}