process KAIJU {

  label 'mem_xlarge'
  label 'cpus_large'
  
  publishDir "$params.results/kaiju/kaiju_1"

  input:
    path db
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.out")
  
  script:
  """
  kaiju -t ${db}/nodes.dmp \
        -f ${db}/kaiju_db_nr_euk.fmi \
        -i ${final_R1} \
        -j ${final_R2} \
        -v \
        -o ${datasetID}.out \
        -z $task.cpus
  """
  
}
