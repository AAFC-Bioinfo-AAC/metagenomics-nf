process METAQUAST {

  label 'cpus_xlarge'
  
  publishDir "$params.results/quast"
  
  input:
    tuple \
      val(type), \
      path(coassembly)

  output:
    path ("${type}_QUAST_coassembly/*")
    
  script:
  """
  metaquast.py ${type}/Coassembly.contigs.fa \
     --max-ref-number 0 \
    -t $task.cpus \
    -o ${type}_QUAST_coassembly
  """

}
