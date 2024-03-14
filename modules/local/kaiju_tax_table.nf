process KAIJU_TAX_TABLE {

  label 'mem_medium'
  
  publishDir "$params.results/kaiju/kaiju_tax_table"

  input:
    path db
    tuple \
      val(datasetID), \
      path(kaiju_out)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.species.summary.tsv")
  
  script:
  """
  kaiju2table -t ${db}/nodes.dmp \
        -n ${db}/names.dmp \
        -r species \
        -o ${datasetID}.species.summary.tsv \
        ${kaiju_out}
  """
  
}