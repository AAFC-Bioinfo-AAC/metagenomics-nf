process KAIJU_FULL_TAX_TABLE {

  label 'mem_medium'
  
  publishDir "$params.results/kaiju/kaiju_full_tax_table"

  input:
    path db
    tuple \
      val(datasetID), \
      path(kaiju_out)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.all_tax.summary.tsv")
  
  script:
  """
  kaiju2table -t ${db}/nodes.dmp \
        -n ${db}/names.dmp \
        -r species \
        -l superkingdom,phylum,class,order,family,genus,species \
        -o ${datasetID}.all_tax.summary.tsv \
        ${kaiju_out}
  """
  
}