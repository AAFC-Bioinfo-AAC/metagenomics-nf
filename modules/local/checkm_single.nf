process CHECKM_SINGLE {

  label 'cpus_xxlarge'
  
  publishDir "$params.results/indiv_assemblies/checkM2_output"
  errorStrategy 'ignore' //sometimes DIAMOND produces no annotation and Checkm2 Exit
  
  input:
    path db
    tuple \
      val(datasetID), \
      path (metabat2_individ_outfiles, stageAs: "indiv_assembled_bins/*")

  output:
      tuple \
        val(datasetID), \
        path("$datasetID/*")
 
  script:
  """
  export HDF5_USE_FILE_LOCKING='FALSE'
  mkdir ${datasetID}
  checkm2 predict \
    --database_path $db \
    --threads 20 \
    -x fa \
    --input indiv_assembled_bins \
    --output-directory ${datasetID}
  """
  
}