process CHECKM {

  label 'cpus_xxlarge'
  errorStrategy 'ignore' //sometimes DIAMOND produces no annotation and Checkm2 Exit

  publishDir "$params.results/coassemblies/checkM2_output"
  
  input:
    tuple \
      val(datasetID), \
      path (metabat2_coassembly_outfiles, stageAs: "Coassembled_bins/*")

  output:
    tuple \
      val(datasetID), \
      path("$datasetID/*")
  
  script:
  """
  export HDF5_USE_FILE_LOCKING='FALSE'
  checkm2 predict \
    --database_path $params.checkm2_db \
    --threads 20 \
    -x fa \
    --input Coassembled_bins \
    --output_directory ${datasetID}
  """
  
}
