process METABAT2_BIN_SINGLE {

  label 'cpus_xlarge'

  publishDir "$params.results/indiv_assemblies/metabat2_bins"
  
  input:
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt"), \
      path(megahit_individual_outfiles, stageAs: "megahit/*")
      
  output:
    tuple \
      val(datasetID), \
      path("${datasetID}/*"), optional: true

  script:
  """
  mkdir ${datasetID}
  metabat2 \
    -i megahit/${datasetID}/${datasetID}.contigs.fa \
    -a ${datasetID}_depth.txt \
    -o ${datasetID}/${datasetID}.individ.bin \
    -t $task.cpus \
    -m 2000 \
    -v
  """
  
}