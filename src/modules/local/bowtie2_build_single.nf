process BOWTIE2_BUILD_SINGLE {

  label 'cpus_xlarge'
  publishDir "$params.results/indiv_assemblies/bwt2_index"

  input:
    tuple \
      val(datasetID), \
      path(megahit_individual_outfiles)

  output:
    tuple \
      val(datasetID), \
      path ("bwt2_index/*")
  
  script:
  """
  mkdir bwt2_index
  bowtie2-build \
  ${datasetID}/${datasetID}.contigs.fa \
    bwt2_index/${datasetID} \
    --quiet \
    --threads $task.cpus
  """
  
}