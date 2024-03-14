process MEGAHIT_SINGLE {
  
  label 'mem_xlarge'
  label 'cpus_xlarge'
  
  publishDir "$params.results/indiv_assemblies/megahit"

  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path ("${datasetID}")

  script:
  """
  megahit -1 ${final_R1} \
          -2 ${final_R2} \
          -o ${datasetID} \
          --out-prefix ${datasetID} \
          -t $task.cpus \
          --min-contig-len 1000 &&

  # Remove intermediate files (intermediate contigs)
  rm -rf ${datasetID}/intermediate_contigs
  """
  
}