process BOWTIE2_MAP_SINGLE {

  label 'cpus_large'

  publishDir "$params.results/indiv_assemblies/bwt2_output_for_metabat"

  input:
    tuple \
      val(datasetID), \
      path ("bwt2_index/*", stageAs: "index/*"), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_sorted.bam")

  script:
  """
  bowtie2 -x index/${datasetID} \
          -1 ${final_R1} \
          -2 ${final_R2} \
          -p $task.cpus | \
          samtools sort \
          --threads $task.cpus \
          -O BAM -o ${datasetID}_sorted.bam
          
  """
  
}