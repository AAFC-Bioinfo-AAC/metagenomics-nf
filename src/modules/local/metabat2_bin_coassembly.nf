process METABAT2_BIN_COASSEMBLY {

  label 'cpus_xlarge'
  
  publishDir "$params.results/coassemblies/metabat2_bins"
  
  input:
    tuple \
      val(type), \
      path (megahit_coassembly_outfiles, stageAs: "megahit/*"), \
      val(datasetID), \
      path("${datasetID}_depth.txt")
 
  output:
    tuple \
      val(datasetID), \
      path("${datasetID}/*"), optional: true
  
  script:
  """
  mkdir ${datasetID}
  metabat2 \
    -i megahit/Coassembly.contigs.fa \
    -a ${datasetID}_depth.txt \
    -o ${datasetID}/${datasetID}.bin \
    -t $task.cpus \
    -m 2000 \
    -v
  """
  
}