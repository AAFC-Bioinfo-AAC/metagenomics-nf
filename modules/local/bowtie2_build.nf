process BOWTIE2_BUILD {

  publishDir "$params.results/coassemblies/bwt2_index/"

  input:
    tuple \
      val(type), \
      path (megahit_coassembly_outfiles, stageAs: "megahit/*")

  output:
    tuple \
      val(type), \
      path ("${type}/*")
  
  script:
  """
  mkdir ${type}
  bowtie2-build megahit/Coassembly.contigs.fa ${type}/coassembly
  """
  
}