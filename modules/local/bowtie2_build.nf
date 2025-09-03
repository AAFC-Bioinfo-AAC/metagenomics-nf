process BOWTIE2_BUILD {

  publishDir "$params.results/coassemblies/bwt2_index/"

  input:
    tuple \
      val(type), \
      path (megahit_coassembly_outfiles)

  output:
    tuple \
      val(type), \
      path ("${type}_index/*")
  
  script:
  """
  mkdir ${type}_index
  bowtie2-build ${type}/Coassembly.contigs.fa ${type}_index/coassembly
  """
  
}
