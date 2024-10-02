/*
According to the Megahit documentation, options --min-count 2 and
--k-list 21,41,61,81,99 are the default (generic metagenomes)
Here we used a list of 4 kmers as in a this paper :
Vosloo S, Huo L, Anderson CL, Dai Z,
Sevillano M, Pinto A. 2021. Evaluating de novo assembly and binning
strategies for time series drinking water metagenomes. Microbiol Spectr9:e01434-21.
https://doi.org/10.1128/Spectrum.01434-21
*/

process COASSEMBLY {

  label 'mem_xxxlarge'
  label 'cpus_xxxxlarge'
  
  publishDir "$params.results/coassemblies/coassemblies"

  input:
    tuple \
      val(type), \
      path (readsR1, stageAs: "readsR1/*"), \
      path (readsR2, stageAs: "readsR2/*"), \
      val(datasetID)

  output:
    tuple \
    val (type), \
    path ("${type}_Megahit_coassembly/*")

  script: 
  """
  cat $readsR1 > coassembly_R1.fastq.gz
  cat $readsR2 > coassembly_R2.fastq.gz
  megahit -1 coassembly_R1.fastq.gz \
          -2 coassembly_R2.fastq.gz \
          -o ${type}_Megahit_coassembly \
          --out-prefix Coassembly \
          -t $task.cpus \
          --min-contig-len 1000 \
          --min-count 2 \
          --k-list 21,41,61,81,99 &&
  
  # Remove intermediate files (intermediate contigs)
  rm -rf ${type}_Megahit_coassembly/intermediate_contigs
  """

}
