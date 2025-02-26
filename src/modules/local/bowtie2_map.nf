process BOWTIE2_MAP {

  label 'cpus_large'
  
  publishDir "$params.results/coassemblies/bwt2_output_for_metabat"

    input:
    tuple \
    val(type), \
    path(final_R1), \
    path(final_R2), \
    val(datasetID), \
    path(bwt2_index, stageAs: "coassembly/*") //to rename the folder containing the bt2 files


    output:
    tuple \
        val(datasetID), \
        val(type), \
        path("${datasetID}_sorted.bam")

    script:
    """
    bowtie2 -x coassembly/coassembly \
            -1 ${final_R1} \
            -2 ${final_R2} \
            -p $task.cpus | \
            samtools sort \
            --threads $task.cpus \
            -O BAM -o ${datasetID}_sorted.bam

    """

}