process COVERM {

 label 'cpus_xlarge'
 label 'mem_large'
 publishDir "$params.results/coverM"

  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)
    path(dereplicated_genomes, stageAs: "dRep_output/*")
 
  output:
     tuple \
       val(datasetID), \
       path("${datasetID}_coverM_output.txt")
  
  script:
  """
  # Define TMPDIR env variable to avoid Samtools (used by coverm) to write in /tmp in local nodes!
  mkdir -p ./tmp
  export TMPDIR=\$PWD/tmp
    
  coverm genome -1 ${final_R1} \
                -2 ${final_R2} \
                --genome-fasta-directory dRep_output \
                --genome-fasta-extension fa \
                --min-covered-fraction 1 \
                --threads $task.cpus -v \
                --output-file ${datasetID}_coverM_output.txt
  """
}
