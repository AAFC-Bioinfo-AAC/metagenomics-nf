process QUALITY_FILTERING {
    
  label 'cpus_medium'

  publishDir "$params.results/trimmed_reads"

  input: 
    tuple val(datasetID), path(read1), path(read2)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_trimmed_R1.fastq.gz"), \
      path("${datasetID}_trimmed_R2.fastq.gz"), \
      path("${datasetID}_unpaired_R1.fastq.gz"), \
      path("${datasetID}_unpaired_R2.fastq.gz"), \
      path("${datasetID}.html")
  
  script:
  """
    fastp -i $read1 \
        -I $read2 \
        -o ${datasetID}_trimmed_R1.fastq.gz \
        -O ${datasetID}_trimmed_R2.fastq.gz \
        --unpaired1 ${datasetID}_unpaired_R1.fastq.gz \
        --unpaired2 ${datasetID}_unpaired_R2.fastq.gz \
        --length_required 100 \
        --thread=$task.cpus \
        --cut_tail --cut_front \
        --cut_mean_quality 15 \
        --qualified_quality_phred 15 \
        --cut_window_size 4 \
        --detect_adapter_for_pe \
        --report_title="${datasetID}" \
        --html=${datasetID}.html
  """
  
}