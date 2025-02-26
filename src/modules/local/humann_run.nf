process HUMANN_RUN {

  label 'mem_medium'
  label 'cpus_xxlarge'

  publishDir "$params.results/humann/humann_run"
  input:
    path chocophlan_db
    path metaphlan_db
    path uniref_db
    tuple \
      val(datasetID), \
      path(reads)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_humann3_output")

  script:
  """
  # Sometimes humann3 is unable to write in the output
  # e.g. :
  # 40K drwxrwx---. 2 jsbrouard grp_jsbrouard   22 25 aoû 13:10 G2-E1-13_humann3_output
  # Try to create the output folder before humann3 and change permissions
  mkdir -p ${datasetID}_humann3_output &&
  chmod a+rwx ${datasetID}_humann3_output &&

  humann -i $reads -o ${datasetID}_humann3_output \
         --threads $task.cpus \
         --remove-temp-output \
         --metaphlan-options "--bowtie2db ${metaphlan_db} --index mpa_vJan21_CHOCOPhlAnSGB_202103" \
         --nucleotide-database $chocophlan_db \
         --protein-database $uniref_db
  """
  
}