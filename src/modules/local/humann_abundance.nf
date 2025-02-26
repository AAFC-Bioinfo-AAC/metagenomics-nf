process HUMANN_ABUNDANCE {

  publishDir "$params.results/humann/humann_results"

  input:
    path (humann3_output, stageAs: "pathway_abundance_files/*")

  output:
    path ("Joined_pathabundance.tsv")
    path ("Joined_pathabundance_relab.tsv")
    path ("Split_humann_pathwayabundance_tables/*")

  script:
  """
  mkdir tsv
  find -L . -iname "*.tsv" | xargs -i cp {} ./tsv # cp tsv files nested in sub-folders..
  humann_join_tables -i tsv -o Joined_pathabundance.tsv --file_name pathabundance
  
  # Then normalize to relative abundance values
  humann_renorm_table \
    --input Joined_pathabundance.tsv \
    --output Joined_pathabundance_relab.tsv \
    --units relab 

  # Finally split the Joinded_pathabundance.tsv table into "stratified" and
  # "unstratified" tables.
  mkdir Split_humann_pathwayabundance_tables
  humann_split_stratified_table \
    --input Joined_pathabundance_relab.tsv \
    --output Split_humann_pathwayabundance_tables
  """
  
}