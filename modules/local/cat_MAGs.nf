process CAT_MAGs {

  label 'cpus_xxlarge'
  publishDir "$params.results/all_mags"

  input:
    tuple \
      val(datasetID), \
      path(hq_bins_from_individual_assemblies, stageAs: "hq_bins/*"), \
      path(all_bins_from_individual_assemblies, stageAs: "all_bins/*"), \
      path(check_m_i)

    tuple \
      val(datasetID), \
      path(hq_bins_from_coassemblies, stageAs: "hq_bins/*"), \
      path(all_bins_from_from_coassemblies, stageAs: "all_bins/*"), \
      path(check_m_c)
 

  output:
     tuple \
      val(datasetID), \
      path("hq_bins"), \
      path("all_bins")

  script:
  """
  ls -l hq_bins >> tester
  pet


  """
  
}
