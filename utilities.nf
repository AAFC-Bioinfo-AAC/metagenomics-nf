
process clean_work_files {

  cache 'lenient'

  input:
  val(file)

  output:
  val(1), emit: IS_CLEAN

  script:
  """
  $baseDir/src/clean_work_files.sh "${file}"
  """
}


