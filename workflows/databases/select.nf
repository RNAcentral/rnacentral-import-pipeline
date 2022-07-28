nextflow.enable.dsl=2



process check_db_md5 {
  container ''

  input:
    tuple val(db_name), val(remote)

  output:
    path("*.csv")


    """
    wget -O target_file $remote
    echo -n "$db_name," >> latest_md5s.csv && md5 -q target_file >> latest_md5s.csv
    """
}


process make_selection {
  publishDir "$projectDir"

  input:
    path latest_md5s

  output:
    path ("*.config")

  """
  env
  $projectDir/../../bin/rnac scan-imports select-for-import $latest_md5s
  """

}


workflow select {

  Channel.fromPath(params.import_selection_remotes) \
  | splitCsv
  | map { row -> tuple(row[0], row[1])}
  | check_db_md5
  | collectFile
  | make_selection

}


workflow {
  select()
}
