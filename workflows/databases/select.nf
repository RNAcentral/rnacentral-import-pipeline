nextflow.enable.dsl=2



process check_db_md5 {
  container ''

  input:
    tuple val(db_name), val(remote)

  output:
    path("*.csv")


    """
    wget -O target_file $remote
    echo -n "$db_name," >> latest_md5s.csv && md5sum target_file | awk 'BEGIN {fs="[ ]"}; {print \$1}' >> latest_md5s.csv
    """
}


process make_selection {
  publishDir "$projectDir"

  input:
    path latest_md5s

  output:
    path ("*.config")
    path ("$latest_md5s")

  """
  rnac scan-imports select-for-import $latest_md5s
  """

}


process update_tracker_table {
  input:
    path latest_md5s

  """
  rnac scan-imports update-tracker $latest_md5s
  """
}


workflow select {

  Channel.fromPath(params.import_selection_remotes) \
  | splitCsv
  | map { row -> tuple(row[0], row[1])}
  | check_db_md5
  | collectFile
  | ( make_selection & update_tracker_table )

}


workflow {
  select()
}
