process fetch_data {
  when { params.databases.plncdb.run }

  containerOptions "--contain --bind $baseDir"

  output:
  path("data")

  """
  rnac plncdb fetch-data $params.databases.plncdb.urls data
  """
}

process parse_data {
  when { params.databases.plncdb.run }

  input:
  path data

  output:
  path ('*.csv')

  """
  rnac plncdb parse $data
  """
}

workflow plncdb {
  emit: data_files

  main:
  // If data is prefetched, this should use the prefetched data directory
    (params.databases.plncdb.prefetch ? params.databases.plncdb.prefetch : fetch_data) | parse_data | set { data_files }

}
