process fetch_data {
  when { params.databases.plncdb.run }

  executor 'local'

  output:
  path("data")

  """
  rnac plncdb fetch $params.databases.plncdb.remote data
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
    fetch_data | parse_data | set { data_files }

}
