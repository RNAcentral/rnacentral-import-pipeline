process fetch_data {
  when { !params.databases.plncdb.prefetch and params.databases.plncdb.run }

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
  rnac notify step "Data parsing for PLncDB" $params.databases.plncdb.data_path$data
  rnac plncdb parse $params.databases.plncdb.data_path$data
  """
}

workflow plncdb {
  emit: data_files

  main:

  Channel.fromPath("$params.databases.plncdb.data_path/*", type:'dir') \
  | parse_data
  | collect
  | set { data_files }


}
