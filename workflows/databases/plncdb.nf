nextflow.enable.dsl = 2

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

  queue 'short'
  memory { 8.GB * task.attempt }

  errorStrategy 'retry'
  maxRetries 16

  input:
  path data

  output:
  path('*.csv')

  """
  # rnac notify step "Data parsing for PLncDB" $params.databases.plncdb.data_path$data
  rnac plncdb parse $params.databases.plncdb.data_path$data
  """
}

workflow plncdb {
  emit: data_files

  main:
  if( params.databases.plncdb.run ) {
    Channel.fromPath("$params.databases.plncdb.data_path/*", type:'dir') \
    | parse_data \
    | flatten
    | collectFile() {csvfile -> [csvfile.name, csvfile.text]} \
    | set { data_files }
  }
  else {
  Channel.empty() | set { data_files }
  }


}
