process fetch_and_process {
  tag { "$name" }
  input:
  tuple val(name), path(data), val(column_name)

  output:
  path('*.csv')

  script:
  """
  rnac lookup genecards $column_name $data urs-info.pickle
  rnac external $name $data urs-info.pickle .
  """
}

workflow genecards_suite {
  emit: data
  main:
    Channel.fromList([
      ['genecards', params.databases.genecards.remote, params.databases.genecards.column],
      ['malacards', params.databases.malacards.remote, params.databases.malacards.column],
    ]) \
    | filter { name, r, c -> params.databases[name].run } \
    | fetch_and_process \
    | set { data }
}
