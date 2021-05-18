process fetch {
  tag { "$name" }
  input:
  tuple val(name), path(data), val(column_name)

  output:
  tuple val(name), path('data.tsv'), val(column_name)

  script:
  if (data.getExtension() == "gz")
    """
    zcat $data > data.tsv
    """
  else
    """
    cp $data data.tsv
    """

}

process process {
  tag { "$name" }
  memory '4GB'

  input:
  tuple val(name), path(data), val(column_name)

  output:
  path('*.csv')

  script:
  """
  rnac genecards-suite lookup $column_name $data urs-info.pickle
  rnac genecards-suite $name $data urs-info.pickle .
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
    | fetch \
    | process \
    | set { data }
}
