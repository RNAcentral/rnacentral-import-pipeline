process find_urls {
  when { params.databases.pirbase.run }
  executor 'local'

  output:
  path("urls.txt")

  """
  rnac pirbase urls-for $params.databases.pirbase.remote urls.txt
  """
}

process find_known {
  when { params.databases.pirbase.run }

  input:
  path(query)

  output:
  path(known)

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > known
  """
}

process parse_data {
  tag { "${url.split('/').last()}" }
  memory '5GB'
  errorStrategy 'ignore'

  input:
  tuple val(url), path(known)

  output:
  path("*.csv")

  """
  wget -O data.json.gz $url
  gzip -d data.json.gz
  rnac pirbase parse data.json $known .
  """
}

workflow pirbase {
  emit: data_files
  main:
    Channel.fromPath('files/import-data/pirbase/known-md5.sql') | set { query }

    find_urls \
    | splitCsv \
    | map { row -> row[0] } \
    | filter { url -> !url.contains('piRBase_tbe.json') } \
    | combine(find_known(query)) \
    | parse_data \
    | flatten \
    | set { data_files }
}
