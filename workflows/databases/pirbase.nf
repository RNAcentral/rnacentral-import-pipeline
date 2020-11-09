process find_urls {
  executor 'local'

  when { params.databases.pirbase.run }

  output:
  path("urls.txt")

  """
  rnac pirbase urls-for $params.database.pirbase.url urls.txt
  """
}

process fetch_and_parse {
  input:
  val(url)

  output:
  path("*.csv")

  """
  wget -O data.json.gz $url
  gzip -d data.json.gz
  rnac pirbase parse data.json .
  """
}


workflow pirbase {
  emit: data_files
  main:
    | find_urls \
    | spltCsv \
    | map { row -> row[0] } \
    | fetch_and_parse \
    | flatten \
    | set { data_files }
}
