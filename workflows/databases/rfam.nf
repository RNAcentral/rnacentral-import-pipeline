process fetch_families {
  when { params.databases.rfam.run }

  input:
  path(query)

  output:
  path('families.tsv')

  """
  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
    $query > families.tsv
  """
}

process fetch_json {
  when { params.databases.rfam.run }

  output:
  path('*.json')

  """
  find $params.databases.rfam.remote -name '*.json' | xargs -I {} cp {} .
  """
}

process parse {
  input:
  tuple path(json), path(families)

  output:
  path('*.csv')

  """
  rnac rfam parse $json $families .
  """
}


workflow rfam {
  emit: data
  main:
    Channel.fromPath('files/import-data/rfam/families.sql') | set { families_sql }
    families_sql | fetch_families | set { families }

    fetch_json \
    | flatten \
    | combine(families) \
    | parse \
    | set { data }
}
