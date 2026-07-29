process generic {
  tag { "$name" }

  input:
  tuple val(name), path(query)

  output:
  path('*.{csv,parquet}')

  when: params.databases.ensembl?.vertebrates?.run || params.databases.rfam?.run

  script:
  def conn = params.connections.rfam
  """
  mysql \
    --host $conn.host \
    --database $conn.database \
    --user $conn.user \
    --port $conn.port \
    < $query > raw.tsv
  rnac rfam $name raw.tsv
  """
}

workflow rfam {
  main:
    channel.fromPath('files/import-data/rfam/{clans,families}.sql') \
    | map { fn -> [fn.baseName, fn] } \
    | generic \
    | flatten \
    | set { data }
  emit: data
}
