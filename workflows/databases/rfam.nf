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
    < $query > families.tsv
  """
}

process fetch_info {
  when { params.databases.rfam.run }

  input:
  path(query)

  output:
  path("info.tsv")

  """
  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
    < $query > info.tsv
  """
}


process parse {
  maxForks 5

  input:
  tuple val(family), path(info), path(sequence_query)

  output:
  path('*.csv')

  """
  wget -O sequences.fa.gz 'http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/${family}.fa.gz'
  gzip -d sequences.fa.gz

  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
   -e "set @family='$family';\\. $sequence_query" > sequences.tsv

  rnac rfam parse $info sequences.tsv sequences.fa .
  """
}


workflow rfam {
  emit: data
  main:
    Channel.fromPath('files/import-data/rfam/select-families') | set { family_sql }
    Channel.fromPath('files/import-data/rfam/families.sql') | set { info_sql }
    Channel.fromPath('files/import-data/rfam/sequences.sql') | set { sequence_sql }

    info_sql | fetch_info | set { info }

    family_sql \
    | fetch_families \
    | splitCsv(sep='\t') \
    | map { row -> row[0] } \
    | combine(info) \
    | combine(sequence_sql) \
    | parse \
    | set { data }
}
