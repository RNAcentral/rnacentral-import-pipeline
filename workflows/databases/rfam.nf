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

process fetch_families_info {
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

process fetch_sequence_info {
  tag { "$family" }
  maxForks 10

  input:
  tuple val(family), path(sequence_family_query), path(sequence_seed_query)

  output:
  tuple val(family), path('sequences.tsv')

  """
  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
   -e "set @family='$family';\\. $sequence_family_query" > sequences_family.tsv

  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
   -e "set @family='$family';\\. $sequence_seed_query" > sequences_seed.tsv

  cat sequences_family.csv sequences_seed.csv > sequences.csv

  """
}

process parse {
  tag { "$family" }

  input:
  tuple val(family), path(sequence_info), path(families_info)

  output:
  path('*.csv')

  """
  wget -O sequences.fa.gz 'http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/${family}.fa.gz'
  gzip -d sequences.fa.gz

  rnac rfam parse $families_info $sequence_info sequences.fa .
  """
}


workflow rfam {
  emit: data
  main:
    Channel.fromPath('files/import-data/rfam/select-families.sql') | set { family_sql }
    Channel.fromPath('files/import-data/rfam/families.sql') | set { info_sql }
    Channel.fromPath('files/import-data/rfam/sequences_family.sql') | set { sequence_family_sql }
    Channel.fromPath('files/import-data/rfam/sequences_seed.sql') | set { sequence_seed_sql }

    info_sql | fetch_families_info | set { info }

    family_sql \
    | fetch_families \
    | splitCsv(sep: '\t', header: true) \
    | map { row -> row.rfam_acc } \
    | combine(sequence_sql) \
    | combine(sequence_seed_sql)
    | fetch_sequence_info \
    | combine(info) \
    | parse \
    | set { data }
}
