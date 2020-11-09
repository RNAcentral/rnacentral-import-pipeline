process fetch_rfam {
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

process fetch_gencode {
  output:
  path('transcripts.gff')

  """
  rnac gencode urls |\
  xargs -I {} wget -O - {} |\
  gzip -d |\
  awk '{ if ($3 == "transcript") print $0 }' > transcripts.gff
  """
}

process fetch_excludes {
  input:
  path(urls)

  output:
  path('ids')

  """
  cat $urls |\
  xargs -I {} wget -O - {} |\
  zgrep '^>' |\
  grep 'processed_transcript' |\
  cut -d ' ' -f1 |\
  tr -d '>' > ids
  """
}

process find_species {
  executor 'local'

  output:
  path('species.txt')

  """
  curl --list-only $params.ensembl.ftp > species.txt
  """
}

process fetch_species_data {
  tag { "$name" }

  input:
  val(name)

  output:
  path("*.dat.gz')

  """
  wget "${params.ensembl.ftp}/$name/*.dat.gz" .
  """
}

process parse_data {
  tag { "${embl.basename}" }

  input:
  tuple path(embl), path(rfam), path(gencode), path(exclude)

  output:
  path('*.csv')

  """
  zcat $embl > data.dat
  rnac external ensembl data.dat $rfam $gencode $exclude .
  """
}

workflow ensembl {
  emit: data
  main:
    fetch_gencode | set { gencode }
    Channel.fromPath('files/import-data/rfam/families.sql') | fetch_rfam | set { rfam }
    Channel.fromPath('files/import-data/ensembl/exclude-urls.txt') | fetch_excludes | set { excludes}

    fetch_species \
    | splitCsv \
    | map { row -> row[0] } \
    | filter { filename ->
      params.ensembl.data_file.exclude.inject(false) { agg, p -> agg || (filename =~ p) }
    } \
    | fetch_species_data \
    | flatten \
    | combine(rfam, gencode, excludes) \
    | parse_data \
    | set { data }
}
