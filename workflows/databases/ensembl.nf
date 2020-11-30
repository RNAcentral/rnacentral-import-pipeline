process fetch_rfam {
  when { params.databases.ensembl.run }

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
  when { params.databases.ensembl.run }

  output:
  path('transcripts.gff')

  shell:
  '''
  rnac gencode urls |\
  xargs -I {} wget -O - {} |\
  gzip -d |\
  awk '{ if ($3 == "transcript") print $0 }' > transcripts.gff
  '''
}

process fetch_excludes {
  when { params.databases.ensembl.run }

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
  when { params.databases.ensembl.run }

  input:
  val(remote)

  output:
  path('species.txt')

  """
  curl --list-only $remote/ | xargs -I {} echo '$remote,{}' >> species.txt
  """
}

process fetch_species_data {
  tag { "$name" }

  input:
  tuple val(ftp), val(name)

  output:
  path("*.dat.gz")

  """
  wget "${ftp}/$name/*.dat.gz" .
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
    Channel.of(params.databases.ensembl.ftp, params.databases.ensembl.rapid_release.ftp) \
    | find_species \
    | splitCsv \
    | filter { _, filename ->
      params.databases.ensembl.data_file.exclude.inject(false) { agg, p -> agg || (filename =~ p) }
    } \
    | fetch_species_data \
    | flatten \
    | combine(rfam) \
    | combine(gencode) \
    | combine(excludes) \
    | parse_data \
    | set { data }
}
