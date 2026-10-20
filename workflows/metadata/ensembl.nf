process assemblies {
  input:
  path(connections)
  path(query)
  path(examples)
  path(known)

  output:
  path('*.{csv,parquet}')

  when:
  params.databases.ensembl?.vertebrates?.run

  script:
  """
  rnac ensembl assemblies $connections $query $examples $known assemblies.${params.writer_format}
  """
}

process fetch_compara {
  errorStrategy 'retry'
  maxRetries 10

  output:
  path('*.nt.fasta.gz')

  when: params.databases.ensembl?.vertebrates?.run

  script:
  """
  wget '$params.databases.ensembl.compara.remote'
  """
}

process process_compara {
  tag { "${gz.baseName}" }
  memory '17GB'

  input:
  path(gz)

  output:
  path("compara.${params.writer_format}")

  script:
  def out = "compara.${params.writer_format}"
  """
  zcat $gz | rnac ensembl compara - $out
  """
}

process proteins {
  input:
  path(connections)
  path(query)

  output:
  path("proteins.${params.writer_format}")

  when: params.databases.ensembl?.vertebrates?.run || params.databases.tarbase?.run || params.databases.lncbase?.run

  script:
  """
  rnac ensembl proteins $connections $query proteins.${params.writer_format}
  """
}

process coordinate_systems {
  errorStrategy 'retry'
  maxRetries 10

  input:
  path(connections)
  path(query)

  output:
  path("coordinate_systems.${params.writer_format}")

  when: params.databases.ensembl?.vertebrates?.run

  script:
  """
  rnac ensembl coordinate-systems $connections $query coordinate_systems.${params.writer_format}
  """
}

workflow compara {
  main:
    fetch_compara | flatten | process_compara | set { data }
  emit: data
}

workflow ensembl {
  main:
    channel.fromPath('config/databases.json') | set { conn }

    channel.fromPath('files/import-data/ensembl/proteins.sql') | set { protein_sql }
    channel.fromPath('files/import-data/ensembl/coordinate-systems.sql') | set { coordinate_systems_sql }

    channel.fromPath('files/import-data/ensembl/assemblies.sql') | set { assemblies_sql }
    channel.fromPath('files/import-data/ensembl/example-locations.json') | set { examples }
    channel.fromPath('files/import-data/ensembl/known-assemblies.sql') | set { known }

    channel.empty() \
    | mix(
      assemblies(conn, assemblies_sql, examples, known),
      coordinate_systems(conn, coordinate_systems_sql),
      proteins(conn, protein_sql),
      compara(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
