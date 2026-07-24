process assemblies {
  input:
  path(connections)
  path(query)
  path(examples)
  path(known)

  output:
  path('*.csv')

  when:
  params.databases.ensembl?.vertebrates?.run

  script:
  """
  rnac ensembl assemblies $connections $query $examples $known
  """
}

process fetch_compara {
  errorStrategy 'retry'
  maxRetries 10

  output:
  path('*.nt.fasta.gz')

  when: { params.databases.ensembl?.vertebrates?.run }

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
  path('compara.csv')

  script:
  """
  zcat $gz | rnac ensembl compara - compara.csv
  """
}

process proteins {
  input:
  path(connections)
  path(query)

  output:
  path('proteins.csv')

  when: { params.databases.ensembl?.vertebrates?.run || params.databases.tarbase?.run || params.databases.lncbase?.run }

  script:
  """
  rnac ensembl proteins $connections $query proteins.csv
  """
}

process coordinate_systems {
  errorStrategy 'retry'
  maxRetries 10

  input:
  path(connections)
  path(query)

  output:
  path('coordinate_systems.csv')

  when: { params.databases.ensembl?.vertebrates?.run }

  script:
  """
  rnac ensembl coordinate-systems $connections $query coordinate_systems.csv
  """
}

process karyotypes {
  errorStrategy 'retry'
  maxRetries 10

  output:
  path('karyotypes.csv')

  when: { params.databases.ensembl?.vertebrates?.run }

  script:
  """
  rnac ensembl karyotypes karyotypes.csv
  """
}

workflow compara {
  main:
    fetch_compara | flatten | process_compara | set { data }
  emit: data
}

workflow ensembl {
  main:
    Channel.fromPath('config/databases.json') | set { conn }

    Channel.fromPath('files/import-data/ensembl/proteins.sql') | set { protein_sql }
    Channel.fromPath('files/import-data/ensembl/coordinate-systems.sql') | set { coordinate_systems_sql }

    Channel.fromPath('files/import-data/ensembl/assemblies.sql') | set { assemblies_sql }
    Channel.fromPath('files/import-data/ensembl/example-locations.json') | set { examples }
    Channel.fromPath('files/import-data/ensembl/known-assemblies.sql') | set { known }

    Channel.empty() \
    | mix(
      assemblies(conn, assemblies_sql, examples, known),
      coordinate_systems(conn, coordinate_systems_sql),
      proteins(conn, protein_sql),
      karyotypes(),
      compara(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
