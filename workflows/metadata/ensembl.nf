process assemblies {
  when { params.databases.ensembl.run }

  input:
  path(connections)
  path(query)
  path(examples)
  path(known)

  output:
  path('*.csv')

  """
  rnac ensembl assemblies $connections $query $examples $known
  """
}

process fetch_compara {
  when { params.databases.ensembl.run }

  output:
  path('*.nt.fasta.gz')

  """
  wget '$params.metadata.ensembl_compara.inputs.fasta'
  """
}

process process_compara {
  memory '20GB'

  input:
  path(gz)

  output:
  path('compara.csv')

  """
  zcat $gz | rnac ensembl compara - compara.csv
  """
}

process proteins {
  when { params.databases.ensembl.run || params.databases.tarbase.run || params.databases.lncbase.run }

  input:
  path(connections)
  path(query)

  output:
  path('proteins.csv')

  """
  rnac ensembl proteins $connections $query proteins.csv
  """
}

process coordinate_systems {
  when { params.databases.ensembl.run }

  input:
  path(connections)
  path(query)

  output:
  path('coordinate_systems.csv')

  """
  rnac ensembl coordinate_systems $connections $query coordinate_systems.csv
  """
}

process karyotypes {
  memory '8GB'
  when { params.databases.ensembl.run }

  output:
  path('karyotypes.csv')

  """
  rnac ensembl karyotypes karyotypes.csv
  """
}

workflow compara {
  emit: data
  main:
    fetch_compara | flatten | process_compara | set { data }
}

workflow ensembl {
  emit: data
  main:
    Channel.fromPath('config/connections.json') | set { conn }

    Channel.fromPath('files/import-data/ensembl/proteins.sql') | set { protein_sql }
    Channel.fromPath('files/import-data/ensembl/coordinate_systems.sql') | set { coordinate_systems_sql }

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
}
