process find_species {
  when { params.databases.ensembl_protists.run }
  executor 'local'

  output:
  path('species.txt')

  """
  curl --list-only $params.ensembl_protists.ftp > species.txt
  """
}

process fetch_species_data {
  when { params.databases.ensembl_protists.run }
  tag { "$name" }

  input:
  val(name)

  output:
  path("*.dat.gz")

  """
  wget "${params.ensembl_protists.ftp}/$name/*.dat.gz" .
  """
}

process parse_data {
  tag { "${embl.basename}" }

  input:
  path(embl)

  output:
  path('*.csv')

  """
  zcat $embl > data.dat
  rnac external ensembl_protists data.dat .
  """
}

workflow ensembl_protists {
  emit: data
  main:
    find_species \
    | splitCsv \
    | map { row -> row[0] } \
    | filter { filename ->
      params.ensembl_protists.exclude.inject(false) { agg, p -> agg || (filename =~ p) }
    } \
    | fetch_species_data \
    | flatten \
    | parse_data \
    | set { data }
}
