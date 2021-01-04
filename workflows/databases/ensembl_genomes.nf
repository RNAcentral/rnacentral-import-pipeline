process fetch_metadata {
  when { params.databases.ensembl.vertebrates.run }

  input:
  path(query)

  output:
  path('families.tsv')

  """
  mysql \
    --host ${params.connections.rfam.host} \
    --port ${params.connections.rfam.port} \
    --user ${params.connections.rfam.user} \
    --database ${params.connections.rfam.database} < ${query} > families.tsv
  """
}

process find_species {
  tag { "$division" }
  memory '4GB'

  input:
  tuple val(division), val(remote)

  output:
  path('species.txt')

  """
  rnac ensembl urls-for $division $remote > species.txt
  """
}

process fetch_species_data {
  tag { "$species" }
  memory '4GB'
  errorStrategy 'retry'
  maxRetries 10
  maxForks 5

  input:
  tuple val(division), val(species), val(dat_path), val(gff_path)

  output:
  tuple val(division), path("*.dat"), path("${species}.gff")

  """
  wget '$dat_path'
  wget '$gff_path'
  zgrep '^#' *.gff3.gz | grep -v '^###\$' > ${species}.gff
  zcat *.gff3.gz | awk '{ if (\$3 !~ /CDS/) { print \$0 } }' >> ${species}.gff
  gzip -d *.gz
  """
}

process parse_data {
  tag { "${embl.baseName}" }
  memory '6GB'

  input:
  tuple val(division), path(embl), path(gff), path(rfam)

  output:
  path('*.csv')

  """
  rnac ensembl parse --family-file $rfam $division $embl $gff
  """
}

workflow ensembl_genomes {
  emit: data
  main:
    Channel.fromPath('files/import-data/rfam/families.sql') | set { rfam }

    Channel.fromList([
      'plants',
      'fungi',
      'protists',
      'metazoa',
      'vertebrates',
    ]) \
    | filter { division -> params.databases.ensembl[division].run } \
    | map { division -> [division, params.databases.ensembl[division].ftp_host] } \
    | find_species \
    | splitCsv \
    | filter { division, species, data_files, gff_path ->
      !params.ensembl[division].exclude.any { p -> species =~ p }
    } \
    | fetch_species_data \
    | flatMap { division, data_files, gff_file ->
      (data_files instanceof ArrayList) ? data_files.collect { [division, it, gff_file] } : [[division, data_files, gff_file]]
    } \
    | parse_data \
    | set { data }
}
