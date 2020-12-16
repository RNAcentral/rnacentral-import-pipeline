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

process find_urls {
  when { params.databases.ensembl.vertebrates.run }
  memory '20GB'

  output:
  path('species.txt')

  """
  rnac ensembl urls-for vertebrates $params.databases.ensembl.vertebrates.ftp_host > species.txt
  """
}

process fetch_species_data {
  tag { "$name" }
  errorStrategy 'retry'
  maxRetries 10
  maxForks 5

  input:
  tuple val(name), val(dat_path), val(gff_path)

  output:
  tuple val(name), path("*.dat"), path("${name}.gff")

  """
  wget '$dat_path'
  wget '$gff_path'
  zgrep '^#' *.gff3.gz | grep -v '^###\$' > ${name}.gff
  zcat *.gff3.gz | awk '{ if (\$3 !~ /CDS/) { print \$0 } }' >> ${name}.gff
  gzip -d *.gz
  """
}

process parse_data {
  tag { "${embl.baseName}" }
  memory { 5.GB }

  input:
  tuple val(name), path(embl), path(gff), path(rfam)

  output:
  path('*.csv')

  """
  rnac ensembl parse vertebrates --family-file $rfam $embl $gff .
  """
}

workflow ensembl {
  emit: data
  main:
    Channel.fromPath('files/import-data/rfam/families.sql') | set { rfam }

    find_urls \
    | splitCsv \
    | filter { name, dat_url, gff_url ->
      !params.databases.ensembl.vertebrates.exclude.any { p -> name.toLowerCase() =~ p }
    } \
    | fetch_species_data \
    | flatMap { name, dat_files, gff_file ->
      (dat_files instanceof ArrayList) ? dat_files.collect { [name, it, gff_file] } : [[name, dat_files, gff_file]]
    } \
    | combine(fetch_metadata(rfam)) \
    | parse_data \
    | set { data }
}
