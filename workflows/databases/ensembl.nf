process fetch_metadata {
  input:
  path(query)

  output:
  path('families.tsv')

  when: params.databases.ensembl?.run

  script:
  """
  mysql \
    --host ${params.connections.rfam.host} \
    --port ${params.connections.rfam.port} \
    --user ${params.connections.rfam.user} \
    --database ${params.connections.rfam.database} < ${query} > families.tsv
  """
}

process find_urls {
  memory '4GB'

  output:
  path('species.csv')

  when: params.databases.ensembl?.run

  script:
  """
  rnac ensembl urls-for --location '${params.databases.ensembl.species_json_url}' species.csv
  """
}

process fetch_species_data {
  tag { species }
  errorStrategy { task.exitStatus == 8 && task.attempt <= 10 ? 'retry' : 'ignore' }
  maxRetries 10
  maxForks 10

  input:
  tuple val(species), val(taxid), val(embl_url), val(gff_url)

  output:
  tuple val(species), path("${species}.dat"), path("${species}.gff")

  script:
  """
  wget -O genes.embl.gz '${embl_url}'
  wget -O genes.gff3.gz '${gff_url}'

  gzip -dc genes.embl.gz > ${species}.dat

  zgrep '^#' genes.gff3.gz | grep -v '^###\$' > ${species}.gff
  zcat genes.gff3.gz | awk '{ if (\$3 !~ /CDS/) { print \$0 } }' >> ${species}.gff
  """
}

process parse_data {
  tag { "${embl.baseName}" }
  memory { 6.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(species), path(embl), path(gff), path(rfam)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac ensembl parse --family-file $rfam $embl $gff .
  """
}

workflow ensembl {
  main:
    channel.fromPath('files/import-data/rfam/families.sql') | fetch_metadata | set { families }

    find_urls() \
    | splitCsv \
    | filter { species, taxid, embl_url, gff_url ->
      !params.databases.ensembl.exclude.any { p -> species.toLowerCase() =~ p }
    } \
    | fetch_species_data \
    | combine(families) \
    | parse_data \
    | set { data }

  emit: data
}
