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
  // species.json is on the EBI FTP, which refuses connections now and then
  errorStrategy { sleep(task.attempt * 60000 as long); return 'retry' }
  maxRetries 3

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
  // /nfs is only visible from a datamover node; the FTP it mirrors refuses most parallel downloads.
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  maxForks 30

  input:
  tuple val(species), val(taxid), val(embl_url), val(gff_url), val(assembly_id)

  output:
  tuple val(species), val(assembly_id), path("${species}.dat"), path("${species}.gff")

  script:
  def embl_path = embl_url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  def gff_path = gff_url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  """
  cp '${embl_path}' genes.embl.gz
  cp '${gff_path}' genes.gff3.gz

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
  tuple val(species), val(assembly_id), path(embl), path(gff), path(rfam)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac ensembl parse --assembly-id '${assembly_id}' --family-file $rfam $embl $gff .
  """
}

workflow ensembl {
  main:
    channel.fromPath('files/import-data/rfam/families.sql') | fetch_metadata | set { families }

    find_urls() \
    | splitCsv \
    | filter { species, taxid, embl_url, gff_url, assembly_id ->
      !params.databases.ensembl.exclude.any { p -> species.toLowerCase() =~ p }
    } \
    | filter { species, taxid, embl_url, gff_url, assembly_id ->
      !params.databases.ensembl.species || species.toLowerCase() in params.databases.ensembl.species
    } \
    | fetch_species_data \
    | combine(families) \
    | parse_data \
    | set { data }

  emit: data
}
