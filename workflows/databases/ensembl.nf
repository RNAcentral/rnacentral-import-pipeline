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

// Only the genesets whose checksums moved since the last release are fetched and parsed.
// See docs/incremental-parsing.md.
process ensembl_source_diff {
  input:
  path(species)

  output:
  path('to_fetch.csv'), emit: to_fetch
  path('scanned.txt'), emit: scanned
  path('sources.csv'), emit: sources

  when: params.databases.ensembl?.run

  script:
  def force = params.force_full_import ? '--force-full' : ''
  def only = params.databases.ensembl.species.collect { s -> "--species '${s}'" }.join(' ')
  """
  rnac ensembl source-diff $force $only $species to_fetch.csv scanned.txt sources.csv
  """
}

// Single reduce: retire what a re-parsed geneset no longer yields.
process ensembl_delta_diff {
  memory { 4.GB * task.attempt }
  errorStrategy { task.attempt < 3 ? 'retry' : 'finish' }

  input:
  path('accessions*.manifest')
  path(scanned)

  output:
  path('manifest.csv'), emit: manifest
  path('deletions.csv'), emit: deletions

  when: params.databases.ensembl?.run

  script:
  """
  # No manifest files at all when every geneset was unchanged; the diff still runs so
  # a geneset that has left the list is retired.
  cat accessions*.manifest > manifest.csv 2>/dev/null || : > manifest.csv
  rnac ensembl delta-diff manifest.csv $scanned deletions.csv
  """
}

process fetch_species_data {
  tag { species }
  // /nfs is only visible from a datamover node; the FTP it mirrors refuses most parallel downloads.
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  maxForks 30

  input:
  tuple val(species), val(taxid), val(embl_url), val(gff_url), val(assembly_id), val(signature)

  output:
  tuple val(species), val(assembly_id), val(embl_url), val(signature), path("${species}.dat"), path("${species}.gff"), optional: true

  script:
  def embl_path = embl_url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  def gff_path = gff_url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  """
  # species.json lists genesets that were never published, some with no EMBL file.
  if [ ! -e '${embl_path}' ] || [ ! -e '${gff_path}' ]; then
    echo "Skipping ${species}: no genes.embl.gz or genes.gff3.gz at ${embl_path}" >&2
    exit 0
  fi

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
  tuple val(species), val(assembly_id), val(source), val(signature), path(embl), path(gff), path(rfam)

  output:
  path('*.{csv,parquet}'), emit: data
  path('accessions.manifest'), emit: accessions

  script:
  """
  rnac ensembl parse --assembly-id '${assembly_id}' --source '${source}' --signature '${signature}' --family-file $rfam $embl $gff .
  """
}

workflow ensembl {
  main:
    channel.fromPath('files/import-data/rfam/families.sql') | fetch_metadata | set { families }

    // Filtered before the diff so an excluded species is never recorded as imported.
    find_urls() \
    | splitCsv \
    | filter { species, taxid, embl_url, gff_url, assembly_id ->
      !params.databases.ensembl.exclude.any { p -> species.toLowerCase() =~ p }
    } \
    | map { row -> row.join(',') } \
    | collectFile(name: 'species.csv', newLine: true) \
    | ensembl_source_diff

    ensembl_source_diff.out.to_fetch \
    | splitCsv \
    | fetch_species_data \
    | combine(families) \
    | parse_data

    parse_data.out.accessions | collect | ifEmpty([]) | set { accessions }
    ensembl_delta_diff(accessions, ensembl_source_diff.out.scanned)

    // manifest.csv, deletions.csv and sources.csv ride the generic wiring in import-data.nf.
    parse_data.out.data \
    | mix(
      ensembl_delta_diff.out.manifest,
      ensembl_delta_diff.out.deletions,
      ensembl_source_diff.out.sources,
    ) \
    | set { data }

  emit: data
}
