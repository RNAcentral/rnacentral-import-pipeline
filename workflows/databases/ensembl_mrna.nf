include { fetch_species_data } from './ensembl'

process find_mrna_urls {
  memory '4GB'
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"

  output:
  path('species.csv')

  when: params.databases.ensembl_coding?.run

  script:
  """
  rnac ensembl urls-for --kind mrna --location '${params.databases.ensembl.nfs_path}/species.json' species.csv
  """
}

// Only genesets whose checksums moved since the last coding release are fetched and
// parsed; ENSEMBL_CODING keeps its own signatures. See docs/incremental-parsing.md.
process coding_source_diff {
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"

  input:
  path(species)

  output:
  path('to_fetch.csv'), emit: to_fetch
  path('scanned.txt'), emit: scanned
  path('sources.csv'), emit: sources

  when: params.databases.ensembl_coding?.run

  script:
  def force = params.force_full_import ? '--force-full' : ''
  """
  rnac ensembl source-diff --database ENSEMBL_CODING --ftp-url '${params.databases.ensembl.ftp_url}' --nfs-path '${params.databases.ensembl.nfs_path}' $force $species to_fetch.csv scanned.txt sources.csv
  """
}

process coding_delta_diff {
  memory { 4.GB * task.attempt }
  errorStrategy { task.attempt < 3 ? 'retry' : 'finish' }

  input:
  path('accessions*.manifest')
  path(scanned)

  output:
  path('manifest.csv'), emit: manifest
  path('deletions.csv'), emit: deletions

  when: params.databases.ensembl_coding?.run

  script:
  """
  cat accessions*.manifest > manifest.csv 2>/dev/null || : > manifest.csv
  rnac ensembl delta-diff --database ENSEMBL_CODING manifest.csv $scanned deletions.csv
  """
}

process parse_mrna {
  tag { species }
  memory { 6.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(species), val(_assembly_id), val(source), val(signature), path(embl), path(gff)

  output:
  path('*.{csv,parquet}'), emit: data
  path('accessions.manifest'), emit: accessions

  script:
  """
  rnac ensembl parse-mrna --source '${source}' --signature '${signature}' $embl $gff .
  """
}

// Homology pairs cross species and its load replaces every coding row, so it needs
// every selected species' GFF each run, not only the genesets that changed.
process fetch_gff {
  tag { species }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  maxForks 30

  input:
  tuple val(species), val(gff_url)

  output:
  path("${species}.gff"), optional: true

  script:
  def gff_path = gff_url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  """
  if [ ! -e '${gff_path}' ]; then
    echo "Skipping ${species}: no genes.gff3.gz at ${gff_path}" >&2
    exit 0
  fi
  zcat '${gff_path}' | awk -F'\\t' '\$3 == "gene" || \$3 == "mRNA"' > ${species}.gff
  """
}

process fetch_homology {
  tag { species }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  maxForks 10

  input:
  tuple val(species), val(url)

  output:
  path("${species}.homology.tsv.gz")

  script:
  def homology_path = url.replace(params.databases.ensembl.ftp_url, params.databases.ensembl.nfs_path)
  """
  cp '${homology_path}' ${species}.homology.tsv.gz
  """
}

process mrna_homology {
  memory '4GB'

  input:
  path(homologies)
  path(gffs)

  output:
  path("compara.${params.writer_format}")

  script:
  def homology_args = homologies.collect { "--homology $it" }.join(' ')
  def gff_args = gffs.collect { "--gff $it" }.join(' ')
  """
  rnac ensembl mrna-homology $homology_args $gff_args compara.${params.writer_format}
  """
}

workflow ensembl_mrna {
  main:
    find_mrna_urls() \
    | splitCsv \
    | filter { row -> row[0].toLowerCase() in params.databases.ensembl_coding.species } \
    | set { selected }

    // Filtered before the diff, so a species taken off the list is gone and retired.
    selected \
    | map { row -> row.join(',') } \
    | collectFile(name: 'species.csv', newLine: true) \
    | coding_source_diff

    coding_source_diff.out.to_fetch \
    | splitCsv \
    | map { species, taxid, embl_url, gff_url, assembly_id, _homology_url, signature -> [species, taxid, embl_url, gff_url, assembly_id, signature] } \
    | fetch_species_data \
    | parse_mrna

    parse_mrna.out.accessions | collect | ifEmpty([]) | set { accessions }
    coding_delta_diff(accessions, coding_source_diff.out.scanned)

    selected \
    | filter { row -> row[5] } \
    | map { species, _taxid, _embl_url, _gff_url, _assembly_id, homology_url -> [species, homology_url] } \
    | fetch_homology \
    | set { homologies }

    selected \
    | map { species, _taxid, _embl_url, gff_url, _assembly_id, _homology_url -> [species, gff_url] } \
    | fetch_gff \
    | set { gffs }

    mrna_homology(homologies.collect(), gffs.collect())

    // manifest.csv, deletions.csv and sources.csv ride the generic wiring in import-data.nf.
    parse_mrna.out.data \
    | mix(
      mrna_homology.out,
      coding_delta_diff.out.manifest,
      coding_delta_diff.out.deletions,
      coding_source_diff.out.sources,
    ) \
    | set { data }

  emit: data
}
