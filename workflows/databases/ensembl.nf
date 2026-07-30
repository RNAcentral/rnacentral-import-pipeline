process fetch_metadata {
  input:
  path(query)

  output:
  path('families.tsv')

  when: params.databases.ensembl?._any?.run

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
  memory '20GB'

  input:
  val(division)

  output:
  path('species.txt')

  when: params.databases.ensembl[division]?.run

  script:
    """
    rnac ensembl urls-for $division ${params.databases.ensembl[division].ftp_host} species.txt
    """
}

process fetch_species_data {
  tag { "$species" }
  errorStrategy { task.exitStatus == 8 && task.attempt <= 10 ? 'retry' : 'ignore' }
  maxRetries 10
  maxForks 10

  input:
  tuple val(division), val(species), val(dat_path), val(gff_path)

  output:
  tuple val(division), path("*.dat"), path("${species}.gff")

  script:
  """
  lowercase_assembly_url() {
    PYTHONPATH="${workflow.launchDir}" python -c "from rnacentral_pipeline.databases.ensembl.url_helpers import lowercase_assembly_in_url; print(lowercase_assembly_in_url('\$1'))"
  }

  resolve_ftp_urls() {
    PYTHONPATH="${workflow.launchDir}" python -c "import sys; from rnacentral_pipeline.databases.ensembl.url_helpers import resolve_ftp_urls; sys.stdout.write('\\n'.join(resolve_ftp_urls(sys.argv[1])))" "\$1"
  }

  resolve_ftp_urls '$dat_path' > dat_urls.txt
  test -s dat_urls.txt
  xargs -n 1 wget < dat_urls.txt
  wget '$gff_path' || wget "\$(lowercase_assembly_url '$gff_path')"
  zgrep '^#' *.gff3.gz | grep -v '^###\$' > ${species}.gff
  zcat *.gff3.gz | awk '{ if (\$3 !~ /CDS/) { print \$0 } }' >> ${species}.gff
  gzip -d *.gz
  """
}

process parse_data {
  tag { "${embl.baseName}" }
  memory { 6.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(division), path(embl), path(gff), path(rfam)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac ensembl parse $division --family-file $rfam $embl $gff .
#  rnac ensembl pseudogenes $division $embl ensembl-pseudogenes.csv
  """
}

workflow ensembl {
  main:
    channel.fromPath('files/import-data/rfam/families.sql') | set { rfam }

    channel.fromList([
      'plants',
      'fungi',
      'protists',
      'metazoa',
      'vertebrates',
    ]) \
    | filter { division -> params.databases.ensembl[division]?.run } \
    | find_urls \
    | splitCsv \
    | filter { division, species, _dat_url, _gff_url ->
      !params.databases.ensembl[division].exclude.any { p -> species.toLowerCase() =~ p }
    } \
    | fetch_species_data \
    | flatMap { division, dat_files, gff_file ->
      (dat_files instanceof ArrayList) ? dat_files.collect { f -> [division, f, gff_file] } : [[division, dat_files, gff_file]]
    } \
    | combine(fetch_metadata(rfam)) \
    | parse_data \
    | set { data }

  emit: data
}
