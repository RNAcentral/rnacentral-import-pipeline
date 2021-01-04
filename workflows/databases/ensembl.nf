process find_urls {
  memory '20GB'

  output:
  tuple val(division), path('species.txt')

  """
  rnac ensembl urls-for $division $params.databases.ensembl.vertebrates.ftp_host > species.txt
  """
}

process fetch_species_data {
  tag { "$species" }
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
  memory { '6GB' }

  input:
  tuple val(divsion), path(embl), path(gff), path(rfam)

  output:
  path('*.csv')

  """
  rnac ensembl parse $division --family-file $rfam $embl $gff .
  """
}

workflow ensembl {
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
    | find_urls \
    | splitCsv \
    | filter { division, species, dat_url, gff_url ->
      !params.databases.ensembl[division].exclude.any { p -> species.toLowerCase() =~ p }
    } \
    | fetch_species_data \
    | flatMap { division, dat_files, gff_file ->
      (dat_files instanceof ArrayList) ? dat_files.collect { [division, it, gff_file] } : [[division, dat_files, gff_file]]
    } \
    | combine(fetch_metadata(rfam)) \
    | parse_data \
    | set { data }
}
