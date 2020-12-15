process find_species {
  tag { "$division" }

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
  tag { "${embl.basename}" }

  input:
  tuple val(division), path(embl), path(gff)

  output:
  path('*.csv')

  """
  rnac ensembl genomes parse $division $embl $gff
  """
}

workflow ensembl_genomes {
  emit: data
  main:
    Channel.fromList([
      'plants',
      'fungi',
      'protists',
      'metazoa',
    ]) \
    | filter { division -> params.databases.ensembl[division].run } \
    | map { division -> [name, params.databases.ensembl[division].ftp_host] } \
    | find_species \
    | splitCsv \
    | filter { division, species, data_files, gff_path ->
      params.ensembl[division].exclude.any { p -> species =~ p }
    } \
    | fetch_species_data \
    | flatMap { division, data_files, gff_file ->
      (data_files instanceof ArrayList) ? data_files.collect { [division, it, gff_file] } : [[division, data_files, gff_file]]
    } \
    | parse_data \
    | set { data }
}
