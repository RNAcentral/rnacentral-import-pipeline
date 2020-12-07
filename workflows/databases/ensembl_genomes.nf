process find_species {
  tag { "$name" }

  input:
  tuple val(name), val(remote)

  output:
  path('species.txt')

  """
  rnac ensembl genomes urls-for $name $remote > species.txt
  """
}

process fetch_species_data {
  tag { "$species" }

  input:
  tuple val(name), val(species), val(dat_path), val(gff_path)

  output:
  tuple val(name), path("*.dat"), path('*.gff')

  """
  wget '$dat_path'
  wget '$gff_path'
  gzip -d *.gz
  """
}

process parse_data {
  tag { "${embl.basename}" }

  input:
  tuple val(name), path(embl), path(gff)

  output:
  path('*.csv')

  """
  rnac ensembl genomes parse $name $embl $gff
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
    | filter { name -> params.databases.ensembl[name].run } \
    | map { name -> [name, params.databases.ensembl[name].ftp_host] } \
    | find_species \
    | splitCsv \
    | filter { name, species, data_path, gff_path ->
      params.ensembl[name].exclude.any { p -> species =~ p }
    } \
    | fetch_species_data \
    | flatMap { name, dat_files, gff_file ->
      dat_files.collect { [name, it, gff_file] }
    } \
    | parse_data \
    | set { data }
}
