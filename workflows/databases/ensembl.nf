process fetch_rfam {
  when { params.databases.ensembl.run }

  input:
  path(query)

  output:
  path('families.tsv')

  """
  mysql \
    --host $params.connections.rfam.host \
    --port $params.connections.rfam.port \
    --user $params.connections.rfam.user \
    --database $params.connections.rfam.database \
    $query > families.tsv
  """
}

process fetch_gencode {
  when { params.databases.ensembl.run }

  output:
  path('transcripts.gff')

  shell:
  '''
  rnac gencode urls |\
  xargs -I {} wget -O - {} |\
  gzip -d |\
  awk '{ if ($3 == "transcript") print $0 }' > transcripts.gff
  '''
}

process fetch_excludes {
  when { params.databases.ensembl.run }

  input:
  path(urls)

  output:
  path('ids')

  """
  cat $urls |\
  xargs -I {} wget -O - {} |\
  zgrep '^>' |\
  grep 'processed_transcript' |\
  cut -d ' ' -f1 |\
  tr -d '>' > ids
  """
}

process find_urls {
  executor 'local'
  when { params.databases.ensembl.run }

  input:
  val(remote)

  output:
  path('species.txt')

  """
  rnac ensembl vertebrates urls-for $remote > species.txt
  """
}

process fetch_species_data {
  tag { "$name" }

  input:
  tuple val(name), val(dat_path), val(gff_path)

  output:
  tuple val(name), path("*.dat"), path('*.gff')

  """
  wget '$dat_path'
  wget '$gff_path'
  gzip -d *.gz
  """
}

process parse_data {
  tag { "$name" }

  input:
  tuple val(name), path(embl), path(gff), path(rfam), path(gencode), path(exclude)

  output:
  path('*.csv')

  """
  rnac ensembl vertebrates parse $embl $gff $rfam $gencode $exclude .
  """
}

workflow ensembl {
  emit: data
  main:
    fetch_gencode | set { gencode }
    Channel.fromPath('files/import-data/rfam/families.sql') | fetch_rfam | set { rfam }
    Channel.fromPath('files/import-data/ensembl/exclude-urls.txt') | fetch_excludes | set { excludes}

    Channel.of(params.databases.ensembl.ftp) \
    | find_urls \
    | splitCsv \
    | filter { name, dat_url, gff_url ->
      params.databases.ensembl.data_file.exclude.any { p -> name =~ p }
    } \
    | fetch_species_data \
    | map { name, dat_files, gff_file -> 
      dat_files.collect { [name, it, gff_file] }
    } \
    | combine(rfam) \
    | combine(gencode) \
    | combine(excludes) \
    | parse_data \
    | set { data }
}
