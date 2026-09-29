include { fetch_species_data } from './ensembl'

process find_mrna_urls {
  memory '4GB'

  output:
  path('species.csv')

  when: params.databases.ensembl_mrna?.run

  script:
  """
  rnac ensembl urls-for --kind mrna --location '${params.databases.ensembl.species_json_url}' species.csv
  """
}

process parse_mrna {
  tag { species }
  memory { 6.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(species), path(embl), path(gff)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac ensembl parse-mrna $embl $gff .
  """
}

process fetch_homology {
  tag { species }
  maxForks 10

  input:
  tuple val(species), val(url)

  output:
  path("${species}.homology.tsv.gz")

  script:
  """
  wget -O ${species}.homology.tsv.gz '${url}'
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
    | filter { row -> row[0].toLowerCase() in params.databases.ensembl_mrna.species } \
    | set { selected }

    selected \
    | map { species, taxid, embl_url, gff_url, _homology_url -> [species, taxid, embl_url, gff_url] } \
    | fetch_species_data \
    | set { fetched }

    selected \
    | filter { row -> row[4] } \
    | map { species, _taxid, _embl_url, _gff_url, homology_url -> [species, homology_url] } \
    | fetch_homology \
    | set { homologies }

    fetched | parse_mrna | set { parsed }

    mrna_homology(
      homologies.collect(),
      fetched.map { _species, _embl, gff -> gff }.collect(),
    )

    parsed.mix(mrna_homology.out) | set { data }

  emit: data
}
