process fetch_annotation {
  tag { species.annotation }
  memory '2GB'
  errorStrategy 'terminate'

  input:
  val species

  output:
  tuple val(species), path("${species.annotation}.txt")

  when: { params.databases.circpedia?.run }

  script:
  """
  wget \
    "${params.databases.circpedia.annotation_base_url}${species.annotation}.txt.zip" \
    -O "${species.annotation}.txt.zip"
  unzip "${species.annotation}.txt.zip"
  """
}

process fetch_fasta {
  tag { species.fasta }
  memory '2GB'
  errorStrategy 'terminate'

  input:
  tuple val(species), path(annotation_file)

  output:
  tuple val(species), path(annotation_file), path("${species.fasta}.fa")

  when: { params.databases.circpedia?.run }

  script:
  """
  wget --no-check-certificate \
    "${params.databases.circpedia.fasta_base_url}${species.fasta}.fa.zip" \
    -O "${species.fasta}.fa.zip"
  unzip "${species.fasta}.fa.zip"
  """
}

process parse_data {
  tag { annotation_file.name }
  memory '8GB'

  input:
  tuple val(species), path(annotation_file), path(fasta_file)

  output:
  path('*.csv')

  script:
  """
  rnac circpedia parse \
    --assembly ${species.assembly} \
    $annotation_file \
    $fasta_file \
    .
  """
}

workflow circpedia {
  main:
    channel.fromList(params.databases.circpedia.species)
    | fetch_annotation
    | fetch_fasta
    | parse_data
    | set { data }

  emit: data
}
