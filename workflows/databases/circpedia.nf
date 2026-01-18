process fetch_data {
  when { params.databases.circpedia.run }
  memory '4GB'

  output:
  tuple path('circpedia_annotation.txt'), path('circpedia_sequences.fa')

  """
  # Fetch CIRCpedia V3 annotation file
  wget --no-check-certificate $params.databases.circpedia.remote.annotation -O circpedia_annotation.txt || \
  scp $params.databases.circpedia.remote.annotation circpedia_annotation.txt

  # Fetch CIRCpedia V3 FASTA sequence file
  wget --no-check-certificate $params.databases.circpedia.remote.fasta -O circpedia_sequences.fa || \
  scp $params.databases.circpedia.remote.fasta circpedia_sequences.fa

  # If the downloaded files are compressed, uncompress them
  if [[ -f circpedia_annotation.txt.gz ]]; then
    gunzip circpedia_annotation.txt.gz
  fi
  if [[ -f circpedia_sequences.fa.gz ]]; then
    gunzip circpedia_sequences.fa.gz
  fi
  """
}

process parse_data {
  tag { "$annotation_file.name" }
  memory '8GB'

  input:
  tuple path(annotation_file), path(fasta_file), path(taxonomy)

  output:
  path('*.csv')

  """
  # Parse CIRCpedia data with both annotation and FASTA files
  # Use assembly ID if specified in config, otherwise omit
  if [ -n "$params.databases.circpedia.assembly" ]; then
    rnac circpedia parse $taxonomy $annotation_file $fasta_file . --assembly $params.databases.circpedia.assembly
  else
    rnac circpedia parse $taxonomy $annotation_file $fasta_file .
  fi
  """
}

workflow circpedia {
  take: taxonomy
  emit: data
  main:
    fetch_data \
    | combine(taxonomy) \
    | parse_data \
    | set { data }
}
