process fetch_data {
  when { params.databases.circpedia.run }
  memory '4GB'

  output:
  path('circpedia_*.csv')

  """
  # Fetch CIRCpedia V3 data
  # The remote location should be configured in databases.config
  wget --no-check-certificate $params.databases.circpedia.remote -O circpedia_data.csv || \
  scp $params.databases.circpedia.remote circpedia_data.csv

  # If the downloaded file is compressed, uncompress it
  if [[ -f circpedia_data.csv.gz ]]; then
    gunzip circpedia_data.csv.gz
  elif [[ -f circpedia_data.tar.gz ]]; then
    tar xzf circpedia_data.tar.gz
  fi

  # Rename to standard format
  # If multiple files, they should already match the pattern circpedia_*.csv
  if [[ ! -f circpedia_*.csv ]]; then
    mv circpedia_data.csv circpedia_all.csv 2>/dev/null || true
  fi
  """
}

process parse_data {
  tag { "$csv_file.name" }
  memory '8GB'

  input:
  tuple path(csv_file), path(taxonomy)

  output:
  path('*.csv')

  """
  # Parse CIRCpedia data
  # Use assembly ID if specified in config, otherwise omit
  if [ -n "$params.databases.circpedia.assembly" ]; then
    rnac circpedia parse $taxonomy $csv_file . --assembly $params.databases.circpedia.assembly
  else
    rnac circpedia parse $taxonomy $csv_file .
  fi
  """
}

workflow circpedia {
  take: taxonomy
  emit: data
  main:
    fetch_data \
    | flatten \
    | combine(taxonomy) \
    | parse_data \
    | set { data }
}
