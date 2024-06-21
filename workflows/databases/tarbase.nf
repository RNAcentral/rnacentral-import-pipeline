process fetch_mirbase {
  when { params.databases.tarbase.run }

  input:
  path(query)

  output:
  path("result.fasta")

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > result.json
  json2fasta.py result.json result.fasta
  """
}

process process {
  when { params.databases.tarbase.run }
  tag { "${input.baseName}" }

  input:
  tupel val(input), path(sequences)

  output:
  path("*.csv")

  """
  wget -O data.tsv.gz ${input}
  gzip -d data.tsv.gz
  {
    head -1 data.tsv
    grep -v '^species' data.tsv | sort -k 2,2
  } > sorted.tsv
  rnac tarbase parse sorted.tsv ${sequences} .
  """
}

workflow tarbase {
  emit: data
  main:
    Channel.fromPath("files/import-data/tarbase/query.sql") | set { query }

    query | fetch_mirbase | set { sequences }

    Channel.of(params.tarbase.remotes) \
    | combine(sequences) \
    | process \
    | flatten \
    | set { data }
}
