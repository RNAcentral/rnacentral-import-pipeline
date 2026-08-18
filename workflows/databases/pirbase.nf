process find_urls {
  executor 'local'

  output:
  path("urls.txt")

  script:
  """
  rnac pirbase urls $params.databases.pirbase.remote urls.txt
  """
}

process find_known {
  memory '8GB'

  input:
  path(query)

  output:
  path('known.sqlite')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query \$PGDATABASE > md5s
  rnac pirbase build-known md5s known.sqlite
  """
}

// The full sets are filtered to sequences RNAcentral already holds; without
// that they are tens of millions of piRNAs.
process parse_full {
  tag { code }
  memory '5GB'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(code), val(url), path(known)

  output:
  path('*.{csv,parquet}')

  script:
  """
  wget --timeout=60 --tries=10 --continue -O data.fa.gz '$url'
  gzip -df data.fa.gz
  rnac pirbase parse --known $known $code data.fa .
  """
}

// The gold standard sets are curated and small, so they are imported whole.
// They must not be filtered: their ids barely overlap the full sets, so a new
// gold sequence has nothing to match against.
process parse_gold {
  tag { "${code} gold" }
  memory '2GB'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(code), val(url)

  output:
  path('*.{csv,parquet}')

  script:
  """
  wget --timeout=60 --tries=10 --continue -O gold.fa.gz '$url'
  gzip -df gold.fa.gz
  rnac pirbase parse $code gold.fa .
  """
}

workflow pirbase {
  main:
    if ( params.databases.pirbase?.run ) {
      channel.fromPath('files/import-data/pirbase/known-md5.sql') \
      | find_known \
      | set { known }

      find_urls \
      | splitCsv \
      | branch { row ->
        full: row[1] == 'full'
        gold: row[1] == 'gold'
      } \
      | set { urls }

      urls.full \
      | map { row -> tuple(row[0], row[2]) } \
      | combine(known) \
      | parse_full \
      | set { from_full }

      urls.gold \
      | map { row -> tuple(row[0], row[2]) } \
      | parse_gold \
      | set { from_gold }

      from_full \
      | mix(from_gold) \
      | flatten \
      | set { data_files }
    }
    else {
      channel.empty() | set { data_files }
    }

  emit: data_files
}
