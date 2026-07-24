process merge_and_split_all_publications {
  input:
  path("ref_ids*.csv")

  output:
  path('split-refs/*.csv')

  script:
  """
  set -o pipefail

  mkdir split-refs
  find . -name 'ref_ids*.csv' | xargs cat > all-ids
  split --additional-suffix=".csv" --number l/${params.lookup_publications.maxForks} all-ids split-refs/refs
  """
}

process fetch_publications {
  queue 'datamover'
  memory 8.GB
  container ''

  output:
  path('out')

  when: { params.get('needs_publications', false) }

  script:
  """
  cp /nfs/ftp/public/databases/pmc/PMCOALiteMetadata/PMCOALiteMetadata.tgz .
  tar xvf PMCOALiteMetadata.tgz
  """
}

process lookup_publications {
  memory 6.GB
  maxForks params.lookup_publications.maxForks

  input:
  tuple path(refs), path(pubs)

  output:
  path("references.csv")

  script:
  """
  rnac europepmc stream-lookup --ignore-missing --allow-fallback $pubs $refs references.csv
  """
}

workflow lookup_ref_ids {
  take: ref_ids

  main:
    if (params.get('needs_publications', false)) {
      ref_ids \
      | collect \
      | merge_and_split_all_publications \
      | flatten \
      | combine(fetch_publications()) \
      | lookup_publications \
      | set { publications }
    } else {
      channel.empty() | set { publications }
    }

  emit: publications
}