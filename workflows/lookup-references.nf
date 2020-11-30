process merge_and_split_all_publications {
  input:
  path("ref_ids*.csv")

  output:
  path('split-refs/*.csv')

  """
  set -o pipefail

  mkdir split-refs
  find . -name 'ref_ids*.csv' | xargs cat > all-ids
  split --additional-suffix=".csv" --number l/${params.lookup_publications.maxForks} all-ids split-refs/refs
  """
}

process fetch_publications {
  when { params.will_import_data }

  output:
  path('out')

  """
  curl -L http://europepmc.org/ftp/pmclitemetadata/PMCLiteMetadata.tgz
  tar xvf PMCLiteMetadata.tgz
  """
}

process lookup_publications {
  memory 4.GB
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
  emit: publications

  main:
    ref_ids \
    | collect() \
    | merge_and_split_all_publications \
    | flatten \
    | combine(fetch_publications()) \
    | lookup_publications \
    | set { publications }
}
