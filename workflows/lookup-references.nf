process merge_and_split_all_publications {
  input:
  // Inputs may arrive as either CSV (legacy) or Parquet (writer_format=parquet),
  // depending on what the upstream parsers emitted. The downstream lookup
  // stays text-based, so we decode any parquet inputs to CSV before merging.
  path("ref_ids*.{csv,parquet}")

  output:
  path('split-refs/*.csv')

  """
  set -o pipefail

  mkdir split-refs
  : > all-ids
  for f in ref_ids*.csv; do
    [ -e "\$f" ] && cat "\$f" >> all-ids
  done
  for f in ref_ids*.parquet; do
    [ -e "\$f" ] && parquet-to-csv "\$f" >> all-ids
  done
  split --additional-suffix=".csv" --number l/${params.lookup_publications.maxForks} all-ids split-refs/refs
  """
}

process fetch_publications {
  when { params.needs_publications }
  queue 'datamover'
  memory 8.GB
  container ''

  output:
  path('out')

  """
  cp /nfs/ftp/public/databases/pmc/PMCLiteMetadata/PMCLiteMetadata.tgz .
  tar xvf PMCLiteMetadata.tgz
  """
}

process lookup_publications {
  memory 6.GB
  maxForks params.lookup_publications.maxForks

  input:
  tuple path(refs), path(pubs)

  output:
  path("references.${params.writer_format}")

  script:
  def out = "references.${params.writer_format}"
  """
  rnac europepmc stream-lookup --ignore-missing --allow-fallback $pubs $refs $out
  """
}

workflow lookup_ref_ids {
  take: ref_ids
  emit: publications

  main:
    ref_ids \
    | collect \
    | merge_and_split_all_publications \
    | flatten \
    | combine(fetch_publications()) \
    | lookup_publications \
    | set { publications }
}
