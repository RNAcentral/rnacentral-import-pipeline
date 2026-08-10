process lookup {
  errorStrategy 'retry'
  maxRetries 3

  input:
  // CSV (legacy) or Parquet (writer_format=parquet) emissions of terms; the
  // ontology lookup is text-based, so decode parquet to CSV before merging.
  // Brace globs are illegal as *staging* names: nextflow substitutes the index
  // but leaves {csv,parquet} literal, and the shell expands it into two ln
  // targets. Use the active format, as everywhere else.
  path("terms*.${params.writer_format}")

  output:
  path("ontology_terms.${params.writer_format}")

  script:
  """
  set -o pipefail

  : > merged-terms.csv
  for f in terms*.csv; do
    [ -e "\$f" ] && cat "\$f" >> merged-terms.csv
  done
  for f in terms*.parquet; do
    [ -e "\$f" ] && parquet-to-csv "\$f" >> merged-terms.csv
  done
  sort -u merged-terms.csv >> unique-terms.txt
  rnac ols lookup-terms unique-terms.txt ontology_terms.${params.writer_format}
  """
}


workflow batch_lookup_ontology_information {
  take: term_ids

  main:
    term_ids | collect | lookup | set { term_info }
  emit: term_info
}
