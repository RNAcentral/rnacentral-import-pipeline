process lookup {
  errorStrategy 'retry'
  maxRetries 3

  input:
  // CSV (legacy) or Parquet (writer_format=parquet) emissions of terms; the
  // ontology lookup is text-based, so decode parquet to CSV before merging.
  path('terms*.{csv,parquet}')

  output:
  path('ontology_terms.csv')

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
  rnac ols lookup-terms unique-terms.txt ontology_terms.csv
  """
}


workflow batch_lookup_ontology_information {
  take: term_ids
  emit: term_info

  main:
    term_ids | collect | lookup | set { term_info }
}
