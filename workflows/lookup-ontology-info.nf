process lookup {
  input:
  path('terms*.csv')

  output:
  path('ontology_terms.csv')

  script:
  """
  set -o pipefail

  find . -name 'terms*.csv' | xargs cat | sort -u >> unique-terms.txt
  rnac ols lookup-terms unique-terms.txt ontology_terms.csv
  """
}


workflow batch_lookup_ontology_information {
  take: term_ids
  emit: term_info

  main:
    term_ids | collect | lookup | set { term_info }
}
