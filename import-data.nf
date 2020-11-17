include { lookup_ref_ids } from './workflows/lookup-references'
include { lookup_ontology_info } from './workflows/lookup-ontology-info'
include { parse_databases } from './workflows/parse-databases'
include { parse_metadata } from './workflows/parse-metadata'
include { load_data } from './workflows/load-data'

workflow import_data {
  emit: post_release
  main:
    parse_data() \
    | mix(parse_metadata()) \
    | branch {
      terms: it.name == "terms.csv"
      ref_ids: it.name == "ref_ids.csv"
      csv: true
    } \
    | set { results }

    results.out.terms | lookup_ontology_terms | set { term_info }
    results.out.ref_ids | lookup_ref_ids | set { references }

    results.csv \
    | mix(term_info, references) \
    | load_data \
    | set { post_release }
}

workflow {
  import_data()
}
