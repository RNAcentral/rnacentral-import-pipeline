#!/usr/bin/env nextflow

process quickgo_annotations {
  output:
  file('annotations.csv') into annotations
  file('publications.csv') into publications
  file('pub_med.csv') into pub_maps
  file('terms.csv') into terms

  """
  gzip -d "${params.quickgo.remote}" > quickgo.gpa
  rnac ontologies quickgo quickgo.gpa
  """
}

process fetch_ontology_information {
  input:
  file('terms*.csv') from terms.collect()

  output:
  file('term-info.csv') into term_info

  """
  sort -u terms*.csv | rnac ontologies lookup-terms - term-info.csv
  """
}

process pgload_ontology_annotations {
  input:
  file('term-info.csv') from term_info.collect()
  file('annotations*.csv') from annotations.collect()
  file('publications*.csv') from publications.collect()
  file('pub_map*.csv') from pub_map.collect()

  file(terms_ctl) from Channel.fromPath('files/ontologies/terms.ctl')
  file(ann_ctl) from Channel.fromPath('files/ontologies/annotations.ctl')
  file(pub_ctl) from Channel.fromPath('files/ontologies/publications.ctl')
  file(pub_map_ctl) from Channel.fromPath('files/ontologies/pub-map.ctl')

  """
  find . -name '.ctl' | xargs -I {} cp {} local_{}

  pgloader local_$terms_ctl
  pgloader local_$pub_ctl
  pgloader local_$ann_ctl
  pgloader local_$pub_map_ctl
  """
}
