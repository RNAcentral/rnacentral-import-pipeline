#!/usr/bin/env nextflow

process quickgo_annotations {
  output:
  file('annotations.csv') into quickgo_annotations
  file('publications.csv') into quickgo_publications
  file('publication_mappings.csv') into quickgo_pub_maps
  file('terms.csv') into quickgo_terms

  script:
  """
  fetch "${params.go_annotations.quickgo.remote}" "${params.go_annotations.quickgo.pattern}"
  zcat "${params.go_annotations.quickgo.pattern}" | rnac ontologies quickgo -
  """
}

process pombase_annotations {
  input:
  file(mapping) from Channel.fromPath('files/import-ontologies/mappings/pombase.sql')

  output:
  file('annotations.csv') into pombase_annotations
  file('publications.csv') into pombase_publications
  file('publication_mappings.csv') into pombase_pub_maps
  file('terms.csv') into pombase_terms

  script:
  """
  psql -f $mapping "$PGDATABASE" > mapping.json
  fetch "${params.go_annotations.pombase.remote}" "${params.go_annotations.pombase.pattern}"
  zcat "${params.go_annotations.pombase.pattern}" | rnac ontologies pombase - mapping.json
  """
}

process rfam_ontology_terms {
  input:
  file(sql) from Channel.fromPath('files/ontologies/rfam-ontology-terms.sql')

  output:
  file('terms.csv') into rfam_terms
  file('rfam_mappings.csv') into rfam_mappings

  """
  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $sql > raw.tsv
  rnac ontologies rfam-terms raw.tsv
  """
}

Channel.empty()
  .mix(rfam_terms)
  .mix(quickgo_terms)
  .mix(pombase_terms)
  .set { terms }

process fetch_ontology_information {
  input:
  file('terms*.csv') from terms.collect()

  output:
  file('term-info.csv') into term_info

  """
  sort -u terms*.csv | rnac ontologies lookup-terms - term-info.csv
  """
}

Channel.empty()
  .mix(pombase_annotations)
  .mix(quickgo_annotations)
  .set { annotations }

Channel.empty()
  .mix(pombase_publications)
  .mix(quickgo_publications)
  .set { publications }

Channel.empty()
  .mix(pombase_pub_maps)
  .mix(quickgo_pub_maps)
  .set { pub_maps }

process pgload_ontology_annotations {
  echo true

  input:
  file('term-info.csv') from term_info.collect()
  file('annotations*.csv') from annotations.collect()
  file('publications*.csv') from publications.collect()
  file('publication_mappings*.csv') from pub_maps.collect()
  file('rfam_ontology_mappings*.csv') from rfam_mappings.collect()

  file(terms_ctl) from Channel.fromPath('files/ontologies/terms.ctl')
  file(ann_ctl) from Channel.fromPath('files/ontologies/annotations.ctl')
  file(pub_ctl) from Channel.fromPath('files/ontologies/publications.ctl')
  file(pub_map_ctl) from Channel.fromPath('files/ontologies/pub-map.ctl')
  file(rfam_map_ctl) from Channel.fromPath('files/ontologies/rfam-map.ctl')

  """
  find . -name '*.ctl' | xargs -I {} basename {} | xargs -I {} cp {} local_{}

  pgloader --on-error-stop local_$terms_ctl
  pgloader --on-error-stop local_$pub_ctl
  pgloader --on-error-stop local_$ann_ctl
  pgloader --on-error-stop local_$pub_map_ctl
  pgloader --on-error-stop local_$rfam_map_ctl
  """
}
