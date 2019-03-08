#!/usr/bin/env nextflow

def searches = [
  'mirbase',
  'ensembl'
]

process get_all_references {
  output:
  file("references.db") into all_publications

  """
  fetch europepmc 'http://europepmc.org/ftp/pmclitemetadata/PMCLiteMetadata.tgz' 'out'
  rnac europepmc index-xml out references.db
  """
}

process fetch_raw_publications {
  memory params.text_mining.known_publications.memory

  input:
  file(remote) from Channel.fromPath(params.text_mining.known_publications.remote)

  output:
  file('PMC*') into publication_files mode flatten

  """
  fetch generic '$remote' ignored
  """
}

publication_files
  .combine(searches)
  .map { pubs, name -> [name, pubs] }
  .set { to_mine }

process find_matches { 
  tag { name + ':' + pubs.getName() }

  input:
  set val(name), file(pubs) from to_mine

  output:
  set val(name), file('matches') into matching_publications

  """
  rnac text-mining $name $pubs matches
  """
}

matching_publications
  .map { n, p -> p }
  .collect()
  .set { pmids }

process batch_find_pmids {
  input:
  file('pubs*') from pmids
  file('references.db') from all_publications

  output:
  file('matched-references.csv') into matched_references

  """
  set -o pipefail

  find . -name 'pubs*' |\
  xargs -I {} sort -k,k 1 {} |\
  xargs -I {} rnac europepmc lookup --column 0 --allow-fallback references.db {} - >> matched-references.csv
  """
}

process load_matches {
  input:
  file('matched-references.csv') from matched_references
  file(ctl) from Channel.fromPath('files/text-mining/matched-references.ctl')

  """
  split-and-load $ctl 'matched-references*.csv' ${params.import_data.chunk_size} matched_references
  """
}
