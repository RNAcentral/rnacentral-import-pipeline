#!/usr/bin/env nextflow

def searches = Channel.from([
  'mirbase',
  'ensembl'
])

process get_all_references {
  memory { params.text_mining.get_all_references.directives.memory }

  input:
  val(remote) from Channel.from(params.metadata.europepmc.inputs.data.remote)

  output:
  file("references.db") into all_publications
  file('out/PMC.*.xml') into publication_metadata mode flatten

  """
  fetch europepmc '$remote' 'out'
  rnac europepmc index-xml out references.db
  """
}

process fetch_full_text {
  memory { params.text_mining.fetch_full_text.directives.memory }

  input:
  val(remote) from Channel.from(params.text_mining.fetch_full_text.inputs.remote)

  output:
  file('PMC*') into full_text_files mode flatten

  """
  fetch generic '$remote' ignored
  """
}

full_text_files
  .filter { d -> d.isDirectory() }
  .mix(publication_metadata)
  .into { pattern_publication_files; name_publication_files }

process fetch_names {
  tag { query.getBaseName() }

  input:
  file(query) from Channel.fromPath('files/text-mining/names/*.sql')

  output:
  set val("${query.getBaseName()}"), file('known') into known_names

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > known
  """
}

searches
  .combine(pattern_publication_files)
  .set { to_mine_patterns }

process find_pattern_matches {
  memory { params.text_mining.find_pattern_matches.directives.memory }
  tag { name + ':' + pubs.getName() }

  input:
  set val(name), file(pubs) from to_mine_patterns

  output:
  file('matches') into matching_pattern_publications

  """
  rnac text-mining $name $pubs matches
  """
}

known_names
  .combine(name_publication_files)
  .set { to_mine_names }

process find_name_matches {
  memory { params.text_mining.find_name_matches.directives.memory }
  tag { name + ':' + pubs.getName() }

  input:
  set val(name), file(names), file(pubs) from to_mine_names

  output:
  file('matches') into matching_name_publications

  """
  rnac text-mining $name $names $pubs matches
  """
}

Channel.empty()
  .mix(matching_name_publications, matching_pattern_publications)
  .filter { f -> !f.isEmpty() }
  .collect()
  .set { pmids }

process batch_find_pmids {
  memory { params.text_mining.batch_find_pmids.directives.memory }

  input:
  file('pubs*') from pmids
  file('references.db') from all_publications

  output:
  file('matched-references.csv') into matched_references

  """
  set -o pipefail

  find . -name 'pubs*' |\
  xargs -I {} sort -k 1,1 {} |\
  rnac europepmc lookup --column 0 --allow-fallback references.db - - >> matched-references.csv
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
