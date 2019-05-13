#!/usr/bin/env nextflow

pattern_info = [
  ['ensembl', /ENS\w+[0-9]+(\.[0-9]+)?/, file('files/text-mining/find-ensembl.sql')],
  ['lncipedia', /LINC[0-9]{5}/, file('files/text-mining/find-lncipedia.sql')],
  ['mirbase', /\w{3}-let-7\w-[3|5]p/, file('files/text-mining/find-mirbase.sql')],
  ['mirbase', /\w{3}-mir-[0-9]+\w?/, file('files/text-mining/find-mirbase.sql')],
  ['mirbase', /mir-[0-9]+\w?/, file('files/text-mining/find-mirbase.sql')],
  ['gtrnadb', /tRNA-\w{3}-\w{3}-[0-9]+-[0-9]+/, file('files/text-mining/find-gtrnadb.sql')],
]

name_info = [
  ['hgnc', file('files/text-mining/find-hgnc.sql')],
]

process find_known_publications {
  input:
  file(query) from Channel.fromPath('files/text-mining/known-publications.sql')

  output:
  file('known') into found_publications

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > known
  """
}

process find_known_items {
  input:
  set val(name), val(pattern), file(query) from Channel.from(pattern_info)

  output:
  set val(name), file('known') into known_items

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > known
  """
}

process fetch_raw_publications {
  memory 2.GB

  input:
  file(remote) from Channel.fromPath('/nfs/ftp/pub/databases/pmc/manuscripts/PMC*.txt.tar.gz')

  output:
  file('PMC*') into publication_files mode flatten

  """
  cp $remote batch
  tar xf batch
  """
}


publication_files
  .into { for_patterns; for_names }

for_patterns
  .combine(pattern_info)
  .map { pubs, name, pattern, query -> [name, pubs, pattern] }
  .set { to_search }

for_names
  .combine(name_info)
  .map { pubs, name, query -> [name, pubs, query, query] }
  .set { names_to_search }

process find_names { 
  tag { name + ':' + pubs.getName() }
  errorStrategy 'ignore'

  input:
  set val(name), file(pubs), file(names) from names_to_search

  output:
  set val(name), file('matches') into found_names

  """
  {
  echo match
  find -L $pubs -name '*.txt' |\
  xargs -I {} grep -Howf $names {} |\
  sort -fu |\
  sed 's/:/,/'
  } > matches
  """
}

process find_matches {
  tag { name + ':' + pubs.getName() }
  errorStrategy 'ignore'

  input:
  set val(name), file(pubs), val(pattern) from to_search

  output:
  set val(name), file('matches') into matches

  """
  find -L $pubs -name '*.txt' |\
  xargs -I {} grep -HoiE '$pattern' {} |\
  sort -fu |\
  sed 's/:/,/' > matches
  """
}

matches
  .groupTuple()
  .join(known_items)
  .set { merged_matches }

process merge_matches {
  tag { name }

  input:
  set val(name), file('matches*'), file(possible) from merged_matches

  output:
  set val(name), file('selected-matches') into selected_matches

  """
  set -o pipefail

  {
    echo 'filename,match'
    find . -name 'matches*' | xargs cat
  } > all-matches

  xsv index all-matches
  if [[ "$name" != "mirbase" ]]; then
    xsv join 1 $possible 2 all-matches | xsv select 2,3 > selected-matches
  else
    cp all-matches selected-matches
  fi
  """
}

selected_matches
  .mix(found_names.groupTuple())
  .combine(Channel.fromPath('/nfs/ftp/pub/databases/pmc/manuscripts/filelist.csv'))
  .combine(found_publications)
  .set { to_count }

process find_new_publications {
  tag { name }
  publishDir "$baseDir/text-mining/$name"

  input:
  set val(name), file(selected), file(filename_mapping), file(known_publications) from to_count

  output:
  set val(name), file('new-publications'), file('counts') into __matches

  """
  set -o pipefail

  # Fix the file listing (.xml -> .txt)
  sed 's/.xml/.txt/' $filename_mapping > files

  # Create the required indexes
  xsv index $known_publications
  xsv index files
  xsv index $selected

  # Join matches to PMIDs
  xsv join filename $selected File files | xsv select PMID,match > matched-publications

  # Remove all known PMIDs
  xsv join --left PMID matched-publications pmid $known_publications | xsv search -s 3 '^\$' | xsv select PMID,match > new-publications

  # Count the number of new PMIDs
  xsv select PMID new-publications |\
  sort -u |\
  wc -l > counts
  """
}
