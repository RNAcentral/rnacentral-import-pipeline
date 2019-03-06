#!/usr/bin/env nextflow

pattern_info = [
  ['ensembl', /ENS[A-Z]+\d+(\.\d+)?/ 'files/text-mining/find-ensembl.sql'],
  /* ['lncipedia', /LINC\d{5}/: 'files/text-mining/find-lncipedia.sql'], */
  /* ['mirbase', /\w{3}-let-7\w-[3|5]p/: 'files/text-mining/find-mirbase.sql'], */
  /* ['mirbase', /\w{3}-mir-\d+\w?/: 'files/text-mining/find-mirbase.sql'], */
  /* ['gtrnadb', /\w{3}-\w{3}\d+-\d+/: 'files/text-mining/find-gtrnadb.sql'], */
]

remote_files = Channel.fromPath('/nfs/ftp/pub/databases/pmc/manuscripts/PMC*.txt.tar.gz')

process fetch_raw_publications {
  input:
  file(remote) from remote_files

  output:
  file('*.txt') into publication_files

  """
  cp $remote batch
  tar xf batch
  """
}

publication_files
  .combine(pattern_info)
  .map { pub, name, pattern, query -> [name, pub, pattern] }
  .set { to_search }

process find_matches {
  input:
  set val(name), file(publication), val(pattern) from to_search

  output:
  set val(name), file('matches') into matching_patterns

  """
  grep -Hoiw $pattern $publication | sort -fu | sed 's/:/,/' > matches
  """
}

process find_known {
  input:
  set val(name), val(pattern), file(query) from Channel.from(pattern_info)

  output:
  set val(name), file('known') into known_items

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > known
  """
}

counts
  .collect()
  .combine(known_items)
  .set { merged_counts }

process merge_matches {
  input:
  set val(name), file(possible), file('matches*') from merged_counts

  output:
  set val(name), file('selected-matches') to selected_matches

  """
  set -o pipefail

  {
    echo 'filename,match'
    find . -name 'matches*' | xargs cat
  } > all-matches

  xsv index all-matches
  xsv join 1 $possible 2 all-matches | xsv select 2,3 > selected-matches
  """
}

process find_known_publications {
  input:
  file(query) from known_pubs_query

  output:
  file('known') into found_publications

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > known
  """
}

selected_matches
  .combine(file('/nfs/ftp/pub/databases/pmc/manuscripts/filelist.csv'))
  .join(found_publications)
  .set { to_count }

process find_new_publications {
  publish "$baseDir/text-mining/$name"

  input:
  set val(name), file(known_publications), file(filename_mapping), file(selected) from to_count

  output:
  set val(name), file('new-publications'), file('counts') into __matches

  """
  set -o pipefail

  sed 's|PMC[0-9]\\+XXXXX/||' $known_publications > pubs

  xsv index pubs
  xsv join File pubs filename $selected |\
  xsv select PMID,match |\
  tee new-publications |\
  xsv select PMID |\
  sort -u |\
  wc -l > counts
  """
}
