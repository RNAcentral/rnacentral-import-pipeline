#!/usr/bin/env nextflow

pattern_info = [
  ['ensembl', /ENS[A-Z]+\d+(\.\d+)?/ 'files/text-mining/find-ensembl.sql'],
  /* ['lncipedia', /LINC\d{5}/: 'files/text-mining/find-lncipedia.sql'], */
  /* ['mirbase', /\w{3}-let-7\w-[3|5]p/: 'files/text-mining/find-mirbase.sql'], */
  /* ['mirbase', /\w{3}-mir-\d+\w?/: 'files/text-mining/find-mirbase.sql'], */
  /* ['gtrnadb', /\w{3}-\w{3}\d+-\d+/: 'files/text-mining/find-gtrnadb.sql'], */
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
  input:
  file(remote) from Channel.fromPath('/nfs/ftp/pub/databases/pmc/manuscripts/PMC*.txt.tar.gz')

  output:
  file('**/*.txt') into publication_files

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
  set val(name), file('matches') into matches

  """
  grep -Howi $pattern $publication | sort -fu | sed 's/:/,/' > matches
  """
}

matches
  .groupTuple()
  .join(known_items)
  .set { merged_matches }

process merge_matches {
  input:
  set val(name), file('matches*'), file(possible) from merged_matches

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

  # Create the required indexes
  sed 's|PMC[0-9]\\+XXXXX/||' $filename_mapping > files
  xsv index $known_publications
  xsv index files

  # Join matched publications to known ones
  xsv join File files filename $selected |\
  xsv select PMID,match > matched-publications

  # Remove all known publications
  xsv join --left PMID matched-publications pmid $known_publications |\
  xsv search -s 3 '^$' |\
  xsv select 1,2 > new-publications

  # Count the number of publications
  xsv select PMID |\
  sort -u |\
  wc -l > counts
  """
}
