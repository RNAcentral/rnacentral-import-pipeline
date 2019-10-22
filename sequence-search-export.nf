#!/usr/bin/env nextflow

Channel.fromPath('files/sequence-search/*.sql')
  .map { fn -> [file(fn).baseName, file(fn), ''] }
  .set { simple_queries }

process find_db_to_export {
  input:
  file query from Channel.fromPath('files/ftp-export/sequences/databases.sql')

  output:
  stdout into raw_dbs

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE"
  """
}

raw_dbs
  .splitCsv()
  .combine(Channel.fromPath('files/ftp-export/sequences/database-specific.sql'))
  .map { db, query -> [db, query, "-v db='%${db}%'"] }
  .mix { simple_queries }
  .set { queries }

process query_database {
  tag { name }
  maxForks params.sequence_search.query.max_forks

  input:
  set val(name), file(query), val(parameters) from db_to_export

  output:
  set val(name), file('raw.json') into for_fasta

  """
  psql -v ON_ERROR_STOP=1 -f "$query" $parameters "$PGDATABASE" > raw.json
  """
}

process create_fasta {
  input:
  set val(name), file(json) from for_fasta

  output:
  file("splits/$name*.fasta") into sequences
  file("${name}.hash") into hashes

  ordered = "${name}-ordered.fasta"
  """
  json2fasta.py ${json} ${ordered}
  md5hash ${ordered} > ${name}.hash
  seqkit shuffle --two-pass ${ordered} > ${name}.fasta
  split-sequences \
    --max-file-size ${params.sequence_search.max_file_size} \
    ${name}.fasta splits/
  """
}

process atomic_publish {
  input:
  file(fasta) from sequences.collect()
  file(hash) from hashes.collect()

  script:
  def publish = params.sequence_search.publish
  def remote = (publish.host ? "$publish.host:" : "") + publish.path
  def sequences = 'sequences'
  def compressed = "${sequences}.tar.gz"
  """
  if [[ -n "$publish.host" ]]; then
    ssh "$publish.host" 'mkdir -p $publish.path' || true
  else
    mkdir -p "$publish.path" || true
  fi

  mkdir ${sequences}
  cp ${fasta} ${sequences}
  tar -czvf $compressed $sequences
  cat $hash | md5sum - > ${compressed}.hash
  scp ${compressed} $remote
  """
}
