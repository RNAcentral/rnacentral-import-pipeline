#!/usr/bin/env nextflow

create_memory = params.sequence_search.create_fasta.memory_table

Channel.fromPath('files/sequence-search-export/*.sql')
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
  .mix(simple_queries)
  .map { db, query, param -> [db.toLowerCase().replace(' ', '_'), query, param] }
  .map { db, q, p -> [(db == "tmrna_website" ? "tmrna_web" : db), q, p] }
  .set { queries }

process query_database {
  tag { name }
  maxForks params.sequence_search.max_forks

  input:
  set val(name), file(query), val(parameters) from queries

  output:
  set val(name), file('raw.json') into for_fasta

  """
  psql -v ON_ERROR_STOP=1 -f "$query" $parameters "$PGDATABASE" > raw.json
  """
}

process create_fasta {
  tag { name }
  memory { create_memory.get(name.replaceAll("-", "_"), create_memory.__default) }

  input:
  set val(name), file(json) from for_fasta

  output:
  file("splits/$name*.fasta") into sequences
  file("${name}.hash") into hashes

  script:
  def ordered = "${name}-ordered.fasta"
  """
  json2fasta.py ${json} ${ordered}
  md5sum ${ordered} > ${name}.hash
  seqkit shuffle --two-pass ${ordered} > ${name}.fasta
  esl-seqstat ${name}.fasta > ${name}.seqstat
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
  def sequences = 'sequences-database.fa'
  def compressed = "${sequences}.tar.gz"
  """
  if [[ -n "$publish.host" ]]; then
    ssh "$publish.host" 'mkdir -p $publish.path' || true
  else
    mkdir -p "$publish.path" || true
  fi

  tar -czvf $compressed ${fasta}
  sort $hash | md5sum - > ${compressed}.hash
  scp ${compressed} $remote
  scp ${compressed}.hash $remote
  """
}
