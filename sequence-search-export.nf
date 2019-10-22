#!/usr/bin/env nextflow

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
  .set { db_sequences }

process database_specific_fasta {
  tag { db }
  maxForks params.sequence_search.query.max_forks

  input:
  set val(db), file(query) from db_sequences

  output:
  set val("${db}.fasta"), file('raw.json') into database_json

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" -v db='%${db}%' "$PGDATABASE" > raw.json
  """
}

process fetch_simple {
  tag { query.baseName }
  maxForks params.sequence_search.query.max_forks

  input:
  file(query) from Channel.fromPath('files/sequence-search/*.sql')

  output:
  set val("${query.baseName}.fasta"), file('raw.json') into simple_json

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  """
}

simple_json
  .mix(database_json)
  .set { for_fasta }

process create_fasta {
  input:
  set val(name), file(json) from for_fasta

  output:
  file("splits/$name*.fasta") into sequences
  file("${name}.hash") into hashes

  """
  json2fasta.py ${json} ${name}-ordered.fasta
  md5hash ${name}-ordered.fasta > ${name}.hash
  seqkit shuffle --two-pass ${name}-ordered.fasta > $name
  split-sequences \
    --max-file-size ${params.sequence_search.max_file_size} \
    ${name} splits/
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
