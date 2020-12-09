#!/usr/bin/env nextflow

create_memory = params.sequence_search.create_fasta.memory_table

process find_db_to_export {
  when { params.sequence_search_export.run }

  input:
  path(query)

  output:
  path('dbs.txt')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > dbs.txt
  """
}

process query_database {
  tag { name }
  maxForks params.sequence_search.max_forks

  input:
  tuple val(name), path(query), val(parameters)

  output:
  tuple val(name), path('raw.json')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" $parameters "$PGDATABASE" > raw.json
  """
}

process create_fasta {
  tag { name }
  memory { create_memory.get(name.replaceAll("-", "_"), create_memory.__default) }

  input:
  tuple val(name), path(json)

  output:
  path("splits/$name*.fasta")
  path("${name}.hash")
  path("${name}.seqstat")

  script:
  def ordered = "${name}-ordered.fasta"
  """
  json2fasta ${json} - | rnac ftp-export sequences valid-nhmmer - ${ordered}
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
  path(fasta)
  path(hash)
  path(stats)

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

  tar -czvf $compressed ${fasta} ${stats}
  sort $hash | md5sum - > ${compressed}.hash
  scp ${compressed} $remote
  scp ${compressed}.hash $remote
  """
}

workflow sequence_search_export {
  Channel.fromPath('files/ftp-export/sequences/databases.sql') | set { db_query }
  Channel.fromPath('files/ftp-export/sequences/database-specific.sql') | set { db_specific_query }

  Channel.fromPath('files/sequence-search-export/*.sql') \
  | filter { params.sequence_search_export.run } \
  | map { fn -> [file(fn).baseName, fn, ''] } \
  | set { simple_queries }

  db_query \
  | find_db_to_export \
  | splitCsv() \
  | combine(db_specific_query)
  | map { db, query -> [db, query, "-v db='%${db}%'"] }
  | mix(simple_queries)
  | map { db, query, param -> [db.toLowerCase().replace(' ', '_'), query, param] }
  | map { db, q, p -> [(db == "tmrna_website" ? "tmrna_web" : db), q, p] } \
  | query_data \
  | create_fasta

  create_fasta.out.sequences | collect | set { sequences }
  create_fasta.out.hashes | collect | set { hashes }
  create_fasta.out.stats | collect | set { stats }

  atomic_publish(sequences, hashes, stats)
}

workflow {
  sequence_search_export()
}
