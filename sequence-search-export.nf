#!/usr/bin/env nextflow

nextflow.enable.dsl=2

create_memory = params.sequence_search.create_fasta.memory_table

process find_db_to_export {
  when { params.sequence_search.run }

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
  tuple val(name), val(partition) path(query)

  output:
  tuple val(name), path('raw.json')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v partition=$partition \
    -f "$query" \
    "$PGDATABASE" > raw.json
  """
}

process create_fasta {
  tag { name }
  memory { create_memory.get(name.replaceAll("-", "_"), create_memory.__default) }

  input:
  tuple val(name), path(json)

  output:
  path "splits/$name*.fasta", emit: sequences
  path "${name}.hash", emit: hashes
  path "${name}.seqstat", emit: stats

  script:
  def ordered = "${name}-ordered.fasta"
  """
  json2fasta.py ${json} - | rnac ftp-export sequences valid-nhmmer - ${ordered}
  md5sum ${ordered} > ${name}.hash
  cp ${ordered} ${name}.fasta
  esl-seqstat --dna ${name}.fasta > ${name}.seqstat
  split-sequences \
    --max-file-size ${params.sequence_search.max_file_size} \
    ${name}.fasta splits/
  """
}

process atomic_publish {
  stageInMode 'copy'

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
  | filter { params.sequence_search.run } \
  | map { fn -> [file(fn).baseName, fn, ''] } \
  | set { simple_queries }

  db_query \
  | find_db_to_export \
  | splitCsv \
  | combine(db_specific_query) \
  | map { db, query -> [db, query, "-v db='%${db}%'"] } \
  | mix(simple_queries) \
  | map { db, query, param -> [db.toLowerCase().replace(' ', '_').replace('/', '_'), query, param] } \
  | map { db, q, p -> [(db == "tmrna_website" ? "tmrna_web" : db), q, p] } \
  | query_database \
  | create_fasta

  create_fasta.out.sequences | collect | set { sequences }
  create_fasta.out.hashes | collect | set { hashes }
  create_fasta.out.stats | collect | set { stats }

  atomic_publish(sequences, hashes, stats)
}

workflow {
  sequence_search_export()
}
