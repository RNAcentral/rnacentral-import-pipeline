#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_db_to_export {
  input:
  path(query)

  output:
  path('dbs.txt')

  when: params.export?.sequence_search?.run

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > dbs.txt
  """
}

process query_database {
  tag { "$name" }
  maxForks params.export.sequence_search.max_forks

  input:
  tuple val(name), path(query), val(partition)

  output:
  tuple val(name), path('raw.json')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    $partition \
    -f "$query" \
    "\$PGDATABASE" > raw.json
  """
}

process create_fasta {
  tag { name }
  memory { def cm = params.export.sequence_search.create_fasta.memory_table; cm.get(name.replaceAll("-", "_"), cm.__default) }

  input:
  tuple val(name), path(json)

  output:
  path "splits/$name*.fasta", emit: sequences, optional: true
  path "${name}.hash", emit: hashes
  path "${name}.seqstat", emit: stats

  script:
  def ordered = "${name}-ordered.fasta"
  """
  set -o pipefail
  json2fasta.py ${json} - | rnac ftp-export sequences valid-nhmmer - ${ordered}
  md5sum ${ordered} > ${name}.hash
  cp ${ordered} ${name}.fasta
  mkdir -p splits
  if [ -s ${name}.fasta ]; then
    esl-seqstat --dna ${name}.fasta > ${name}.seqstat
    split-sequences \
      --max-file-size ${params.export.sequence_search.max_file_size} \
      ${name}.fasta splits/
  else
    echo "WARNING: ${name}.fasta is empty - no nhmmer-valid sequences matched; skipping seqstat/split" >&2
    : > ${name}.seqstat
  fi
  """
}

process atomic_publish {
  stageInMode 'copy'
  queue 'datamover'

  input:
  path(fasta)
  path(hash)
  path(stats)

  output:
  val('done')

  script:
  def publish = params.export.sequence_search.publish
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

workflow sequence_search {
  take: _flag
  main:
    channel.fromPath('files/ftp-export/sequences/databases.sql') | set { db_query }
    channel.fromPath('files/ftp-export/sequences/database-specific.sql') | set { db_specific_query }

    channel.fromPath('files/sequence-search-export/*.sql') \
    | filter { params.export.sequence_search.run } \
    | map { fn -> [file(fn).baseName, fn, ''] } \
    | set { simple_queries }

    find_db_to_export(db_query) \
    | splitCsv \
    | combine(db_specific_query) \
    | map { db, full_descr, query -> [db, query, "-v db='%${full_descr}%'"] } \
    | mix(simple_queries) \
    | map { db, query, param -> [db.toLowerCase().replace(' ', '_').replace('/', '_'), query, param] } \
    | map { db, q, p -> [(db == "tmrna_website" ? "tmrna_web" : db), q, p] } \
    | query_database \
    | create_fasta

    create_fasta.out.sequences | collect | set { sequences }
    create_fasta.out.hashes | collect | set { hashes }
    create_fasta.out.stats | collect | set { stats }

    atomic_publish(sequences, hashes, stats) | set { done }
  emit: done
}

workflow {
  sequence_search(channel.of('ready'))
}
