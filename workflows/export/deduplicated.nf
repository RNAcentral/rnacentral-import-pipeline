#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// The deduplicated set (tier-2 representatives of the sequence clustering) is a subset of
// the active sequences, so it is cut from the shared active dump rather than re-queried.
// Downstream, the FASTA and the HuggingFace split use the same writers as the full set.

process fetch_deduplicated_urs {
  input:
    path(query)

  output:
    path('deduplicated_urs.txt')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > deduplicated_urs.txt
  """
}

process filter_dump {
  memory '16 GB'

  input:
    path(active_json)
    path(urs)

  output:
    path('deduplicated_sequences.json'), emit: json

  script:
  """
  set -euo pipefail
  export PYTHONIOENCODING=utf8
  rnac ftp-export sequences filter-dump $active_json $urs deduplicated_sequences.json
  """
}

process deduplicated_fasta {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
    path(dedup_json)

  output:
    path('rnacentral_active_deduplicated.fasta.gz')

  script:
  """
  set -euo pipefail
  export PYTHONIOENCODING=utf8
  json2fasta.py $dedup_json - | gzip > rnacentral_active_deduplicated.fasta.gz
  """
}

workflow deduplicated_sequences {
  take:
    active_json
  main:
    if (params.clustering.export) {
      urs = fetch_deduplicated_urs(channel.fromPath('files/ftp-export/sequences/deduplicated-urs.sql'))
      filter_dump(active_json, urs)
      deduplicated_fasta(filter_dump.out.json)
      dedup_json = filter_dump.out.json
    } else {
      dedup_json = channel.empty()
    }
  emit:
    dedup_json
}
