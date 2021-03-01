process active {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'move'
  when: params.ftp_export.sequences.active.run

  input:
  path(query)
  path('template.txt')

  output:
  path('rnacentral_active.fasta.gz'), emit: 'active'
  path('example.txt'), emit: 'example'
  path('readme.txt'), emit: 'readme'

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta - rnacentral_active.fasta
  head rnacentral_active.fasta > example.txt
  gzip rnacentral_active.fasta
  cp template.txt readme.txt
  """
}

process inactive {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'move'
  when: params.ftp_export.sequences.inactive.run

  input:
  path(query)

  output:
  path('rnacentral_inactive.fasta.gz')

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta - - | gzip > rnacentral_inactive.fasta.gz
  """
}

process species_specific {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'move'
  when: params.ftp_export.sequences.species.run

  input:
  path(query)

  output:
  path('rnacentral_species_specific_ids.fasta.gz')

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta - - | gzip > rnacentral_species_specific_ids.fasta.gz
  """
}

process find_db_to_export {
  input:
  path(query)

  output:
  path('dbs.txt')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > dbs.txt
  """
}

process database_specific {
  tag { db }
  maxForks params.ftp_export.sequences.by_database.max_forks
  publishDir "${params.ftp_export.publish}/sequences/by-database", mode: 'move'
  when: params.ftp_export.sequences.by_database.run

  input:
  tuple val(db), path(query)

  output:
  file('*.fasta')

  script:
  """
  set -o pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" -v db='%${db}%' "$PGDATABASE" > raw.json
  json2fasta raw.json ${db.toLowerCase().replaceAll(' ', '_')}.fasta
  """
}

process extract_nhmmer_valid {
  publishDir "${params.ftp_export.publish}/sequences/.internal/", mode: 'move'
  when: params.ftp_export.sequences.nhmmer.run

  input:
  path(rna)

  output:
  path('rnacentral_nhmmer.fasta')

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  zcat $rna | rnac ftp-export sequences valid-nhmmer - rnacentral_nhmmer.fasta
  """
}

process extract_nhmmer_invalid {
  publishDir "${params.ftp_export.publish}/sequences/.internal/", mode: 'move'
  when: params.ftp_export.sequences.nhmmer.run

  input:
  path(rna)

  output:
  path('rnacentral_nhmmer_excluded.fasta')

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  zcat $rna | rnac ftp-export sequences invalid-nhmmer - rnacentral_nhmmer_excluded.fasta
  """
}


workflow fasta_export {
  Channel.fromPath('files/ftp-export/sequences/active.sql') | set { active_sql }
  Channel.fromPath('files/ftp-export/sequences/readme.txt') | set { readme }
  active(active_sql, readme)

  active.out.active | (extract_nhmmer_valid & extract_nhmmer_invalid)

  Channel.fromPath('files/ftp-export/sequences/inactive.sql') | inactive
  Channel.fromPath('files/ftp-export/sequences/species-specific.sql') | species_specific

  Channel.fromPath('files/ftp-export/sequences/databases.sql') \
  | find_dbs \
  | splitCsv \
  | combine(Channel.fromPath('files/ftp-export/sequences/database-specific.sql')) \
  | database_specific
}
