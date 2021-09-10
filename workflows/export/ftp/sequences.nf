process active {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'copy'
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
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta.py - rnacentral_active.fasta
  head rnacentral_active.fasta > example.txt
  gzip rnacentral_active.fasta
  cp template.txt readme.txt
  """
}

process inactive {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'copy'
  when: params.ftp_export.sequences.inactive.run

  input:
  path(query)

  output:
  path('rnacentral_inactive.fasta.gz')

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta.py - - | gzip > rnacentral_inactive.fasta.gz
  """
}

process species_specific {
  when: params.ftp_export.sequences.species.run

  input:
  path(query)

  output:
  path('rnacentral_species_specific_ids.fasta'), emit: sequences

  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | json2fasta.py - rnacentral_species_specific_ids.fasta 
  """
}

process create_ssi {
  publishDir "${params.ftp_export.publish}/sequences/.internal/", mode: 'copy'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.ssi')

  """
  esl-sfetch --index rnacentral_species_specific_ids.fasta
  """
}

proecss compress_species_fasta {
  publishDir "${params.ftp_export.publish}/sequences/", mode: 'copy'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.gz')

  """
  gzip -k rnacentral_species_specific_ids.fasta
  """
}


process find_dbs {
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
  publishDir "${params.ftp_export.publish}/sequences/by-database", mode: 'copy'
  when: params.ftp_export.sequences.by_database.run

  input:
  tuple val(db), path(query)

  output:
  file('*.fasta')

  script:
  """
  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" -v db='%${db}%' "$PGDATABASE" > raw.json
  json2fasta.py raw.json ${db.toLowerCase().replaceAll(' ', '_').replace('/', '_')}.fasta
  """
}

workflow fasta_export {
  Channel.fromPath('files/ftp-export/sequences/active.sql') | set { active_sql }
  Channel.fromPath('files/ftp-export/sequences/readme.txt') | set { readme }
  active(active_sql, readme)

  Channel.fromPath('files/ftp-export/sequences/inactive.sql') | inactive
  Channel.fromPath('files/ftp-export/sequences/species-specific.sql') \
  | species_specific \
  | (compress_species_fasta & create_ssi)

  Channel.fromPath('files/ftp-export/sequences/databases.sql') \
  | find_dbs \
  | splitCsv \
  | combine(Channel.fromPath('files/ftp-export/sequences/database-specific.sql')) \
  | database_specific
}
