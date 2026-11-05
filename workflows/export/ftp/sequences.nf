process active {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path(active_json)
  path('template.txt')

  output:
  path('rnacentral_active.fasta.gz'), emit: 'active'
  path('example.txt'), emit: 'example'
  path('readme.txt'), emit: 'readme'

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  json2fasta.py "$active_json" rnacentral_active.fasta
  head rnacentral_active.fasta > example.txt
  gzip rnacentral_active.fasta
  cp template.txt readme.txt
  """
}

process inactive {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path(query)

  output:
  path('rnacentral_inactive.fasta.gz')

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" | json2fasta.py - - | gzip > rnacentral_inactive.fasta.gz
  """
}

process species_specific {
  input:
  path(query)

  output:
  path('rnacentral_species_specific_ids.fasta')

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -q -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral_species_specific_ids.fasta
  """
}

process create_ssi {
  publishDir "${params.export.ftp.publish}/sequences/.internal/", mode: 'copy'
  memory '8GB'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.ssi')

  script:
  """
  esl-sfetch --index rnacentral_species_specific_ids.fasta
  """
}

process compress_species_fasta {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.gz')

  script:
  """
  gzip < rnacentral_species_specific_ids.fasta > rnacentral_species_specific_ids.fasta.gz
  """
}


process database_specific {
  tag { db }
  publishDir "${params.export.ftp.publish}/sequences/by-database", mode: 'copy'

  input:
  tuple val(db), path('raw.json')

  output:
  file('*.fasta.gz')

  script:
  """
  export PYTHONIOENCODING=utf8
  json2fasta.py raw.json ${db.toLowerCase().replaceAll(' ', '_').replace('/', '_')}.fasta
  gzip ${db.toLowerCase().replaceAll(' ', '_').replace('/', '_')}.fasta
  """
}

process by_species_readme {
  publishDir "${params.export.ftp.publish}/sequences/by-species/", mode: 'copy'

  input:
  path('template.txt')

  output:
  path('readme.txt')

  script:
  """
  cp template.txt readme.txt
  """
}

process find_species {
  input:
  path(query)

  output:
  path('species.csv')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > species.csv
  """
}

process split_by_species {
  cpus 4
  memory '16 GB'
  time '12h'
  publishDir "${params.export.ftp.publish}/sequences/by-species", mode: 'copy'

  input:
  tuple path(species_csv), path(fasta)

  output:
  path('**/*.fasta.gz')

  script:
  """
  set -euo pipefail

  split_by_species.py ${species_csv} ${fasta} .

  # Compress every per-species file in parallel with stock gzip.
  find . -mindepth 2 -name '*.fasta' -print0 | xargs -0 -r -P ${task.cpus} gzip
  """
}

workflow fasta_export {
  take:
    active_json
    database_dumps
  main:
  channel.fromPath('files/ftp-export/sequences/readme.txt') | set { readme }
  active(active_json, readme)

  channel.fromPath('files/ftp-export/sequences/inactive.sql') | inactive
  channel.fromPath('files/ftp-export/sequences/species-specific.sql') \
  | species_specific \
  | (compress_species_fasta & create_ssi)

  database_dumps | database_specific

  channel.fromPath('files/ftp-export/sequences/species.sql') \
  | find_species \
  | combine(species_specific.out) \
  | split_by_species

  channel.fromPath('files/ftp-export/sequences/by-species-readme.txt') | by_species_readme

  emit:
    active.out.active.mix(
      inactive.out,
      compress_species_fasta.out,
      create_ssi.out,
      database_specific.out,
      split_by_species.out,
      by_species_readme.out
    )
}
