process build_id_mapping {
  publishDir "${params.export.ftp.publish}/id_mapping/", mode: 'copy'

  input:
  path(query)
  path('template.txt')

  output:
  path('id_mapping.tsv.gz'), emit: mapping
  path("example.txt"), emit: 'example'
  path("readme.txt"), emit: 'readme'

  """
  set -euo pipefail

  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw_id_mapping.tsv
  rnac ftp-export id-mapping raw_id_mapping.tsv - | sort -u > id_mapping.tsv
  head id_mapping.tsv > example.txt
  gzip id_mapping.tsv
  cat template.txt > readme.txt
  """
}

process database_mapping {
  publishDir "${params.export.ftp.publish}/id_mapping/database_mappings/", mode: 'copy'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('id_mapping.tsv.gz')

  output:
  path('*.tsv')

  """
  set -euo pipefail

  zcat id_mapping.tsv.gz | awk '{ print >> (tolower(\$2) ".tsv") }'
  """
}

workflow id_mapping {
  Channel.fromPath('files/ftp-export/id-mapping/id_mapping.sql') | set { id_query }
  Channel.fromPath('files/ftp-export/id-mapping/readme.txt') | set { readme_template }

  build_id_mapping(id_query, readme_template)

  build_id_mapping.out.mapping | database_mapping
}
