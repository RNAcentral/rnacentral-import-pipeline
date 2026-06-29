process build_id_mapping_chunk {
  tag { chunk }
  maxForks 4

  input:
  tuple val(chunk), path(query)

  output:
  path("chunk_${chunk}.tsv")

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -v chunk=${chunk} -f "$query" "\$PGDATABASE" > raw_${chunk}.json
  rnac ftp-export id-mapping raw_${chunk}.json - | sort -T . -u > chunk_${chunk}.tsv
  """
}

process merge_id_mapping {
  publishDir "${params.export.ftp.publish}/id_mapping/", mode: 'copy'

  input:
  path('chunk*.tsv')
  path('template.txt')

  output:
  path('id_mapping.tsv.gz'), emit: mapping
  path("example.txt"), emit: 'example'
  path("readme.txt"), emit: 'readme'

  script:
  """
  set -euo pipefail

  sort -T . -m -u chunk*.tsv > id_mapping.tsv
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

  script:
  """
  set -euo pipefail

  zcat id_mapping.tsv.gz | awk '{ print >> (tolower(\$2) ".tsv") }'
  """
}

workflow id_mapping {
  Channel.fromPath('files/ftp-export/id-mapping/id_mapping.sql') | set { id_query }
  Channel.fromPath('files/ftp-export/id-mapping/readme.txt') | set { readme_template }

  Channel.of('0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f') \
  | combine(id_query) \
  | build_id_mapping_chunk \
  | collect \
  | combine(readme_template) \
  | merge_id_mapping

  merge_id_mapping.out.mapping | database_mapping
}
