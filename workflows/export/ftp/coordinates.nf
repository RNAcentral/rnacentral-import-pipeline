process readme {
  publishDir "${params.export.ftp.publish}/genome_coordinates/", mode: 'copy'

  input:
  path(raw)

  output:
  path('readme.txt')

  """
  cp $raw readme.txt
  """
}

process find_jobs {
  input:
  path(query)

  output:
  path('coordinates.txt')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > coordinates.txt
  """
}

process fetch {
  tag { "${assembly}-${species}" }
  maxForks 2

  input:
  tuple val(assembly), val(species), val(taxid), path(query)

  output:
  tuple val(assembly), val(species), path('result.json')

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "$PGDATABASE" > result.json
  """
}

process generate_bed {
  tag { "${assembly}-${species}" }
  publishDir "${params.export.ftp.publish}/genome_coordinates/bed/", mode: 'copy'

  input:
  tuple val(assembly), val(species), path(raw_data)

  output:
  tuple val(assembly), path("${species}.${assembly}.bed.gz")

  """
  set -euo pipefail

  rnac ftp-export coordinates as-bed $raw_data |\
  sort -k1,1 -k2,2n |\
  gzip > ${species}.${assembly}.bed.gz
  """
}

process generate_gff3 {
  tag { "${assembly}-${species}" }
  memory params.export.ftp.coordinates.gff3.memory
  publishDir "${params.export.ftp.publish}/genome_coordinates/gff3", mode: 'copy'

  input:
  tuple val(assembly), val(species), path(raw_data)

  output:
  path("${species}.${assembly}.gff3.gz")

  """
  set -euo pipefail

  rnac ftp-export coordinates as-gff3 $raw_data - |\
  sort -t"`printf '\\t'`" -k1,1 -k4,4n |\
  bgzip > "${species}.${assembly}".gff3.gz
  """
}

process index_gff3 {
  publishDir "${params.export.ftp.publish}/.genome-browser", mode: 'copy'

  input:
  path(gff)

  output:
  path(${gff.baseName}.gz.tbi)
  """
  tabix -p gff $gff
  """
}

workflow export_coordinates {
  Channel.fromPath('files/ftp-export/genome_coordinates/known-coordinates.sql') | set { known }
  Channel.fromPath('files/ftp-export/genome_coordinates/query.sql') | set { query }

  readme(Channel.fromPath('files/ftp-export/genome_coordinates/readme.mkd'))

  known \
  | find_jobs \
  | splitCsv \
  | combine(query) \
  | fetch \
  | filter { _a, _s, fn -> !fn.isEmpty() } \
  | set { coordinates }

  coordinates | generate_bed
  coordinates | generate_gff3 | index_gff3
}
