process dump {
  publishDir "${params.export.ftp.publish}/editing/", mode: 'copy'

  input:
    path query

  output:
    path('rediportal-mapping.csv')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query "\$PGDATABASE" > rediportal-mapping.csv
  """

}


workflow rediportal {
  main:

  channel.fromPath('files/ftp-export/rediportal/editing-locations.sql') | set { query }

  query | dump

  emit: dump.out
}
