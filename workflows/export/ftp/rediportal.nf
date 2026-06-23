process dump {
  publishDir "${params.export.ftp.publish}/editing/", mode: 'copy'

  input:
    path query


  script:
  """
  psql -f $query "\$PGDATABASE" > rediportal-mapping.csv
  """

}


workflow rediportal {

  Channel.fromPath('files/ftp-export/rediportal/editing-locations.sql') | set { query }

  query | dump

}
