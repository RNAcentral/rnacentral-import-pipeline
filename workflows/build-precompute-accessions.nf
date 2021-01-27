process find_partitions {
  output:
  path('paritions.csv')

  """
  psql \
    -c 'COPY (SELECT id, descr from rnc_database) TO STDOUT (FORMAT CSV)' \
  "$PGDATABASE" | \
  awk -F, '{
    print "xref_p"\$1"_not_deleted", \$2
    print "xref_p"\$1"_deleted", \$2
  }' > partitions.csv
  """
}

process setup_accession_table {
  executor 'local'

  input:
  path(sql)

  output:
  val('done')

  """
  psql -v ON_ERROR_STOP=1 -f $sql
  """
}

process build_accession_table {
  maxForks 1

  input:
  tuple val(partition), val(db), path(sql)

  output:
  val('done')

  """
  psql \
  -v ON_ERROR_STOP=1 \
  -v "db_name=$db" \
  -v "partition=$partition" \
  -f $sql \
  "$PGDATABASE"
  """
}

process finalize_accession_table {
  input:
  path(sql)

  output:
  val('done')

  """
  psql -v ON_ERROR_STOP=1 -f $sql
  """
}


workflow build_precompute_accessions {
  take: ready
  emit: built
  main:
    Channel.fromPath('files/precompute/get-accessions/schema.sql') | set { schema }
    Channel.fromPath('files/precompute/get-accessions/insert-chunk.sql') | set { chunk_sql }
    Channel.fromPath('files/precompute/get-accessions/post-insert.sql') | set { post_sql }

    find_partitions | splitCsv | set { partitions }

    ready \
    | combine(schema) \
    | map { _flag, sql -> sql } \
    | setup_accession_table \
    | combine(partitions) \
    | map { _flag, p, descr -> [p, descr] } \
    | combine(chunk_sql) \
    | build_accession_table \
    | collect \
    | combine(post_sql) \
    | finalize_accession_table \
    | set { built }
}
