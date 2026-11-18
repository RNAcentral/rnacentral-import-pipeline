#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process dump_data {

  time '2d'

  output:
  path("prod-dump")

  script:
  prod = params.dbs.production
  """
  export PGPASSWORD='${prod.password}'

  pg_dump \
    -h ${prod.host} -U ${prod.user} -d ${prod.db_name} -j 4 --schema ${prod.schema} -F d -O -f 'prod-dump' \
    --exclude-table=auth_group \
    --exclude-table=auth_group_permissions \
    --exclude-table=auth_permission \
    --exclude-table=auth_user \
    --exclude-table=auth_user_groups \
    --exclude-table=auth_user_user_permissions \
    --exclude-table=corsheaders_corsmodel \
    --exclude-table=django_admin_log \
    --exclude-table=django_content_type \
    --exclude-table=django_migrations \
    --exclude-table=django_session \
    --exclude-table=django_site \
    --exclude-table='load_*' \
    --exclude-table='pipeline_tracking_*' \
    --exclude-table=xref_pk_not_unique
  """
}

process populate_public {

  time '4d'

  input:
  path(dump_file)

  script:
  public_db = params.dbs.pub
  """
  psql -c "set maintenance_work_mem='1GB'" \$PUBLIC
  psql -c 'drop schema if exists rnacen cascade' \$PUBLIC
  psql -c "ALTER ROLE rnacen SET statement_timeout = '2d';" \$PUBLIC

  export PGPASSWORD='${public_db.password}'
  pg_restore -x -h ${public_db.host} -U ${public_db.user} -d ${public_db.db_name} -j 2 $dump_file

  psql -c 'revoke usage on schema rnacen from public' \$PUBLIC
  psql -c 'grant usage on schema rnacen to reader' \$PUBLIC
  psql -c 'grant SELECT on ALL tables in schema rnacen to reader' \$PUBLIC
  psql -c 'analyze' \$PUBLIC
  """
}

process publish_dump {
  publishDir "${params.export.ftp.publish}/database_files/", mode: 'copy'
  time '1d'

  input:
  path(dump)

  output:
  path('database_files.tar')

  script:
  """
  tar -cf database_files.tar -C $dump .
  """
}

workflow update_public {
  take: _ready
  main:
    dump_data | set { dump }
    dump | populate_public
    dump | publish_dump
  emit: publish_dump.out
}

workflow {
  update_public(channel.of('ready'))
}
