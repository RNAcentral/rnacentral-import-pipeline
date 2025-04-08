#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process dump_data {
  container ''
  time '2d'

  output:
  path("prod-dump")

  script:
  prod = params.dbs.production
  """
  export PGPASSWORD='${prod.password}'

  /usr/pgsql-11/bin/pg_dump \
    -h ${prod.host} -U ${prod.user} -d ${prod.db_name} -j 4 --schema ${prod.schema} -F d -O -f 'prod-dump' \
    --exclude-table=auth_group \
    --exclude-table=auth_group_permissions \
    --exclude-table=auth_user \
    --exclude-table=auth_user_groups \
    --exclude-table=auth_user_user_permissions \
    --exclude-table=corsheaders_corsmodel \
    --exclude-table=django_admin_log \
    --exclude-table=django_content_type \
    --exclude-table=django_migrations \
    --exclude-table=django_session \
    --exclude-table=django_site \
    --exclude-table=load_ensembl_analysis_status \
    --exclude-table=load_flybase \
    --exclude-table=load_genome_mapping \
    --exclude-table=load_max_versions \
    --exclude-table=load_md5_collisions \
    --exclude-table=load_md5_new_sequences \
    --exclude-table=load_md5_stats \
    --exclude-table=load_ontology_terms \
    --exclude-table=load_overlaps \
    --exclude-table=load_retro_tmp \
    --exclude-table=load_rnacentral \
    --exclude-table=load_rnacentral_all \
    --exclude-table=load_rnc_accessions \
    --exclude-table=load_rnc_references \
    --exclude-table=load_rnc_related_sequences \
    --exclude-table=load_upi_max_versions \
    --exclude-table=xref_pk_not_unique
  """
}

process populate_public {
  container ''
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
  /usr/pgsql-10/bin/pg_restore -x -h ${public_db.host} -U ${public_db.user} -d ${public_db.db_name} -j 2 $dump_file

  psql -c 'revoke usage on schema rnacen from public' \$PUBLIC
  psql -c 'grant usage on schema rnacen to reader' \$PUBLIC
  psql -c 'grant SELECT on ALL tables in schema rnacen to reader' \$PUBLIC
  psql -c 'analyze' \$PUBLIC
  """
}

workflow update_public {
  dump_data | populate_public
}

workflow {
  update_public()
}
