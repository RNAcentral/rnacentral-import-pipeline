process dump_data {
  container ''

  output:
  path("prod-dump")

  script:
  prod = params.dbs.production
  """
  /nfs/dbtools/pgsql10/usr/pgsql-10/bin/pg_dump \
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
  	--exclude-table=xref_pk_not_unique \
  	-h ${prod.host} -U ${prod.user} -d ${prod.db_name} -j 4 -n ${prod.db_name} -F d -O -f "prod-dump"
  """
}

process restore_public {
  container ''

  input:
  path(dump_file)

  script:
  public_db = params.dbs.public_admin
  """
  PUBLIC='postgres://${public_db.user}:${public_db.password}@${public_db.host}:${public_db.port}/${public_db.db_name}'

  export PGPASSWORD='${public.password}'
  psql -c "set maintenance_work_mem='1GB'" $PUBLIC
  psql -c 'drop schema if exists rnacen cascade' $PUBLIC 

  /nfs/dbtools/pgsql10/usr/pgsql-10/bin/pg_restore -x -h ${public_db.host} -U ${public_db.user} -d ${public_db.db_name} -j 2 $dump_file

  psql -c 'revoke usage on schema rnacen from public' $PUBLIC
  psql -c 'grant usage on schema rnacen to reader' $PUBLIC
  psql -c 'grant SELECT on ALL tables in schema rnacen to reader' $PUBLIC
  psql -c 'analyze' $PUBLIC
  """
}

workflow update_public {
  dump_data | restore_public
}

workflow {
  update_public()
}
