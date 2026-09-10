CREATE OR REPLACE FUNCTION rnc_load_xref.do_checks(p_in_db_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    row record;
    dup_sql text;
BEGIN
  -- assign new id if pk is null
  -- perform rnc_update.verify_xref_id_not_null();

  -- id is a single sequence (xref_pk_seq) and each partition already carries
  -- its own unique index on id, so a duplicate can only be this dbid's ids
  -- colliding with some other, untouched partition -- no need to re-aggregate
  -- the whole table.
  dup_sql := format(
    'select x.id, count(*) as cnt
       from xref x
      where x.id in (
        select id from xref_p%1$s_deleted
        union all
        select id from xref_p%1$s_not_deleted
      )
      group by x.id
     having count(*) > 1',
    p_in_db_id
  );

  -- create MV with no data if not exists
  execute 'drop materialized view if exists xref_pk_not_unique';
  execute 'create materialized view xref_pk_not_unique as (' || dup_sql || ') with no data';

  -- refresh MV (populate with current data)
  execute 'refresh materialized view xref_pk_not_unique';

  -- raise exception if pk is violated across inherited tables
  for row in select id, cnt from xref_pk_not_unique
  loop
      RAISE 'ID: % has % duplicates in tables inheriting xref', row.id, row.cnt USING ERRCODE = 'unique_violation';
  end loop;

END;

$function$
