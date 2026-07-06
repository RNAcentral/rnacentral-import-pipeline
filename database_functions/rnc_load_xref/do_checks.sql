CREATE OR REPLACE FUNCTION rnc_load_xref.do_checks(p_in_db_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    row record;
BEGIN
  -- assign new id if pk is null
  -- perform rnc_update.verify_xref_id_not_null();

  -- create MV with no data if not exists
  execute 'drop materialized view if exists xref_pk_not_unique';
  execute 'create materialized view xref_pk_not_unique as (select id, count (*) as cnt from xref group by id having count (*) > 1) with no data';

  -- refresh MV (populate with current data)
  execute 'refresh materialized view xref_pk_not_unique';

  -- raise exception if pk is violated across inherited tables
  for row in select id, cnt from xref_pk_not_unique
  loop
      RAISE 'ID: % has % duplicates in tables inheriting xref', row.id, row.cnt USING ERRCODE = 'unique_violation';
  end loop;

END;

$function$

