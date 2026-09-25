DROP FUNCTION IF EXISTS rnc_load_xref.do_checks(bigint);

CREATE OR REPLACE FUNCTION rnc_load_xref.do_checks(p_in_db_id bigint, p_in_load_release bigint)
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

  -- ids come from one sequence, so only rows inserted this release can collide;
  -- bounding by their id range gives one index range scan per partition (an IN
  -- list was planned as a probe per id per partition and ran for days). The CTE
  -- is MATERIALIZED so min/max cannot become a bottom-up walk of the id index.
  dup_sql := format(
    'with new_ids as materialized (
        select id from xref_p%1$s_deleted where created = %2$s
        union all
        select id from xref_p%1$s_not_deleted where created = %2$s
     )
     select x.id, count(*) as cnt
       from xref x
      where x.id between (select min(id) from new_ids)
                     and (select max(id) from new_ids)
      group by x.id
     having count(*) > 1',
    p_in_db_id, p_in_load_release
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
