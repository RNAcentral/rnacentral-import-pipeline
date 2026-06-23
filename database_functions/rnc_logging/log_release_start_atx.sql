CREATE OR REPLACE FUNCTION rnc_logging.log_release_start_atx(p_dbid bigint, p_this_release bigint, p_prev_release bigint DEFAULT NULL::bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    l_dbid RNACEN.rnc_release.dbid%TYPE := p_dbid;
    l_this_release RNACEN.rnc_release.id%TYPE := p_this_release;
    l_prev_release RNACEN.rnc_release.id%TYPE := coalesce(coalesce(p_prev_release, release.get_previous_release(l_dbid, l_this_release)), 0);
    l_sql_stmt varchar(4000);
    l_start_time RNACEN.release_stats.start_time%TYPE := clock_timestamp();
BEGIN

with new_stats as (
    SELECT
          l_dbid dbid,
          l_this_release this_release,
          l_start_time start_time,
           COUNT(taxid) ff_taxid_count,
           COUNT(*) ff_loaded_rows           
      FROM rnacen.load_rnacentral
    )
insert INTO rnacen.release_stats as s (
  dbid,
  this_release,
  start_time,
  ff_taxid_nulls,
  ff_loaded_rows      
)
 select
        q.dbid,
        q.this_release,
        q.start_time,
        (q.ff_loaded_rows - q.ff_taxid_count) ff_taxid_nulls,
        q.ff_loaded_rows
   from new_stats as q
   on conflict (this_release) do update
    SET
      dbid = excluded.dbid,
      start_time = excluded.start_time,
      ff_taxid_nulls = excluded.ff_taxid_nulls,
      ff_loaded_rows = excluded.ff_loaded_rows;

  END;

$function$

