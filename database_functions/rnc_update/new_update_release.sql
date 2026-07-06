CREATE OR REPLACE FUNCTION rnc_update.new_update_release(p_in_dbid bigint, p_in_release_id bigint, p_release_type character DEFAULT 'f'::bpchar)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
declare

    c_load cursor for
       select
              dbid,
              id,
              release_type,
              release_date,
              force_load
         from rnacen.rnc_release
        where status = 'L' 
          and dbid = p_in_dbid 
          and id = p_in_release_id
        order by id;

begin

    raise notice 'Launching update of dbid %, release id %', p_in_dbid, p_in_release_id;

    perform rnc_update.move_staging_data(p_in_dbid => p_in_dbid);
    perform rnc_update.load_release(
        p_in_dbid         => p_in_dbid,
        p_in_load_release => p_in_release_id
    );

  end;
$function$

