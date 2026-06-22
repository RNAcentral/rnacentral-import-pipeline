CREATE OR REPLACE FUNCTION rnc_update.mark_as_done(p_in_dbid bigint, p_in_load_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    perform release.set_release_status( p_in_load_release, 'D' );
    perform database.set_current_release( p_in_dbid, p_in_load_release );
    --COMMIT;

  END;

$function$

