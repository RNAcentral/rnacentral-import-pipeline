CREATE OR REPLACE FUNCTION rnc_logging.log_release_start(p_dbid bigint, p_this_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
        v_query     text;

BEGIN
        -- v_query := 'SELECT true FROM RNC_LOGGING.log_release_start_atx ( ' || quote_nullable(p_dbid) || ',' || quote_nullable(p_this_release) || ' )';
        -- PERFORM * FROM pg_background_result(pg_background_launch(v_query)) AS p (ret boolean);
        PERFORM RNC_LOGGING.log_release_start_atx ( p_dbid , p_this_release );
END;
$function$

