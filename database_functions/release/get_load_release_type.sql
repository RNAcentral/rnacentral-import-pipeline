CREATE OR REPLACE FUNCTION release.get_load_release_type(in_dbid bigint)
 RETURNS character
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE
    v_has_prior boolean;
BEGIN
    -- First load of a database is FULL (bootstraps the xref partition); every later
    -- load is INCREMENTAL. prepare_releases('F') can still force FULL for all.
    SELECT EXISTS (
        SELECT 1 FROM RNACEN.rnc_release WHERE dbid = in_dbid
    ) INTO v_has_prior;

    IF v_has_prior THEN
        RETURN 'I';
    ELSE
        RETURN 'F';
    END IF;
END;
$function$
