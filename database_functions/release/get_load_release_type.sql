CREATE OR REPLACE FUNCTION release.get_load_release_type(in_dbid bigint)
 RETURNS character
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE
    v_has_prior boolean;
    v_is_delta  boolean := false;
BEGIN
    -- First load of a database is FULL (bootstraps the xref partition). After that a
    -- delta-parsed database (one that has an import manifest) uses a DELTA load, so
    -- absent rows are never retired; any other database uses an INCREMENTAL load.
    -- prepare_releases('F') can still force FULL for all.
    SELECT EXISTS (
        SELECT 1 FROM RNACEN.rnc_release WHERE dbid = in_dbid
    ) INTO v_has_prior;

    IF NOT v_has_prior THEN
        RETURN 'F';
    END IF;

    IF to_regclass('rnacen.pipeline_tracking_import') IS NOT NULL THEN
        SELECT EXISTS (
            SELECT 1
            FROM rnacen.pipeline_tracking_import m
            JOIN rnacen.rnc_database d ON d.descr = m.database
            WHERE d.id = in_dbid
        ) INTO v_is_delta;
    END IF;

    IF v_is_delta THEN
        RETURN 'D';
    ELSE
        RETURN 'I';
    END IF;
END;
$function$
