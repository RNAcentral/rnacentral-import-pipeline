CREATE OR REPLACE FUNCTION release.get_load_release_type(in_dbid bigint)
 RETURNS character
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE
    v_has_prior boolean;
    v_is_delta  boolean := false;
    v_descr     text;
BEGIN
    -- First load is FULL (bootstraps the partition). A manifest-tracked database
    -- gets DELTA; everything else stays FULL, since its parser resubmits the whole
    -- dataset anyway (see docs/incremental-xref-loading.md).
    SELECT EXISTS (
        SELECT 1 FROM RNACEN.rnc_release WHERE dbid = in_dbid
    ) INTO v_has_prior;

    IF NOT v_has_prior THEN
        RETURN 'F';
    END IF;

    -- HGNC is deliberately kept FULL, as it was imported before delta existed.
    -- FULL retires by absence, so the parser (rnac hgnc map) is pinned to a full
    -- parse to match. To run HGNC as a delta, lift both pins together.
    SELECT descr INTO v_descr FROM RNACEN.rnc_database WHERE id = in_dbid;
    IF v_descr = 'HGNC' THEN
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
        RETURN 'F';
    END IF;
END;
$function$
