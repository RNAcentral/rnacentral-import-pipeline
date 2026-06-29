CREATE OR REPLACE FUNCTION rnc_load_xref.load_max_versions_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    EXECUTE 'TRUNCATE TABLE load_max_versions';

    INSERT INTO load_max_versions
    SELECT DISTINCT
      l.AC,
      L.DBID,
      MAX(L.MAX_VERSION_I) MAX_VERSION_I
    FROM load_upi_max_versions l
    GROUP BY
      L.AC,
      L.DBID;

    --COMMIT;

    execute 'analyze load_max_versions';

END;

$function$

