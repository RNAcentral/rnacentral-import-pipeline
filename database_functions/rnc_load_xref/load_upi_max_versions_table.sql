CREATE OR REPLACE FUNCTION rnc_load_xref.load_upi_max_versions_table(p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'TRUNCATE TABLE load_upi_max_versions';

    INSERT INTO load_upi_max_versions
      SELECT DISTINCT
        L.IN_AC,
        L.IN_DBID,
        MAX(coalesce(PREVIOUS_XREF.VERSION_I, 0)) MAX_VERSION_I, -- VERSION_I set to 0 at first
        PREVIOUS_XREF.UPI
      FROM rnacen.xref previous_xref RIGHT OUTER JOIN load_retro_tmp l ON (PREVIOUS_XREF.DBID = L.IN_DBID AND PREVIOUS_XREF.AC = L.IN_AC)
     WHERE L.COMPARABLE_PROT_UPI IS NOT NULL AND L.IN_DBID = p_in_dbid   -- outer join, left side can be NULL
     -- outer join, left side can be NULL
       AND NOT EXISTS (
            SELECT 1
            FROM RNACEN.XREF X
            WHERE X.AC      = L.IN_AC
            AND X.DBID    = L.IN_DBID
            AND X.DELETED = 'N'
          )
      GROUP BY
        L.IN_AC,
        L.IN_DBID,
        PREVIOUS_XREF.UPI;

    --COMMIT;

    execute 'analyze load_upi_max_versions';

  END;

$function$

