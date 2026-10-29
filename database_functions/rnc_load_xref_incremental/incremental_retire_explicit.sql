CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_retire_explicit(p_in_dbid bigint, p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    v_descr text;
BEGIN
    -- DELTA-mode deletion: retire only the accessions the parser staged in
    -- load_deletions (scoped to this dbid by database name), never rows that are
    -- merely absent -- that is what keeps "absent rows are never deleted" true.
    UPDATE rnacen.xref u
    SET deleted = 'Y',
        last    = COALESCE(p_previous_release, u.last)
    WHERE u.dbid = p_in_dbid
      AND u.deleted = 'N'
      AND EXISTS (
            SELECT 1
            FROM load_deletions d
            JOIN rnacen.rnc_database db ON db.descr = d.database
            WHERE db.id = p_in_dbid
              AND d.accession = u.ac
          );

    -- WormBase accessions are '<ENA record id>:WORMBASE:<gene>', and their deletions
    -- arrive in ENA's list under the bare record id.
    SELECT descr INTO v_descr FROM rnacen.rnc_database WHERE id = p_in_dbid;
    IF v_descr = 'WORMBASE' THEN
        UPDATE rnacen.xref u
        SET deleted = 'Y',
            last    = COALESCE(p_previous_release, u.last)
        WHERE u.dbid = p_in_dbid
          AND u.deleted = 'N'
          AND EXISTS (
                SELECT 1
                FROM load_deletions d
                WHERE d.database = 'ENA'
                  AND d.accession = regexp_replace(u.ac, ':WORMBASE:[^:]+$', '')
              );
    END IF;
END;
$function$
