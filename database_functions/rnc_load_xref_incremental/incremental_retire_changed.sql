CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_retire_changed(p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- In-place analogue of populate_pel_tables1 (deleted='Y' branch, case 2): retire an
    -- active row whose version changed (the row moves to the _deleted partition). The
    -- replacement was inserted by incremental_new_versions; the last-guard excludes
    -- this release's inserts.
    UPDATE rnacen.xref u
    SET deleted = 'Y',
        last    = COALESCE(p_previous_release, u.last)
    FROM load_retro_tmp l
    WHERE u.deleted = 'N'
      AND u.ac   = l.in_ac
      AND u.dbid = l.in_dbid
      AND u.urs  = l.comparable_prot_upi
      AND l.comparable_prot_upi IS NOT NULL
      AND u.last < l.in_load_release
      AND NOT (u.version = l.in_version OR (u.version IS NULL AND l.in_version IS NULL));
END;
$function$
