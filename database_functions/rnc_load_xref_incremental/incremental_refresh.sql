CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_refresh(p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- In-place analogue of populate_pel_tables1 (deleted='N' branch, case 1): refresh
    -- an unchanged active row's last/taxid (the common case on a reload). The
    -- last < in_load_release guard skips rows inserted earlier this release.
    UPDATE rnacen.xref u
    SET last  = l.in_load_release,
        taxid = COALESCE(l.in_taxid, u.taxid)
    FROM load_retro_tmp l
    WHERE u.deleted = 'N'
      AND u.ac   = l.in_ac
      AND u.dbid = l.in_dbid
      AND u.urs  = l.comparable_prot_upi
      AND l.comparable_prot_upi IS NOT NULL
      AND u.last < l.in_load_release
      AND (u.version = l.in_version OR (u.version IS NULL AND l.in_version IS NULL));
END;
$function$
