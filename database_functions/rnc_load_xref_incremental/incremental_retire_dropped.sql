CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_retire_dropped(p_in_dbid bigint, p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- In-place analogue of populate_pel_tables4's second arm (case 5): mark inactive an
    -- active (ac, urs) that vanished from the load -- the step the old Oracle path lacked.
    -- Its first arm (carry deleted history forward) has no analogue: leaving those rows
    -- untouched is the whole win. NOT EXISTS spares this release's inserts (all in-load).
    UPDATE rnacen.xref u
    SET deleted = 'Y',
        last    = COALESCE(p_previous_release, u.last)
    WHERE u.dbid = p_in_dbid
      AND u.deleted = 'N'
      AND NOT EXISTS (
            SELECT 1
            FROM load_retro_tmp l
            WHERE l.comparable_prot_upi IS NOT NULL
              AND l.in_ac   = u.ac
              AND l.in_dbid = u.dbid
              AND l.comparable_prot_upi = u.urs
          );
END;
$function$
