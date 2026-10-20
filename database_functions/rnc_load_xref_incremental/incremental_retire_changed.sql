-- The one-arg signature this replaces must go, or CREATE OR REPLACE leaves it
-- behind as an overload and the loader keeps calling the unpruned version.
DROP FUNCTION IF EXISTS rnc_load_xref_incremental.incremental_retire_changed(bigint);

CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_retire_changed(p_in_dbid bigint, p_previous_release bigint)
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
    WHERE
      -- Constant on the partition key. Joining xref on l.in_dbid alone compares
      -- the key to another table's column, which prunes nothing, so the plan
      -- reads all 59 partitions. load_retro_tmp holds only this dbid, so this
      -- changes no rows.
      u.dbid = p_in_dbid
      AND u.deleted = 'N'
      AND u.ac   = l.in_ac
      AND u.dbid = l.in_dbid
      AND u.urs  = l.comparable_prot_upi
      AND l.comparable_prot_upi IS NOT NULL
      AND u.last < l.in_load_release
      AND NOT (u.version = l.in_version OR (u.version IS NULL AND l.in_version IS NULL));
END;
$function$
