-- CREATE OR REPLACE cannot rename an input parameter, and this one goes from
-- p_previous_release to p_in_dbid, so the old function has to be dropped first.
DROP FUNCTION IF EXISTS rnc_load_xref_incremental.incremental_new_versions(bigint);

CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_new_versions(p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- In-place analogue of populate_pel_tables2 (case C): insert a fresh active
    -- generation for an accession whose sequence version changed (same urs); the old
    -- row is retired by incremental_retire_changed. version_i carries over unchanged
    -- (the join forces urs equality), exactly as in the full path.
    INSERT INTO rnacen.xref (
        ac, dbid, version, version_i, urs,
        created, last, deleted, taxid, timestamp, userstamp
    )
    SELECT
        l.in_ac,
        l.in_dbid,
        l.in_version,
        CASE WHEN x.urs != l.comparable_prot_upi THEN x.version_i + 1
             ELSE x.version_i END,
        l.comparable_prot_upi,
        l.in_load_release,               -- created
        l.in_load_release,               -- last
        'N',
        l.in_taxid,
        clock_timestamp(),
        current_user
    FROM load_retro_tmp l
    JOIN rnacen.xref x
      ON  x.ac   = l.in_ac
      AND x.dbid = l.in_dbid
      AND x.urs  = l.comparable_prot_upi
    WHERE
      -- Constant on the partition key. Joining xref on l.in_dbid alone compares
      -- the key to another table's column, which prunes nothing, so the plan
      -- reads all 59 partitions. load_retro_tmp holds only this dbid, so this
      -- changes no rows.
      x.dbid = p_in_dbid
      AND l.comparable_prot_upi IS NOT NULL
      AND x.deleted = 'N'
      AND x.last < l.in_load_release
      AND NOT (
            x.last < l.in_load_release
        AND (x.version = l.in_version OR (x.version IS NULL AND l.in_version IS NULL))
      );
END;
$function$
