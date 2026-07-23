CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental_new_accessions(p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- In-place analogue of populate_pel_tables3 (cases D + Gap A), retargeted to
    -- rnacen.xref. Uses load_upi_max_versions / load_max_versions, populated from the
    -- original xref before any mutation.

    -- Case D: accessions with no current active row (brand new or previously deleted).
    INSERT INTO rnacen.xref (
        ac, dbid, version, version_i, upi,
        created, last, deleted, taxid, timestamp, userstamp
    )
    SELECT
        t.in_ac,
        t.in_dbid,
        t.in_version,
        CASE
            WHEN t.max_previous_xref_version_i = 0 THEN 1   -- new xref
            WHEN lumv.upi = t.comparable_prot_upi           -- same UPI, keep version_i
                THEN t.max_previous_xref_version_i
            ELSE t.max_previous_xref_version_i + 1          -- updated UPI, bump version_i
        END,
        t.comparable_prot_upi,
        t.in_load_release,
        t.in_load_release,
        'N',
        t.in_taxid,
        clock_timestamp(),
        current_user
    FROM load_upi_max_versions lumv
    RIGHT OUTER JOIN (
        SELECT
            l.in_ac,
            l.in_dbid,
            l.in_version,
            lmv.max_version_i AS max_previous_xref_version_i,
            l.comparable_prot_upi,
            l.in_load_release,
            l.in_taxid
        FROM load_retro_tmp l, load_max_versions lmv
        WHERE l.comparable_prot_upi IS NOT NULL
          AND l.in_dbid = p_in_dbid
          AND lmv.ac    = l.in_ac
          AND lmv.dbid  = l.in_dbid
    ) t
      ON  lumv.ac            = t.in_ac
      AND lumv.max_version_i = t.max_previous_xref_version_i
      AND lumv.dbid          = t.in_dbid;

    -- Gap A: a new sequence variant for an accession that already has active rows but
    -- none with this upi. version_i = max(version_i) + 1, matching populate_pel_tables3.
    INSERT INTO rnacen.xref (
        ac, dbid, version, version_i, upi,
        created, last, deleted, taxid, timestamp, userstamp
    )
    SELECT
        l.in_ac,
        l.in_dbid,
        l.in_version,
        COALESCE(mv.max_version_i, 0) + 1,
        l.comparable_prot_upi,
        l.in_load_release,
        l.in_load_release,
        'N',
        l.in_taxid,
        clock_timestamp(),
        current_user
    FROM load_retro_tmp l
    LEFT JOIN (
        SELECT ac, MAX(version_i) AS max_version_i
        FROM rnacen.xref
        WHERE dbid = p_in_dbid
        GROUP BY ac
    ) mv ON mv.ac = l.in_ac
    WHERE l.in_dbid = p_in_dbid
      AND l.comparable_prot_upi IS NOT NULL
      AND EXISTS (
            SELECT 1 FROM rnacen.xref x
            WHERE x.ac = l.in_ac AND x.dbid = p_in_dbid AND x.deleted = 'N'
          )
      AND NOT EXISTS (
            SELECT 1 FROM rnacen.xref x
            WHERE x.ac = l.in_ac AND x.dbid = p_in_dbid AND x.deleted = 'N'
              AND x.upi = l.comparable_prot_upi
          );
END;
$function$
