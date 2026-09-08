CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables3(p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    INSERT INTO RNACEN.XREF_PEL_NOT_DELETED(
        AC,
        DBID,
        VERSION,
        VERSION_I,
        URS,
        CREATED,
        LAST,
        DELETED,
        TAXID,
        TIMESTAMP,
        USERSTAMP
      )
    SELECT
      T.IN_AC,
      T.IN_DBID,
      T.IN_VERSION,
      CASE
        WHEN T.MAX_PREVIOUS_XREF_VERSION_I = 0 -- new xref
          THEN 1
        WHEN LUMV.UPI = T.COMPARABLE_PROT_UPI  -- same UPI, keep old VERSION_I
          THEN T.MAX_PREVIOUS_XREF_VERSION_I
        ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1 -- updated URS, update VERSION_I
      END,
      T.COMPARABLE_PROT_UPI UPI,
      T.IN_LOAD_RELEASE CREATED,
      T.IN_LOAD_RELEASE as "last",
      'N' DELETED,
      T.IN_TAXID,
      clock_timestamp() as "timestamp",
      USER USERSTAMP
    FROM load_upi_max_versions lumv RIGHT OUTER JOIN (
        SELECT
          L.IN_AC,
          L.IN_DBID,
          L.IN_VERSION,
          LMV.MAX_VERSION_I MAX_PREVIOUS_XREF_VERSION_I,
          L.COMPARABLE_PROT_UPI,
          L.IN_LOAD_RELEASE,
          L.IN_TAXID
        FROM
          LOAD_RETRO_TMP L,
          LOAD_MAX_VERSIONS LMV
        WHERE
          L.COMPARABLE_PROT_UPI IS NOT NULL
        AND L.IN_DBID = p_in_dbid
        AND LMV.AC           = L.IN_AC
        AND LMV.DBID         = L.IN_DBID
      ) t ON (LUMV.AC = T.IN_AC AND LUMV.MAX_VERSION_I = T.MAX_PREVIOUS_XREF_VERSION_I AND LUMV.DBID = T.IN_DBID);

    -- Gap A: new sequences for accessions that ALREADY have active rows.
    -- The block above only covers accessions with NO current active row (that is
    -- baked into load_upi_max_versions_table's NOT EXISTS deleted='N' filter). Now
    -- that populate_pel_tables1/2 join on (ac, urs), an incoming sequence whose
    -- (ac, comparable_prot_upi) is not currently active has no other insert path,
    -- so handle it here. Disjoint from the block above (EXISTS active row) and from
    -- tables1/2 (which require a MATCHING active urs). On a no-change reload the
    -- inner NOT EXISTS is always false, so this inserts nothing.
    -- version_i choice: next generation for the accession, i.e. max(version_i)+1,
    -- matching this function's existing "updated URS -> max+1" rule. REVIEW if
    -- multi-sequence accessions need different version_i semantics.
    INSERT INTO RNACEN.XREF_PEL_NOT_DELETED(
        AC, DBID, VERSION, VERSION_I, URS, CREATED, LAST, DELETED, TAXID, TIMESTAMP, USERSTAMP
      )
    SELECT
      L.IN_AC,
      L.IN_DBID,
      L.IN_VERSION,
      COALESCE(mv.max_version_i, 0) + 1,
      L.COMPARABLE_PROT_UPI,
      L.IN_LOAD_RELEASE,
      L.IN_LOAD_RELEASE,
      'N',
      L.IN_TAXID,
      clock_timestamp(),
      USER
    FROM LOAD_RETRO_TMP L
    LEFT JOIN (
        SELECT AC, MAX(VERSION_I) max_version_i
        FROM RNACEN.XREF
        WHERE DBID = p_in_dbid
        GROUP BY AC
      ) mv ON mv.AC = L.IN_AC
    WHERE L.IN_DBID = p_in_dbid
    AND L.COMPARABLE_PROT_UPI IS NOT NULL
    AND EXISTS (
          SELECT 1 FROM RNACEN.XREF X
          WHERE X.AC = L.IN_AC AND X.DBID = p_in_dbid AND X.DELETED = 'N'
        )
    AND NOT EXISTS (
          SELECT 1 FROM RNACEN.XREF X
          WHERE X.AC = L.IN_AC AND X.DBID = p_in_dbid AND X.DELETED = 'N'
          AND X.URS = L.COMPARABLE_PROT_UPI
        );

    --COMMIT;

  END;

$function$
