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
        UPI,
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
        ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1 -- updated UPI, update VERSION_I
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

    --COMMIT;

  END;

$function$

