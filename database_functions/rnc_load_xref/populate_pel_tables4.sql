CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables4(p_in_dbid bigint, v_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
/*
    INSERT INTO XREF_PEL_DELETED(
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        TIMESTAMP,
        USERSTAMP,
        AC,
        VERSION,
        LAST,
        DELETED,
        TAXID
      )
    SELECT
      x.DBID,
      x.CREATED,
      x.UPI,
      x.VERSION_I,
      x.TIMESTAMP,
      x.USERSTAMP,
      x.AC,
      x.VERSION,
      CASE X.DELETED
        WHEN 'N'
        THEN
          coalesce(v_previous_release, X.LAST) -- retire updated entries
        ELSE                               -- retired xref that weren't updated
          X.LAST
      END as "last",
      'Y' DELETED,
      x.taxid
    FROM RNACEN.XREF X
    WHERE X.dbid = p_in_dbid
    AND (
        x.deleted = 'Y'
      OR NOT EXISTS (
          SELECT 1
          FROM LOAD_RETRO_TMP L
          WHERE L.COMPARABLE_PROT_UPI IS NOT NULL
          AND X.AC                 = L.IN_AC
          AND X.DBID               = L.IN_DBID
        )
      );
*/

    INSERT INTO XREF_PEL_DELETED(
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        TIMESTAMP,
        USERSTAMP,
        AC,
        VERSION,
        LAST,
        DELETED,
        TAXID
      )
    SELECT
      x.DBID,
      x.CREATED,
      x.UPI,
      x.VERSION_I,
      x.TIMESTAMP,
      x.USERSTAMP,
      x.AC,
      x.VERSION,
      CASE X.DELETED
        WHEN 'N'
        THEN
          coalesce(v_previous_release, X.LAST) -- retire updated entries
        ELSE                               -- retired xref that weren't updated
          X.LAST
      END as "last",
      'Y' DELETED,
      x.taxid
    FROM RNACEN.XREF X
    WHERE X.dbid = p_in_dbid
    AND x.deleted = 'Y'
    UNION ALL
    SELECT
      x.DBID,
      x.CREATED,
      x.UPI,
      x.VERSION_I,
      x.TIMESTAMP,
      x.USERSTAMP,
      x.AC,
      x.VERSION,
      CASE X.DELETED
        WHEN 'N'
        THEN
          coalesce(v_previous_release, X.LAST) -- retire updated entries
        ELSE                               -- retired xref that weren't updated
          X.LAST
      END as "last",
      'Y' DELETED,
      x.taxid
    FROM RNACEN.XREF X
    WHERE X.dbid = p_in_dbid
    AND x.deleted = 'N'
    AND NOT EXISTS (
          SELECT 1
          FROM LOAD_RETRO_TMP L
          WHERE L.COMPARABLE_PROT_UPI IS NOT NULL
          AND X.AC                 = L.IN_AC
          AND X.DBID               = L.IN_DBID
        )
    ;

    --COMMIT;

END;

$function$

