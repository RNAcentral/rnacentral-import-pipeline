CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental2(p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN


    -- inserts new xref
    -- for each deleted = 'Y'
    -- in the above update
    INSERT
      /*+ PARALLEL (XREF 4) */
    INTO
      RNACEN.xref(
        ac,
        dbid,
        version,
        version_i,
        urs,

        created,
        last,
        deleted,
        taxid
      )
    SELECT
      /*+ PARALLEL (L 4) PARALLEL (X 4) */
      l.IN_AC,
      l.IN_DBID,
      l.IN_VERSION,
      CASE
        WHEN(

            X.UPI != L.COMPARABLE_PROT_UPI
          )
        OR (
            X."VERSION" != L.IN_VERSION
          )
        OR (
            X."VERSION"    IS NULL
          AND L.IN_VERSION IS NOT NULL
          )
        OR (

            X."VERSION"    IS NOT NULL
          AND L.IN_VERSION IS NULL
          )
        THEN X.VERSION_I + 1
        ELSE X.VERSION_I
      END,
      COMPARABLE_PROT_UPI,
      IN_LOAD_RELEASE CREATED,
      IN_LOAD_RELEASE as last,
      'N' DELETED,
      IN_TAXID
    FROM
      RNACEN.XREF X,

      LOAD_RETRO_TMP L
    WHERE
      X.AC        = L.IN_AC
    AND X.DBID    = L.IN_DBID
    AND X.DELETED = 'Y'
    AND x.last    = p_previous_release
    AND l.COMPARABLE_PROT_UPI IS NOT NULL
    AND NOT(
        X.UPI = L.COMPARABLE_PROT_UPI
      AND (
          X."VERSION" = L.IN_VERSION

        OR (
            X."VERSION"    IS NULL
          AND L.IN_VERSION IS NULL
          )
        )
      );

    --COMMIT;

  END;

$function$
