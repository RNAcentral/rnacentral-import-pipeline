CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental1(p_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
BEGIN
    /*
    -- updates last, deleted and taxid
    -- of existing xrefs
    */
    UPDATE
      /*+ PARALLEL (U 4) */
      RNACEN.XREF U
    SET(
        U.LAST,
        U.DELETED,
        u.taxid

      )
      =
      (
        SELECT
          /*+ PARALLEL (L 4) */
          CASE
            WHEN(
                U.URS = L.COMPARABLE_PROT_UPI
              AND (
                  U."VERSION" = L.IN_VERSION
                OR (
                    U."VERSION"    = NULL
                  AND L.IN_VERSION = NULL
                  )
                )
              )
            THEN L.IN_LOAD_RELEASE
            ELSE p_previous_release
          END,
          CASE
            WHEN(
                U.URS = L.COMPARABLE_PROT_UPI

              AND (
                  U."VERSION" = L.IN_VERSION
                OR (
                    U."VERSION"    = NULL
                  AND L.IN_VERSION = NULL
                  )
                )
              )
            THEN 'N'
            ELSE 'Y'
          END,

          CASE
            WHEN(
                U.URS = L.COMPARABLE_PROT_UPI
              AND (
                  U."VERSION" = L.IN_VERSION
                OR (
                    U."VERSION"    = NULL
                  AND L.IN_VERSION = NULL
                  )
                )

              )
            THEN coalesce(L.IN_TAXID, U.TAXID)
            ELSE U.TAXID
          END
        FROM
          load_retro_tmp l
        WHERE
          u.deleted                = 'N'
        AND u.AC                   = l.IN_AC
        AND U.DBID                 = L.IN_DBID
        AND U.last                 < L.IN_LOAD_RELEASE
        AND l.COMPARABLE_PROT_UPI IS NOT NULL
      )

    WHERE
      U.DELETED = 'N'
    AND EXISTS (
        SELECT
          /*+ PARALLEL (X 4) */
          1
        FROM
          LOAD_RETRO_TMP x
        WHERE
          u.AC                     = x.IN_AC
        AND U.DBID                 = X.IN_DBID
        AND U.last                 < X.IN_LOAD_RELEASE

        AND x.COMPARABLE_PROT_UPI IS NOT NULL
      );

    --COMMIT;

  END;

$function$
