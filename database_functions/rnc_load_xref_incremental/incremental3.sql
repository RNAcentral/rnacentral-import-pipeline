CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.incremental3()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    -- inserts a new xref
    -- for dbid, ac where xref:
    -- a) does not exists at all (VERSION_I = 1)
    -- b) does not exists because has been deleted
    --    1) and was having the same UPI, Version (VERSION_I = OLD_VERSION_I)
    --    2) and was having a dferent UPI, VERSION (VERSION_I = OLD_VERSION_I
    -- + 1)
    INSERT
      /*+ PARALLEL (XREF 4) */

    INTO
      RNACEN.xref(
        ac,
        dbid,
        version,
        version_i,
        upi,
        created,
        last,
        deleted,
        taxid
      )
    SELECT
      /*+ PARALLEL (T 4) */
      T.IN_AC,
      T.IN_DBID,
      T.IN_VERSION,
      CASE T.COMPARABLE_PROT_UPI
        WHEN (
            SELECT DISTINCT (PREVIOUS_XREF.UPI)
            FROM
              RNACEN.XREF PREVIOUS_XREF
            WHERE
              PREVIOUS_XREF.AC          = T.IN_AC
            AND PREVIOUS_XREF.VERSION_I = T.MAX_PREVIOUS_XREF_VERSION_I
            AND PREVIOUS_XREF.DBID      = T.IN_DBID
          )
        THEN T.MAX_PREVIOUS_XREF_VERSION_I
        ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1
      END,
      T.COMPARABLE_PROT_UPI UPI,
      T.IN_LOAD_RELEASE CREATED,
      T.IN_LOAD_RELEASE as LAST,
      'N' DELETED,
      T.IN_TAXID TAXID

    FROM
      (
        SELECT
          /*+ PARALLEL (L 4) */
          L.IN_AC,
          L.IN_DBID,
          L.IN_VERSION,
          NVL (
          (
            SELECT
              /*+ PARALLEL (PREVIOUS_XREF 4) */
              MAX(PREVIOUS_XREF.VERSION_I)
            FROM

              RNACEN.XREF PREVIOUS_XREF
            WHERE
              PREVIOUS_XREF.AC     = L.IN_AC
            AND PREVIOUS_XREF.DBID = L.IN_DBID
          )
          , 0 ) MAX_PREVIOUS_XREF_VERSION_I,
          L.COMPARABLE_PROT_UPI,
          L.IN_LOAD_RELEASE,
          L.IN_TAXID
        FROM
          LOAD_RETRO_TMP L
        WHERE
          L.COMPARABLE_PROT_UPI IS NOT NULL

        AND NOT EXISTS (
            SELECT
              /*+ PARALLEL (X 4) */
              1
            FROM
              RNACEN.XREF X
            WHERE
              X.AC        = L.IN_AC
            AND X.DBID    = L.IN_DBID
            AND X.DELETED = 'N'
            AND X.LAST   <= L.IN_LOAD_RELEASE
              -- using also = in X.LAST <= L.IN_LOAD_RELEASE

              -- to eliminate from the selection
              -- the xref inserted right above after updating
              --  active xref exists with another protein/version
              -- MUST have been updated/inserted above
              -- nothing else to insert here
          )
      ) T;
    --COMMIT;

  END;

$function$

