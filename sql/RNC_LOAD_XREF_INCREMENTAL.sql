set define off

create or replace PACKAGE RNC_LOAD_XREF_INCREMENTAL AS

  PROCEDURE load_xref_incremental(p_previous_release IN RNACEN.rnc_release.ID%TYPE,
                                  p_in_dbid          IN RNACEN.rnc_database.id%TYPE);

END RNC_LOAD_XREF_INCREMENTAL;
/
create or replace PACKAGE BODY RNC_LOAD_XREF_INCREMENTAL AS


  PROCEDURE incremental1(
    p_previous_release IN RNACEN.rnc_release.ID%TYPE)
  IS

  BEGIN
    /*
    -- updates last, deleted and taxid
    -- of existing xrefs
    */
    UPDATE
      /*+ PARALLEL (U 4) */
      RNACEN.XREF U
    SET
      (
        U.LAST,
        U.DELETED,
        u.taxid

      )
      =
      (
        SELECT
          /*+ PARALLEL (L 4) */
          CASE
            WHEN
              (
                U.UPI = L.COMPARABLE_PROT_UPI
              AND
                (
                  U."VERSION" = L.IN_VERSION
                OR

                  (
                    U."VERSION"    IS NULL
                  AND L.IN_VERSION IS NULL
                  )
                )
              )
            THEN L.IN_LOAD_RELEASE
            ELSE p_previous_release
          END "LAST",
          CASE
            WHEN
              (
                U.UPI = L.COMPARABLE_PROT_UPI

              AND
                (
                  U."VERSION" = L.IN_VERSION
                OR
                  (
                    U."VERSION"    IS NULL
                  AND L.IN_VERSION IS NULL
                  )
                )
              )
            THEN 'N'
            ELSE 'Y'
          END "DELETED",

          CASE
            WHEN
              (
                U.UPI = L.COMPARABLE_PROT_UPI
              AND
                (
                  U."VERSION" = L.IN_VERSION
                OR
                  (
                    U."VERSION"    IS NULL
                  AND L.IN_VERSION IS NULL
                  )
                )

              )
            THEN NVL (L.IN_TAXID, U.TAXID)
            ELSE U.TAXID
          END "TAXID"
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
    AND EXISTS
      (
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

    COMMIT;

  end incremental1;


  PROCEDURE incremental2(
    p_previous_release IN RNACEN.rnc_release.ID%TYPE)
  IS
  BEGIN


    -- inserts new xref
    -- for each deleted = 'Y'
    -- in the above update
    INSERT
      /*+ PARALLEL (XREF 4) */
    INTO
      RNACEN.xref
      (
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
      /*+ PARALLEL (L 4) PARALLEL (X 4) */
      l.IN_AC,
      l.IN_DBID,
      l.IN_VERSION,
      CASE
        WHEN
          (

            X.UPI != L.COMPARABLE_PROT_UPI
          )
        OR
          (
            X."VERSION" != L.IN_VERSION
          )
        OR
          (
            X."VERSION"    IS NULL
          AND L.IN_VERSION IS NOT NULL
          )
        OR
          (

            X."VERSION"    IS NOT NULL
          AND L.IN_VERSION IS NULL
          )
        THEN X.VERSION_I + 1
        ELSE X.VERSION_I
      END VERSION_I_NEW,
      COMPARABLE_PROT_UPI,
      IN_LOAD_RELEASE CREATED,
      IN_LOAD_RELEASE LAST,
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
    AND NOT
      (
        X.UPI = L.COMPARABLE_PROT_UPI
      AND
        (
          X."VERSION" = L.IN_VERSION

        OR
          (
            X."VERSION"    IS NULL
          AND L.IN_VERSION IS NULL
          )
        )
      );

    COMMIT;

  END incremental2;



  PROCEDURE incremental3
  IS
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
      RNACEN.xref
      (
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
        WHEN
          (
            SELECT
              /*+ PARALLEL (PREVIOUS_XREF 4) */
              UNIQUE(PREVIOUS_XREF.UPI)
            FROM
              RNACEN.XREF PREVIOUS_XREF

            WHERE
              PREVIOUS_XREF.AC          = T.IN_AC
            AND PREVIOUS_XREF.VERSION_I = T.MAX_PREVIOUS_XREF_VERSION_I
            AND PREVIOUS_XREF.DBID      = T.IN_DBID
          )
        THEN T.MAX_PREVIOUS_XREF_VERSION_I
        ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1
      END VERSION_I,
      T.COMPARABLE_PROT_UPI UPI,
      T.IN_LOAD_RELEASE CREATED,
      T.IN_LOAD_RELEASE LAST,
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

        AND NOT EXISTS
          (
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
    COMMIT;

  END incremental3;


  /*

    Main entry point.
  */
  PROCEDURE load_xref_incremental(
    p_previous_release IN RNACEN.rnc_release.ID%TYPE,
    p_in_dbid IN RNACEN.rnc_database.id%TYPE
  )
  IS
  BEGIN

    /*
    -- this is an INCREMENTAL release
    -- this is NOT a FULL release
    -- better to AVOID Partition Exchange Loading

    -- updates last, deleted and taxid
    -- of existing xrefs
    */
    /*
    Unfortunatley when using parallel with conventional DML we are hitting a
    bug. I am disabling parall DML but leaving PARALLEL hints in place so that
    the query part of the DML statements should still be able to use it.
    Bug 9831227 : Parallel DML fails with ORA-12839 or ORA-600
    Description
     you have a parent table and child table with a parent referential
    constraint
    then running Parallel DML again this may fail with an unexpected ORA-12839
    or even fail with an ORA-600.

    Workaround
    Setting "_disable_parallel_conventional_load" = true can help avoid this.
    */
    EXECUTE IMMEDIATE 'ALTER SESSION DISABLE PARALLEL DML';

    incremental1(p_previous_release);

    incremental2(p_previous_release);

    incremental3();

    /*
    Re-enable parall DML previously disabled because of the

    Bug 9831227 : Parallel DML fails with ORA-12839 or ORA-600
    Description
     you have a parent table and child table with a parent referential
    constraint
    then running Parallel DML again this may fail with an unexpected ORA-12839
    or even fail with an ORA-600.
    Workaround
    Setting "_disable_parallel_conventional_load" = true can help avoid this.
    */
    EXECUTE IMMEDIATE 'ALTER SESSION FORCE PARALLEL DML';

  END;


END RNC_LOAD_XREF_INCREMENTAL;
/
set define on
