create or replace
PACKAGE BODY RNC_LOAD_XREF AS

  /*
  * Package for importing cross-references from the staging table.
  * PEL - Partition Exchange Loading
  * Two temporary tables are populated and inserted into the main XREF table
  * to optimize performance.
  */

  /*
    Truncate temporary tables, which requires disabling constraints.
  */
  PROCEDURE prepare_pel_tables IS
  BEGIN

    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_DELETED DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_NOT_DELETED DROP STORAGE';
    BEGIN
      -- INCAPSULATED TO IGNORE ANY 'DOES NOT EXISTS ERROR'
      -- WHEN RESTARTING THE PROCEDURE AFTER UNPREDICTABLE EXCEPTIONS
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$AC
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$AC"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$CREATED
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$CREATED"';
      --  DDL for Index XREF_PEL_DELETED$UPI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$UPI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$TAXID
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$TAXID"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$AC$UPPER
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER"';
      --------------------------------------------------------
      --  Constraints for Table XREF
      --------------------------------------------------------
      /*
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_CK1";
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DBID" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("CREATED" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("LAST" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("UPI" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("VERSION_I" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DELETED" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("TIMESTAMP" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("USERSTAMP" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("AC" NOT NULL DISABLE NOVALIDATE);
      */
      --------------------------------------------------------
      --  Ref Constraints for Table XREF
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK1"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK2"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK3"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK4"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$AC
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$CREATED
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED"';
      --  DDL for Index XREF_PEL_NOT_DELETED$UPI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$TAXID
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$AC$UPPER
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER"';
      --------------------------------------------------------
      --  Constraints for Table XREF
      --------------------------------------------------------
      /*
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_CK1";
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DBID" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("CREATED" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("LAST" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("UPI" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("VERSION_I" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DELETED" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("TIMESTAMP" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("USERSTAMP" NOT NULL DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("AC" NOT NULL DISABLE NOVALIDATE);
      */
      --------------------------------------------------------
      --  Ref Constraints for Table XREF
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK1"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK2"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK3"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK4"';

    EXCEPTION
    WHEN OTHERS THEN
      NULL;
    END;

  END prepare_pel_tables;

  /*
    Select existing xrefs with unchanged UPI and external db VERSION
    OR retire old entries (updated UPI / changed VERSION).
    New entres are inserted in populate_pel_tables2.
  */
  PROCEDURE populate_pel_tables1 (
    v_previous_release RNACEN.rnc_release.id%TYPE)
  IS
  BEGIN

    INSERT --+ APPEND PARALLEL (XREF_PEL_DELETED 2) PARALLEL (
      -- XREF_PEL_NOT_DELETED 2)
      ALL
    WHEN DELETED = 'Y' THEN
    INTO
      XREF_PEL_DELETED
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
      VALUES
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
    WHEN DELETED = 'N' THEN
    INTO
      XREF_PEL_NOT_DELETED
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
      VALUES
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
    SELECT
      /*+ PARALLEL (X 2) PARALLEL (L 2) */
      x.DBID,
      x.CREATED,
      x.UPI,
      x.VERSION_I,
      x."TIMESTAMP",
      x.USERSTAMP,
      x.AC,
      x."VERSION",
      -- x.LAST
      CASE
        WHEN
          (
            x.LAST  < l.in_load_release
          AND X.UPI = L.COMPARABLE_PROT_UPI
          AND
            (
              X."VERSION" = L.IN_VERSION
            OR
              (
                X."VERSION"    IS NULL
              AND L.IN_VERSION IS NULL
              )
            )
          )
        THEN L.IN_LOAD_RELEASE
        ELSE NVL (v_previous_release, x.LAST)
      END "LAST",
      -- x.DELETED
      CASE
        WHEN
          (
            x.LAST  < l.in_load_release
          AND X.UPI = L.COMPARABLE_PROT_UPI
          AND
            (
              X."VERSION" = L.IN_VERSION
            OR
              (
                X."VERSION"    IS NULL
              AND L.IN_VERSION IS NULL
              )
            )
          )
        THEN 'N'
        ELSE 'Y'
      END "DELETED",
      -- x.taxid
      CASE
        WHEN
          (
            x.LAST  < l.in_load_release
          AND X.UPI = L.COMPARABLE_PROT_UPI
          AND
            (
              X."VERSION" = L.IN_VERSION
            OR
              (
                X."VERSION"    IS NULL
              AND L.IN_VERSION IS NULL
              )
            )
          )
        THEN NVL (L.IN_TAXID, x.TAXID)
        ELSE x.TAXID
      END "TAXID"
    FROM
      LOAD_RETRO_TMP l,
      XREF x
    WHERE
        x.AC                   = l.IN_AC
    AND x.DBID                 = l.IN_DBID
    AND l.comparable_prot_upi IS NOT NULL
    AND x.deleted              = 'N';

    COMMIT;

  END populate_pel_tables1;

  /*
    Insert new xrefs, which supersede xrefs retired in populate_pel_tables1.
    Increment VERSION_I when the accession and dbid from an existing xref
    are matched to a different UPI.
  */
  PROCEDURE populate_pel_tables2 IS
  BEGIN
    INSERT --+ APPEND PARALLEL (XREF_PEL_NOT_DELETED 8)
    INTO
      XREF_PEL_NOT_DELETED
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
    SELECT
      /*+ PARALLEL (X 4) PARALLEL (L 4) */
      l.IN_DBID,
      l.IN_LOAD_RELEASE CREATED,
      l.COMPARABLE_PROT_UPI UPI,
      CASE
        WHEN
          (
            x.UPI != l.COMPARABLE_PROT_UPI
          )
        THEN x.VERSION_I + 1
        ELSE X.VERSION_I
      END VERSION_I,
      SYSDATE "TIMESTAMP",
      USER USERSTAMP,
      L.IN_AC AC,
      l.in_version "VERSION",
      l.IN_LOAD_RELEASE "LAST",
      'N' DELETED,
      IN_TAXID TAXID
    FROM
      LOAD_RETRO_TMP l,
      XREF x
    WHERE
        x.AC                   = L.IN_AC
    AND x.DBID                 = L.IN_DBID
    AND l.COMPARABLE_PROT_UPI IS NOT NULL
    AND NOT -- this condition differentiates this procedure from populate_pel_tables1
      (
        x.LAST  < l.IN_LOAD_RELEASE
      AND x.UPI = l.COMPARABLE_PROT_UPI
      AND
        (
          X."VERSION" = L.IN_VERSION
        OR
          (
            X."VERSION"    IS NULL
          AND L.IN_VERSION IS NULL
          )
        )
      )
    AND x.deleted = 'N';

    COMMIT;

  END populate_pel_tables2;

  /*
    -- Insert a new xref
    -- for dbid, ac where xref:
    -- a) does not exist at all (VERSION_I = 0 will become 1 when inserting
    --    notice outer join to previous_xref
    -- b) does not exist because has been deleted
    --    1) and had the same UPI, Version (VERSION_I = OLD_VERSION_I)
    --    2) and had a different UPI, VERSION (VERSION_I = OLD_VERSION_I + 1)
  */

  /*
    Create a temporary table with maximum internal versions for each
    ac/db/upi combination.
  */
  PROCEDURE load_upi_max_versions_table(
    p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE
  ) IS
  BEGIN

    EXECUTE IMMEDIATE 'TRUNCATE TABLE load_upi_max_versions DROP STORAGE';

    INSERT /*+ APPEND */
    INTO load_upi_max_versions
      SELECT DISTINCT
        /*+ PARALLEL (PREVIOUS_XREF 4) PARALLEL (L 4) */
        L.IN_AC,
        L.IN_DBID,
        MAX(NVL (PREVIOUS_XREF.VERSION_I, 0)) MAX_VERSION_I, -- VERSION_I set to 0 at first
        PREVIOUS_XREF.UPI
      FROM
        LOAD_RETRO_TMP L,
        RNACEN.XREF PREVIOUS_XREF
      WHERE
        L.COMPARABLE_PROT_UPI IS NOT NULL
        AND L.IN_DBID = p_in_dbid
        AND PREVIOUS_XREF.DBID (+)  = L.IN_DBID -- outer join, left side can be NULL
        AND PREVIOUS_XREF.AC   (+)  = L.IN_AC   -- outer join, left side can be NULL
        AND NOT EXISTS
          (
            SELECT
              /*+ PARALLEL (X 4) */
              1
            FROM
              RNACEN.XREF X
            WHERE
                X.AC      = L.IN_AC
            AND X.DBID    = L.IN_DBID
            AND X.DELETED = 'N'
          )
      GROUP BY
        L.IN_AC,
        L.IN_DBID,
        PREVIOUS_XREF.UPI;

    COMMIT;

    BEGIN
      DBMS_STATS.GATHER_TABLE_STATS ( OWNNAME => 'RNACEN', TABNAME => 'load_upi_max_versions', ESTIMATE_PERCENT => 10 );
    END;

  END load_upi_max_versions_table;

  /*
    Create a temporary table with maximum internal versions for each
    ac/db combination.
  */
  PROCEDURE load_max_versions_table IS
  BEGIN
    EXECUTE IMMEDIATE 'TRUNCATE TABLE load_max_versions DROP STORAGE';

    INSERT /*+ APPEND */ INTO load_max_versions
    SELECT DISTINCT
      /*+ PARALLEL (L 4) */
      l.AC,
      L.DBID,
      MAX(L.MAX_VERSION_I) MAX_VERSION_I
    FROM load_upi_max_versions l
    GROUP BY
      L.AC,
      L.DBID;

    COMMIT;

    BEGIN
      DBMS_STATS.GATHER_TABLE_STATS ( OWNNAME => 'RNACEN', TABNAME => 'load_max_versions', ESTIMATE_PERCENT => 10 );
    END;

  END load_max_versions_table;

  /*

  */
  PROCEDURE populate_pel_tables3 (
    p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE
  ) IS
  BEGIN

    INSERT --+ APPEND PARALLEL (XREF_PEL_NOT_DELETED 4)
    INTO
      RNACEN.XREF_PEL_NOT_DELETED
      (
        AC,
        DBID,
        "VERSION",
        VERSION_I,
        UPI,
        CREATED,
        "LAST",
        "DELETED",
        TAXID,
        "TIMESTAMP",
        USERSTAMP
      )
    SELECT
      /*+ PARALLEL (T 4) PARALLEL (LUMV 4) */
      T.IN_AC,
      T.IN_DBID,
      T.IN_VERSION,
      CASE
        WHEN T.MAX_PREVIOUS_XREF_VERSION_I = 0 -- new xref
          THEN 1
        WHEN LUMV.UPI = T.COMPARABLE_PROT_UPI  -- same UPI, keep old VERSION_I
          THEN T.MAX_PREVIOUS_XREF_VERSION_I
        ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1 -- updated UPI, update VERSION_I
      END VERSION_I,
      T.COMPARABLE_PROT_UPI UPI,
      T.IN_LOAD_RELEASE CREATED,
      T.IN_LOAD_RELEASE LAST,
      'N' DELETED,
      T.IN_TAXID,
      SYSDATE "TIMESTAMP",
      USER USERSTAMP
    FROM
      LOAD_UPI_MAX_VERSIONS lumv,
      (
        SELECT
          /*+ PARALLEL (L 4) PARALLEL (LMV 4) */
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
      ) T
    WHERE
        LUMV.AC (+)            = T.IN_AC
    AND LUMV.MAX_VERSION_I (+) = T.MAX_PREVIOUS_XREF_VERSION_I
    AND LUMV.DBID (+)          = T.IN_DBID;

    COMMIT;

  END populate_pel_tables3;

  /*
    Existing xrefs not present in the new release for which the old entries
    must be retired or those that were already retired and weren't updated.
  */
  PROCEDURE populate_pel_tables4 (
    p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE,
    v_previous_release IN RNACEN.rnc_release.ID%TYPE
  ) IS
  BEGIN
    INSERT --+ APPEND PARALLEL (XREF_PEL_DELETED 4)
    INTO
      XREF_PEL_DELETED
      (
        DBID,
        CREATED,
        UPI,
        VERSION_I,
        "TIMESTAMP",
        USERSTAMP,
        AC,
        "VERSION",
        "LAST",
        "DELETED",
        "TAXID"
      )
    SELECT
      /*+ PARALLEL (X 4) */
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
          NVL (v_previous_release, X.LAST) -- retire updated entries
        ELSE                               -- retired xref that weren't updated
          X.LAST
      END LAST,
      'Y' DELETED,
      x.taxid
    FROM
      RNACEN.XREF X
    WHERE
      X.dbid = p_in_dbid
    AND
      (
        x.deleted = 'Y'
      OR NOT EXISTS
        (
          SELECT
            /*+ PARALLEL (L 4) */
            1
          FROM
            LOAD_RETRO_TMP L
          WHERE
            L.COMPARABLE_PROT_UPI IS NOT NULL
          AND X.AC                 = L.IN_AC
          AND X.DBID               = L.IN_DBID
        )
      );

    COMMIT;

  END populate_pel_tables4;

  /*
    Prepare temporary PEL tables and perform the exchange.
  */
  PROCEDURE do_pel_exchange(
    p_in_db_id IN RNACEN.rnc_database.ID%TYPE)
  IS
  BEGIN
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_DELETED$AC
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$AC" ON "RNACEN"."XREF_PEL_DELETED" ("AC") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$AC" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_DELETED$CREATED
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$CREATED" ON "RNACEN"."XREF_PEL_DELETED" ("CREATED") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$CREATED" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_DELETED$UPI
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$UPI" ON "RNACEN"."XREF_PEL_DELETED" ("UPI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$UPI" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_DELETED$TAXID
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$TAXID" ON "RNACEN"."XREF_PEL_DELETED" ("TAXID") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$TAXID" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_DELETED$AC$UPPER
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER" ON "RNACEN"."XREF_PEL_DELETED" (UPPER("AC")) PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER" NOPARALLEL';
    --------------------------------------------------------
    --  Constraints for Table XREF
    --------------------------------------------------------
    /*
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_CK1" CHECK (deleted IN ('Y', 'N')) ENABLE VALIDATE;
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DBID" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("CREATED" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("LAST" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("UPI" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("VERSION_I" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DELETED" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("TIMESTAMP" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("USERSTAMP" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("AC" NOT NULL ENABLE VALIDATE);
    */
    --------------------------------------------------------
    --  Ref Constraints for Table XREF
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK1" FOREIGN KEY ("CREATED") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK2" FOREIGN KEY ("DBID") REFERENCES "RNACEN"."RNC_DATABASE" ("ID") ENABLE VALIDATE';

    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK3" FOREIGN KEY ("LAST") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK4" FOREIGN KEY ("UPI") REFERENCES "RNACEN"."RNA" ("UPI") ENABLE VALIDATE';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_NOT_DELETED$AC
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("AC") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_NOT_DELETED$CREATED
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("CREATED") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_NOT_DELETED$UPI
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("UPI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_NOT_DELETED$TAXID
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("TAXID") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID" NOPARALLEL';
    --------------------------------------------------------
    --  DDL for Index XREF_PEL_NOT_DELETED$AC$UPPER
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER" ON "RNACEN"."XREF_PEL_NOT_DELETED" (UPPER("AC")) PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
    EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER" NOPARALLEL';
    --------------------------------------------------------
    --  Constraints for Table XREF
    --------------------------------------------------------
    /*
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_CK1" CHECK (deleted IN ('Y', 'N')) ENABLE VALIDATE;
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DBID" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("CREATED" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("LAST" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("UPI" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("VERSION_I" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DELETED" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("TIMESTAMP" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("USERSTAMP" NOT NULL ENABLE VALIDATE);
    ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("AC" NOT NULL ENABLE VALIDATE);
    */
    --------------------------------------------------------
    --  Ref Constraints for Table XREF
    --------------------------------------------------------
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK1" FOREIGN KEY ("CREATED") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK2" FOREIGN KEY ("DBID") REFERENCES "RNACEN"."RNC_DATABASE" ("ID") ENABLE VALIDATE';
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK3" FOREIGN KEY ("LAST") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
    EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK4" FOREIGN KEY ("UPI") REFERENCES "RNACEN"."RNA" ("UPI") ENABLE VALIDATE';

    DBMS_STATS.GATHER_TABLE_STATS (OWNNAME => 'RNACEN', TABNAME =>
      'XREF_PEL_DELETED', ESTIMATE_PERCENT => DBMS_STATS.AUTO_SAMPLE_SIZE,
      METHOD_OPT => 'FOR ALL COLUMNS SIZE AUTO', DEGREE => 8, CASCADE => TRUE);
    DBMS_STATS.GATHER_TABLE_STATS (OWNNAME => 'RNACEN', TABNAME =>
      'XREF_PEL_NOT_DELETED', ESTIMATE_PERCENT => DBMS_STATS.AUTO_SAMPLE_SIZE,
      METHOD_OPT => 'FOR ALL COLUMNS SIZE AUTO', DEGREE => 4, CASCADE => TRUE);

    -- perform partition exchange, main PEL statements
    EXECUTE IMMEDIATE 'ALTER TABLE XREF EXCHANGE SUBPARTITION XREF_P' || TO_CHAR (p_in_db_id ) || '_NOT_DELETED WITH TABLE XREF_PEL_NOT_DELETED INCLUDING INDEXES WITHOUT VALIDATION';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF EXCHANGE SUBPARTITION XREF_P' || TO_CHAR ( p_in_db_id )|| '_DELETED WITH TABLE XREF_PEL_DELETED INCLUDING INDEXES WITHOUT VALIDATION';

  END do_pel_exchange;

  /*
    Main entry point.
  */
  PROCEDURE load_xref(
    p_previous_release IN RNACEN.rnc_release.ID%TYPE,
    p_in_dbid IN RNACEN.rnc_database.id%TYPE
  )
  IS
  BEGIN
    prepare_pel_tables;
    populate_pel_tables1(p_previous_release);
    populate_pel_tables2;

    load_upi_max_versions_table(p_in_dbid);
    load_max_versions_table;
    populate_pel_tables3(p_in_dbid);
    populate_pel_tables4(p_in_dbid, p_previous_release);

    do_pel_exchange(p_in_dbid);
  END;

END RNC_LOAD_XREF;